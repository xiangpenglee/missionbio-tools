# -*- coding: utf-8 -*-
"""
mission bio single-cell pipeline code
written by ben 11.25.2018

"""

# modules
from __future__ import division
import os
import os.path
import csv
from itertools import product, combinations, izip
import json
import subprocess

class TapestriSample(object):
    # class for storing metadata for each tapestri sample (one tube, run, etc...)

    def __init__(self,
                 sample_num,
                 r1,
                 r2,
                 output_folder):

        self.sample_num = sample_num            # number identifying sample (or tube)

        self.r1 = r1                      # R1 fastq path
        self.r2 = r2                      # R2 fastq path

        self.r1_trimmed = output_folder + '/' + str(sample_num) + '-trimmed_R1.fastq'     # trimmed R1 file
        self.r2_trimmed = output_folder + '/' + str(sample_num) + '-trimmed_R2.fastq'     # trimmed R2 file

        self.cell_cutadapt = output_folder + '/' + str(sample_num) + '_cell_barcode_cutadapt.txt'     # cutadapt report (cell barcodes)
        self.barcode_counts = output_folder + '/' + str(sample_num) + '_barcode_counts.tsv'           # sample barcode counts
        self.read_id_json = output_folder + '/' + str(sample_num) + '_barcodes.json'                  # json of read ids

    def filter_valid_reads(self,
                           r1_start,
                           mb_barcodes,
                           bar_ind_1,
                           bar_ind_2):
        # filter r1 files to only keep reads with correct barcode structure in r1

        r1_in = self.r1

        cmd = 'python3 /usr/local/bin/cutadapt' \
              ' -a r1_start=%s' \
              ' -j 16 -O 8 -e 0.2 %s' \
              % (r1_start,
                 r1_in)

        trim_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        bar_count = 0       # total count of all barcodes
        barcodes = {}       # counts for each barcode
        read_id_dict = {}   # dict for storing read ids

        # iterate through all Read 1 records
        for line in trim_process.stdout:

            # R1
            header_1 = line.strip()                             # header
            id_1 = header_1.split(' ')[0][1:]                   # read id
            seq_1 = trim_process.stdout.next().strip()          # sequence string
            trim_process.stdout.next()                          # +
            trim_process.stdout.next()                          # qual

            # find barcodes and check that they are a valid MB barcode
            check = check_seq(seq_1, bar_ind_1, bar_ind_2, mb_barcodes)

            if check == 'fail':
                continue

            else:
                barcode = check[0] + check[1] + '-' + str(self.sample_num)

                # add barcode to dict
                try:
                    barcodes[barcode] += 1

                except KeyError:
                    barcodes[barcode] = 1

                read_id_dict[id_1] = barcode
                bar_count += 1

                # print counter
                if bar_count % 1e6 == 0:
                    print '%d valid trimmed pairs observed.' % bar_count

        # save barcode list to file
        barcodes_tsv = open(self.barcode_counts, 'w')
        for bar in barcodes:
            barcodes_tsv.write('%s\t%d\n' % (bar, barcodes[bar]))
        barcodes_tsv.close()

        print 'Exporting barcodes to JSON...'

        # export barcodes to json file
        json_export(read_id_dict, self.read_id_json)

        print 'Barcode correction complete...'

    def barcode_reads(self,
                      r1_start,
                      r1_end,
                      r2_end,
                      r1_min_len,
                      r2_min_len,
                      paired_end):
        # for valid reads, add barcode header to fastq file and trim

        r1_in = self.r1
        r2_in = self.r2

        r1_out = open(self.r1_trimmed, 'w')
        r2_out = open(self.r2_trimmed, 'w')

        print 'Importing JSON...'

        read_id_dict = json_import(self.read_id_json)

        print 'Beginning adapter cutting...'

        if paired_end == 'Paired':

            cmd = 'python3 /usr/local/bin/cutadapt' \
                  ' -g %s' \
                  ' -a %s' \
                  ' -A %s' \
                  ' --interleaved -j 16 -n 3 -O 8 -e 0.2 %s %s' \
                  ' 2> %s' \
                  % (r1_start,
                     r1_end,
                     r2_end,
                     r1_in,
                     r2_in,
                     self.cell_cutadapt)

        else:

            cmd = 'python3 /usr/local/bin/cutadapt' \
                  ' -g %s' \
                  ' -a %s' \
                  ' -j 16 -n 2 -O 8 -e 0.2 %s' \
                  ' 2> %s' \
                  % (r1_start,
                     r1_end,
                     r1_in,
                     self.cell_cutadapt)

        trim_process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)

        bar_count = 0  # total count of all barcodes

        # iterate through all Read 1 records
        for line in trim_process.stdout:

            if paired_end == 'Paired':

                # R1
                header_1 = line.strip()
                id_1 = header_1.split(' ')[0][1:]
                seq_1 = trim_process.stdout.next().strip()
                trim_process.stdout.next()
                qual_1 = trim_process.stdout.next().strip()

                # R2
                header_2 = trim_process.stdout.next().strip()
                id_2 = header_1.split(' ')[0][1:]
                seq_2 = trim_process.stdout.next().strip()
                trim_process.stdout.next()
                qual_2 = trim_process.stdout.next().strip()

                assert id_1 == id_2, 'Read IDs do not match! Check input FASTQ files.'

                try:
                    cell_barcode = read_id_dict[id_1]
                except KeyError:
                    continue

                if len(seq_1) < r1_min_len or len(seq_2) < r2_min_len:
                    continue

                else:
                    # add barcode to header
                    id = '@' + id_1 + '_' + cell_barcode
                    header_1 = id# + ' ' + header_1.split(' ')[1]
                    header_2 = id# + ' ' + header_2.split(' ')[1]

                    # write to output fastq files
                    r1_out.write('%s\n%s\n+\n%s\n' % (header_1, seq_1, qual_1))
                    r2_out.write('%s\n%s\n+\n%s\n' % (header_2, seq_2, qual_2))

                    bar_count += 1

                    # print counter
                    if bar_count % 1e6 == 0:
                        print '%d valid trimmed pairs observed.' % bar_count

            else:

                # R1
                header_1 = line.strip()
                id_1 = header_1.split(' ')[0][1:]
                seq_1 = trim_process.stdout.next().strip()
                trim_process.stdout.next()
                qual_1 = trim_process.stdout.next().strip()

                try:
                    cell_barcode = read_id_dict[id_1]
                except KeyError:
                    continue

                if len(seq_1) < r1_min_len:
                    continue

                else:
                    # add barcode to header
                    id = '@' + id_1 + '_' + cell_barcode
                    header_1 = id# + ' ' + header_1.split(' ')[1]

                    # write to output fastq files
                    r1_out.write('%s\n%s\n+\n%s\n' % (header_1, seq_1, qual_1))

                    bar_count += 1

                    # print counter
                    if bar_count % 1e6 == 0:
                        print '%d valid trimmed pairs observed.' % bar_count

        print '%d total valid trimmed pairs saved to file.' % bar_count

        r1_out.close()
        r2_out.close()

def json_import(filename):
    # imports json data into python. If json file does not exist, returns empty {}

    if not os.path.exists(filename):
        json_obj = {}
    else:
        with open(filename) as f:
            json_obj = json.load(f)

    return json_obj

def json_export(json_obj, filename, overwrite=True, update=False):
    # exports a json object to file. If overwrite is off, writing will fail. Can also update an existing json

    if os.path.exists(filename) and not overwrite:
        print 'File exists. Will not overwrite. Exiting...'
        raise SystemExit

    elif update:
        if not os.path.exists(filename):
            old_json = {}
        else:
            old_json = json_import(filename)

        json_obj.update(old_json)

    with open(filename, 'w') as out:
        json.dump(json_obj, out)

def load_barcodes(barcode_file, max_dist, check_barcodes=False):
    # loads barcodes from csv and checks all pairwise distances to ensure error correction will work
    # returns a dictionary of barcodes with their descriptions

    # load barcodes from csv
    reader = csv.reader(open(barcode_file, 'r'))
    barcodes = {}
    for barcode, desc in reader:
        barcodes[barcode] = desc

    # optional distance check for barcodes (can turn off once validated)
    if check_barcodes:

        # check all barcodes have same length
        lengths = map(len, barcodes.keys())
        if len(set(lengths)) != 1:
            print 'Barcodes must all be same length! Exiting...'
            raise SystemExit

        # check pairwise hamming distances
        dist_req = 2 * max_dist + 1         # for this max_dist, need this distance between all barcodes
        pairs = list(combinations(barcodes.keys(), 2))
        for pair in pairs:
            if hd(pair[0], pair[1]) < dist_req:
                print 'Error: The edit distance between barcodes %s and %s is less than %d.\n' \
                      'An error correction of %d bases will not work.' % (pair[0], pair[1], dist_req, max_dist)

    return barcodes

def generate_hamming_dict(barcodes):
    # return dictionary of strings within 1 hamming distance of barcodes

    # key: uncorrected barcode; value: corrected barcode
    barcode_dict = {}

    for barcode in barcodes:
        barcode_dict[barcode] = barcode
        hd_1 = sorted(hamming_circle(barcode, 1))
        for uncorrected_bar in hd_1:
            barcode_dict[uncorrected_bar] = barcode

    return barcode_dict

# from https://codereview.stackexchange.com/a/88919
def hamming_circle(s, n):
    # generate strings over alphabet whose hamming distance from s is exactly n

    alphabet = 'ATCG'
    s = s.upper()

    for positions in combinations(range(len(s)), n):

        for replacements in product(range(len(alphabet) - 1), repeat=n):
            cousin = list(s)

            for p, r in zip(positions, replacements):
                if cousin[p] == alphabet[r]:
                    cousin[p] = alphabet[-1]

                else:
                    cousin[p] = alphabet[r]

            yield ''.join(cousin)

def hd(str1, str2):
    # compute Hamming distance between two strings of equal length
    assert len(str1) == len(str2)
    return sum(c1 != c2 for c1, c2 in izip(str1, str2))

def check_seq(seq, bar_ind_1, bar_ind_2, barcodes):
    # checks a sequence for valid barcodes
    # outputs a list of raw and corrected information if valid, 'fail' otherwise

    try:
        bar_1 = ''.join([seq[i] for i in bar_ind_1])
        bar_2 = ''.join([seq[i] for i in bar_ind_2])

    except IndexError:
        return 'fail'

    # try and error correct the barcode
    corr_barcode_1 = correct_barcode(barcodes, bar_1)
    corr_barcode_2 = correct_barcode(barcodes, bar_2)

    if (corr_barcode_1 == 'invalid') or (corr_barcode_2 == 'invalid'):
        return 'fail'

    return [corr_barcode_1, corr_barcode_2]

def correct_barcode(barcodes, raw_barcode):
    # attempts to correct a raw barcode using list of valid barcodes
    # return: raw barcode, corrected barcode

    # check if barcode is in hamming dictionary, if not, return invalid
    try:
        return barcodes[raw_barcode]
    except KeyError:
        return 'invalid'
