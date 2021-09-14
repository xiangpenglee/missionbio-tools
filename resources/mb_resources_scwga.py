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
from itertools import product, combinations, izip, izip_longest
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

    def barcode_reads(self,
                      r1_start,
                      r1_end,
                      r2_start,
                      r2_end,
                      r1_min_len,
                      r2_min_len,
                      bar_ind_1,
                      bar_ind_2,
                      mb_barcodes,
                      paired_end):
        # for valid reads, add barcode header to fastq file and trim

        # input files
        r1_in = self.r1
        r2_in = self.r2

        # output files
        r1_out = open(self.r1_trimmed, 'w')
        r2_out = open(self.r2_trimmed, 'w')

        # generate info file for reads (to identify valid barcodes)
        info_cmd = 'cutadapt -g %s -O 8 -e 0.2 %s -j 8 --info-file - -o /dev/null --quiet' % (r1_start, r1_in)
        info_file = subprocess.Popen(info_cmd, stdout=subprocess.PIPE, shell=True)

        # generate list of seqs
        if paired_end:
            trim_cmd = 'cutadapt -g %s -a %s -G %s -A %s -O 8 -e 0.2 %s %s -j 8 -n 3 --interleaved -o - 2> %s' % (r1_start,
                                                                                                                  r1_end,
                                                                                                                  r2_start,
                                                                                                                  r2_end,
                                                                                                                  r1_in,
                                                                                                                  r2_in,
                                                                                                                  self.cell_cutadapt)

        else:
            trim_cmd = 'cutadapt -g %s -a %s -O 8 -e 0.2 %s -j 8 -n 2 -o - 2> %s' % (r1_start,
                                                                                      r1_end,
                                                                                      r1_in,
                                                                                      self.cell_cutadapt)

        trim_file = subprocess.Popen(trim_cmd, stdout=subprocess.PIPE, shell=True)

        valid_reads = 0     # valid, barcoded read count
        total_reads = 0     # total count of all reads
        too_short = 0
        barcodes = {}       # counts for each barcode

        # iterate through info file (barcodes) and trim file (reads)
        for info_line, trim_line in izip(info_file.stdout, trim_file.stdout):

            assert info_line.split('\t')[0] == trim_line.strip()[1:], 'Cluster IDs do not match!'

            total_reads += 1

            # no valid adapter
            if info_line.split('\t')[1] == '-1':
                trim_file.stdout.next()
                trim_file.stdout.next()
                trim_file.stdout.next()

                if paired_end:
                    trim_file.stdout.next()
                    trim_file.stdout.next()
                    trim_file.stdout.next()
                    trim_file.stdout.next()

                continue

            # contains valid adapter
            else:
                bar_seq = info_line.split('\t')[4]

                # find barcodes and check that they are a valid MB barcode
                check = check_seq(bar_seq, bar_ind_1, bar_ind_2, mb_barcodes)

                # not a valid barcode
                if check == 'fail':
                    trim_file.stdout.next()
                    trim_file.stdout.next()
                    trim_file.stdout.next()

                    if paired_end:
                        trim_file.stdout.next()
                        trim_file.stdout.next()
                        trim_file.stdout.next()
                        trim_file.stdout.next()

                    continue

                # valid barcode
                else:
                    barcode = check[0] + check[1] + '-' + str(self.sample_num)

                    # R1 from trimmed file
                    header_1 = trim_line.strip()
                    id_1 = header_1.split(' ')[0][1:]
                    seq_1 = trim_file.stdout.next().strip()
                    trim_file.stdout.next()
                    qual_1 = trim_file.stdout.next().strip()

                    if paired_end:

                        # R2 from trimmed file
                        header_2 = trim_file.stdout.next().strip()
                        id_2 = header_1.split(' ')[0][1:]

                        assert id_1 == id_2, 'Cluster IDs in interleaved input do not match!'

                        seq_2 = trim_file.stdout.next().strip()
                        trim_file.stdout.next()
                        qual_2 = trim_file.stdout.next().strip()

                    # check reads for length
                    if len(seq_1) < r1_min_len:
                        too_short += 1
                        continue

                    if paired_end:
                        if len(seq_2) < r2_min_len:
                            too_short += 1
                            continue

                    # barcode is valid, add to dict
                    try:
                        barcodes[barcode] += 1

                    except KeyError:
                        barcodes[barcode] = 1

                    # add barcoded headers and reads to file
                    id = '@' + id_1 + '_' + barcode
                    header_1 = id# + ' ' + header_1.split(' ')[1]
                    r1_out.write('%s\n%s\n+\n%s\n' % (header_1, seq_1, qual_1))

                    if paired_end:
                        header_2 = id# + ' ' + header_2.split(' ')[1]
                        r2_out.write('%s\n%s\n+\n%s\n' % (header_2, seq_2, qual_2))

                    valid_reads += 1

                    # print counter
                    if valid_reads % 1e6 == 0:
                        print '%d valid trimmed pairs saved to file.' % valid_reads

        # save barcode list to file
        barcodes_tsv = open(self.barcode_counts, 'w')
        for bar in barcodes:
            barcodes_tsv.write('%s\t%d\n' % (bar, barcodes[bar]))
        barcodes_tsv.close()

    def count_barcodes(self,
                       r1_start,
                       bar_ind_1,
                       bar_ind_2,
                       mb_barcodes):
    # count barcodes from fastq file

        # for valid reads, add barcode header to fastq file and trim

        # input files
        r1_in = self.r1

        # generate info file for reads (to identify valid barcodes)
        info_cmd = 'cutadapt -g %s -O 8 -e 0.2 %s -j 16 --info-file - -o /dev/null --quiet' % (r1_start, r1_in)
        info_file = subprocess.Popen(info_cmd, stdout=subprocess.PIPE, shell=True)

        barcodes = {}  # counts for each barcode

        # iterate through info file (barcodes) and trim file (reads)
        for info_line in info_file.stdout:

            # no valid adapter
            if info_line.split('\t')[1] == '-1':
                continue

            # contains valid adapter
            else:
                bar_seq = info_line.split('\t')[4]

                # find barcodes and check that they are a valid MB barcode
                check = check_seq(bar_seq, bar_ind_1, bar_ind_2, mb_barcodes)

                # not a valid barcode
                if check == 'fail':
                    continue

                # valid barcode
                else:
                    barcode = check[0] + check[1] + '-' + str(self.sample_num)

                    # add barcode to dict
                    try:
                        barcodes[barcode] += 1

                    except KeyError:
                        barcodes[barcode] = 1

        # save barcode list to file
        barcodes_tsv = open(self.barcode_counts, 'w')
        for bar in barcodes:
            barcodes_tsv.write('%s\t%d\n' % (bar, barcodes[bar]))
        barcodes_tsv.close()

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
