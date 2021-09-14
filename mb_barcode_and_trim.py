'''
Mission Bio Barcode and Trim
Ben Demaree - 10.4.2019

A simple script for pre-processing of Mission Bio data.

'''

import os
import subprocess
import sys
import argparse
from multiprocessing import Process

# add mission bio cell processing script
sys.path.append(os.path.join(sys.path[0], 'resources'))
import mb_resources

def wait(processes):
    # waits for processes to finish
    return [process.communicate() for process in processes]

def file_summary(samples):
    # print summary of all samples identified

    for sample in samples:

        s = vars(sample)

        for item in s:
            print item,': ' ,s[item]
        print '\n'

    print '%d samples identified.\n' % len(samples)

def generate_samples(R1_files, R2_files, output_folder, sample_label):
    # store sample filenames in Sample objects

    samples = []        # list of sample objects

    for i in range(len(R1_files)):
        # assign samples to objects

        r1 = R1_files[i]
        r2 = R2_files[i]

        sample_num = i + 1
        sample_label = sample_label + '-' + str(sample_num)

        samples.append(mb_resources.TapestriSample(sample_num, r1, r2, output_folder))

    return samples


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''
    
    Mission Bio Barcode and Trim
    Ben Demaree - 10.4.2019
    
    A simple script for pre-processing of Mission Bio data.
    
    Performs barcode error correction and adds cell barcode to fastq headers in the following format:
    
    ---
    @M00179:262:000000000-G4F3T:1:1101:16518:1378_CTCTGCGGTCACTAGACT-1
    CCGACTCGAGCGACTGATGTGGTGAGAGCAAACACAGAGAGAATGGATAAGAAATACTGAATGGGCCGGGCGTGGTGGCTCACGCCTGTAATCCCAGCACTT
    +
    GHGGGGGGGGGGGGGGFHHH4D?F3GHHHHHHHGHGHHHHGAFHHHFHHHBHHGHHHHGFFHF30GHGGGGGGGGGEGHHHHHHHGHFHFEHHHHHGHHGHH
    ---

    ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('sample_label', type=str,
                        help='Label for this sample')

    parser.add_argument('input_fastq_folder', type=str,
                        help='Input folder containing raw FASTQ files (must have fastq.gz extension)')

    parser.add_argument('output_folder', type=str,
                        help='Output folder for processed files (must exist)')

    parser.add_argument('--chem_version', default='V2', choices=['V1', 'V2'],
                        help='Chemistry version (V1 or V2) (default: V2)')

    parser.add_argument('--r1_min_len', type=int, default=30,
                        help='Minimum Read 1 length after trimming (default: 30)')

    parser.add_argument('--r2_min_len', type=int, default=30,
                        help='Minimum Read 2 length after trimming (default: 30)')

    parser.add_argument('--single', action='store_true', default=False,
                        help='Option to process single-end sequencing run (default: paired-end)')

    args = parser.parse_args()  # parse arguments

    sample_label = args.sample_label
    output_folder = os.path.abspath(os.path.expanduser(args.output_folder))
    input_fastq_folder = os.path.abspath(os.path.expanduser(args.input_fastq_folder))
    chem_version = args.chem_version
    r1_min_len = args.r1_min_len
    r2_min_len = args.r2_min_len
    single = args.single
    if single:
        paired_end = 'Single'
    else:
        paired_end = 'Paired'

    # check that folders exist
    if not os.path.isdir(output_folder):
        raise OSError(output_folder + ' is not a directory. Exiting...')

    if not os.path.isdir(input_fastq_folder):
        raise OSError(input_fastq_folder + ' is not a directory. Exiting...')

    fastq_files = [input_fastq_folder + '/' + f for f in os.listdir(input_fastq_folder) if '.fastq.gz' in f]
    fastq_files.sort()

    ### experiment properties

    barcode_counts = output_folder + '/' + sample_label + '_barcode_counts.tsv'
    cell_barcode_cutadapt = output_folder + '/' + sample_label + '_cell_barcode_cutadapt.txt'
    fastq_out_r1 = output_folder + '/' + sample_label + '_R1.fastq'
    fastq_out_r2 = output_folder + '/' + sample_label + '_R2.fastq'

    # filtering parameters for cutadapt based on chemistry version
    if chem_version == 'V1':
        r1_start = 'CGATGACG'
        r1_end = 'CTGTCTCTTATACACATCT'
        r2_end = 'CGTCATCG'
        bar_ind_1, bar_ind_2 = range(8), range(-8, 0)
        cell_barcode_csv = sys.path[0] + '/resources/v1_barcodes.csv'

    elif chem_version == 'V2':
        r1_start = 'GTACTCGCAGTAGTC'
        r1_end = 'CTGTCTCTTATACACATCT'
        r2_end = 'GACTACTGCGAGTAC'
        bar_ind_1, bar_ind_2 = range(9), range(-9, 0)
        cell_barcode_csv = sys.path[0] + '/resources/v2_barcodes.csv'

    print '''
####################################################################################
# Step 1: get input file names and store in TapestriSample objects
####################################################################################
    '''

    if paired_end == 'Paired':

        # get all fastq filenames
        fastq_files.sort()
        R1_files = [fastq_files[i] for i in range(len(fastq_files)) if i % 2 == 0]
        R2_files = [fastq_files[i] for i in range(len(fastq_files)) if i % 2 == 1]

        assert len(R1_files) == len(R2_files), 'Number of R1 files does not match number of R2 files!'

        # check filenames
        for i in range(len(R1_files)):
            assert '_R1_' in R1_files[i], 'Bad R1 filename: %s' % R1_files[i]
            assert '_R2_' in R2_files[i], 'Bad R2 filename: %s' % R2_files[i]
            assert R1_files[i].split('_')[0] == R2_files[i].split('_')[0], 'Filename mismatch!'

        # store sample info in Sample objects
        samples = generate_samples(R1_files, R2_files, output_folder, sample_label)

        # display sample summary
        file_summary(samples)

    else:
        R1_files = [fastq_files[i] for i in range(len(fastq_files))]

        # check filenames
        for i in range(len(R1_files)):
            assert '_R1_' in R1_files[i], 'Bad R1 filename: %s' % R1_files[i]

        # store sample info in Sample objects
        samples = generate_samples(R1_files, R1_files, output_folder, sample_label)

        # display sample summary
        file_summary(samples)

    print '''
####################################################################################
# Step 2: filter reads for cell barcode, perform error correction, trim reads
####################################################################################
    '''

    # load mission bio barcode csv file
    barcodes = mb_resources.load_barcodes(cell_barcode_csv, 1, False)

    # generate hamming dictionary for error correction
    barcodes = mb_resources.generate_hamming_dict(barcodes)

    print 'Barcode sequences loaded into dictionary.\n'

    # for panel reads, filter reads with valid barcode structure and export to new fastq
    print 'Extracting barcodes from raw fastq files...\n'

    # process cell barcodes
    process_barcodes = []
    for sample in samples:
        p = Process(
            target=sample.filter_valid_reads,
            args=(r1_start,
                  barcodes,
                  bar_ind_1,
                  bar_ind_2))
        process_barcodes.append(p)
        p.start()

    # wait for processes to finish
    for p in process_barcodes:
        p.join()

    # cut adapters from reads and add barcodes to header
    # for each sample, split file into chunks and send to separate processes

    cut_adapters = []
    for sample in samples:
        p = Process(
            target=sample.barcode_reads,
            args=(r1_start,
                  r1_end,
                  r2_end,
                  r1_min_len,
                  r2_min_len,
                  paired_end))
        cut_adapters.append(p)
        p.start()

    # wait for processes to finish
    for p in cut_adapters:
        p.join()

    # combine fastq files and zip
    out_r1 = [s.r1_trimmed for s in samples]
    out_r2 = [s.r2_trimmed for s in samples]

    wait([subprocess.Popen('cat %s > %s' % (' '.join(out_r1), fastq_out_r1), shell=True)])
    wait([subprocess.Popen('cat %s > %s' % (' '.join(out_r2), fastq_out_r2), shell=True)])

    wait([subprocess.Popen('gzip -f %s' % fastq_out_r1, shell=True)])
    wait([subprocess.Popen('gzip -f %s' % fastq_out_r2, shell=True)])

    # combine cutadapt reports
    cell_cutadapt = [s.cell_cutadapt for s in samples]
    wait([subprocess.Popen('cat %s > %s' % (' '.join(cell_cutadapt), cell_barcode_cutadapt), shell=True)])

    # combine barcode count files
    bc = [s.barcode_counts for s in samples]
    wait([subprocess.Popen('cat %s > %s' % (' '.join(bc), barcode_counts), shell=True)])

    # delete temporary files
    [os.remove(s.r1_trimmed) for s in samples]
    [os.remove(s.r2_trimmed) for s in samples]
    [os.remove(s.cell_cutadapt) for s in samples]
    [os.remove(s.barcode_counts) for s in samples]
    [os.remove(s.read_id_json) for s in samples]
