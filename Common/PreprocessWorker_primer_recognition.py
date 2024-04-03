#####################################################################################################
# Author    YHLee
# Date      2021-04-26 ~
# Editor
# Revised
# Note      user module for multiprocessing in pre-processing script.
#           cover the case of light chain & RT primers binding head of C genes
#           version 1.1
#           applying MiXCR, not UMI processing
#
#####################################################################################################

import csv
import logging
import multiprocessing
import os
import re
from copy import deepcopy
import tempfile

from multiprocessing import Process
from datetime import datetime
from Common.MultiprocessLogger import QueueHandler
from Common.SeqUtil import Helper, is_nucleotide_or_protein
from Common.FileUtil import read_file, read_file_new, write_file
from Common.Enum import ChainType
from Common.ProcessUtil import run_igblast_new

# for memory usage test
import psutil

# New import for blast_run
from Bio.Blast.Applications import NcbiblastnCommandline


class PreprocessWorker(object):
    """ Do multiple works for preprocessing.
        step 1. Error correction
                1a) primer recognition -> return dict
                1b) apply MiXCR -> return error_correction file
        step 2. Isotype annotation (by BLAST)
        step 3. Ig-related Information annotation (by IgBLAST)
        step 4. functional reads filtration
    """

    def __init__(self, work_dir, sample_name_list, infile_suffix='', outfile_suffix='',
                 fwd_primers=None, rev_primers=None, multiprocess_log_queue=None, chain_type=None, file_type='tsv', debug=False):
        # Setup working directory
        if not os.path.isdir(work_dir):
            logging.info('Not valid working directory path...')
            raise RuntimeError
        else:
            self.dir = work_dir

        self.sample_list = sample_name_list
        self.log_queue = multiprocess_log_queue

        self.infile_suffix = infile_suffix
        self.outfile_suffix = outfile_suffix

        self.fwd_primers = fwd_primers
        self.rev_primers = rev_primers

        self.chain_type = chain_type

        if file_type in ['tsv', 'csv']:
            self.file_type = file_type
        else:
            logging.info('Not supported file type... Must be csv or tsv')
            raise ValueError
        self.delim = ',' if file_type=='csv' else '\t'

        self.debug = debug

        self.initial_RAM_used = psutil.virtual_memory().used / 1024 ** 3


    def run(self, num_process):
        set_num_process = max(5, num_process)
        if set_num_process > 30:
            set_num_process = 30                # max number of process is limited to 30 (Junholab server)
        manager = multiprocessing.Manager()

        # for the RAM usage check
        max_RAM_used = psutil.virtual_memory().used / 1024 ** 3

        chunk_size_max = 20000000

        # copy the sample_list. Considering that whole steps will be done at once...
        samples = self.sample_list.copy()

        while True:
            # Do works according to step
            if len(samples) == 0:
                break
            self.sample_name_to_put = samples.pop(0)
            sample_dir = os.path.join(self.dir, self.sample_name_to_put)


            start_time_main = datetime.now()
            logging.info('Started Error Correction for %s', self.sample_name_to_put)

            in_dir = sample_dir
            out_dir = os.path.join(sample_dir, '1_error_correction')
            if not os.path.exists(out_dir):
                os.mkdir(out_dir)
            target_file = os.path.join(in_dir, '%s.csv' % self.sample_name_to_put)

            # do primer recognition
            if True:
                RAM_used_queue = manager.Queue()  # RAM_used queue
                start_time = datetime.now()

                # get chunk size according to num_process and put it to the queue
                file_line_len = self.get_file_lines(target_file) - 1        # exclude header line (csv file)
                chunk_size = int(file_line_len / set_num_process) + 1
                if chunk_size > file_line_len:
                    chunk_size = file_line_len
                elif chunk_size > chunk_size_max:
                    chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                chunk_start_pos_list = self.extract_chunk_start_pos(file_line_len, chunk_size)


                # 1. do primer recognition with fwd primers
                if True:
                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    sub_step = 'a'
                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step a.", set_num_process)
                    extract_workers = []
                    for i in range(set_num_process):
                        p = ErrorCorrecter(sub_step, result_queue, self.log_queue, RAM_used_queue, in_queue=target_queue, chunk_size=chunk_size,
                                           input_file=target_file, _dir=in_dir, fwd_primer_list=self.fwd_primers, rev_primer_list=self.rev_primers)

                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    primer_recognized_dict = self.merge_chunk_dict_type2(out_result_queue)

                    # temporary file write for the comparison
                    if self.debug:
                        outfile_dir = sample_dir
                        output_file = os.path.join(outfile_dir, '%s_%s_fwd_primer_recognized.csv' % (self.sample_name_to_put, self.outfile_suffix))

                        w_header = ['full_NT', 'readcount']
                        with open(output_file, 'w') as handle:
                            write_file(handle, w_header, list(primer_recognized_dict.values()))

                    logging.info('a. Fwd primer recognition was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                # 2. do primer recognition with rev primers
                if True:
                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    sub_step = 'b'
                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step b.", set_num_process)
                    extract_workers = []
                    for i in range(set_num_process):
                        p = ErrorCorrecter(sub_step, result_queue, self.log_queue, RAM_used_queue, in_queue=target_queue, chunk_size=chunk_size,
                                           input_file=target_file, _dir=in_dir, fwd_primer_list=self.fwd_primers, rev_primer_list=self.rev_primers)

                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    primer_recognized_dict = self.merge_chunk_dict_type2(out_result_queue)

                    # temporary file write for the comparison
                    if self.debug:
                        outfile_dir = sample_dir
                        output_file = os.path.join(outfile_dir, '%s_%s_rev_primer_recognized.csv' % (self.sample_name_to_put, self.outfile_suffix))

                        w_header = ['full_NT', 'readcount']
                        with open(output_file, 'w') as handle:
                            write_file(handle, w_header, list(primer_recognized_dict.values()))

                    logging.info('b. Rev primer recognition was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))


            logging.info("Finished Error Correction for %s", self.sample_name_to_put)
            logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time_main))


        logging.info("Primer recognition for all files are finished")
        logging.info('Max RAM used : %.2f GB' % (max_RAM_used - self.initial_RAM_used))


    ### work functions for the ProcessWorker class
    def get_file_lines(self, target):
        count = 0
        with open(target, 'r') as handle:
            for line in handle:
                count += 1
        return count

    def get_key_len(self, target_dict):
        count = len(target_dict)
        return count

    def extract_chunk_start_pos(self, target_len, chunk_size):
        chunk_start_list = []

        i = 0
        while True:
            if chunk_size * i > target_len:
                break

            line_start = chunk_size * i
            i += 1

            chunk_start_list.append(line_start)

        logging.info("Size of divided chunks for target %s: %s" % (self.sample_name_to_put, chunk_size))

        return chunk_start_list


    # special function to merge consensus sequence chunk dictionaries & merge isotype chunk dictionaries.
    def merge_chunk_dict_type2(self, in_result_queue, rc_col_name='readcount'):
        # dict form
        # {key1: {}, key2: {}, ...}

        result_dict = {}
        while in_result_queue.qsize() > 0:
            result_chunk = in_result_queue.get()
            for k in result_chunk:
                if k not in result_dict:
                    result_dict[k] = result_chunk[k]
                else:
                    result_dict[k][rc_col_name] += result_chunk[k][rc_col_name]

        return result_dict



class ErrorCorrecter(Process):
    """ Do error correction (UMI processing)

    1a. primer recognition : file -> dict
    1b. applying MiXCR: dict -> file (no MP)
    """

    def __init__(self, sub_step, out_queue, log_queue, RAM_used_queue, in_queue='', chunk_size='', input_file='', input_dict={}, _dir=None,
                 fwd_primer_list=None, rev_primer_list=None):
        Process.__init__(self)

        self.input_file = input_file

        self.out_chunk = {}

        self.sub_step = sub_step
        self.out_queue = out_queue
        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.in_queue = in_queue
        self.chunk_size = chunk_size

        self.input_dict = input_dict
        self.dir = _dir
        self.fwd_primer_list = fwd_primer_list
        self.rev_primer_list = rev_primer_list


    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        # logging.info("Started run")

        while True:
            # Get the work from the queue and expand the tuple
            chunk_start_pos = self.in_queue.get()
            if chunk_start_pos is None:
                self.in_queue.task_done()
                break

            if self.sub_step == 'a':
                self.primer_recognition_f(chunk_start_pos, self.input_file)
            elif self.sub_step == 'b':
                self.primer_recognition_r(chunk_start_pos, self.input_file)

            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

            self.in_queue.task_done()

            # save chunk data into result queue
            self.out_queue.put(self.out_chunk)


        # logging.info("Finished run")


    # worker functions
    def primer_recognition_f(self, line_start, input_file_name):
        sHelper = Helper()

        # set the direction of primers into sense strand (plus/plus strand) (fwd, rev)
        forwardPrimerList = self.fwd_primer_list


        # add forward primer sequences containing 1 mismatch
        forwardPrimerMatchList = []
        for forwardPrimer in forwardPrimerList:
            forwardPrimerMatchList.append(forwardPrimer)
        for forwardPrimer in forwardPrimerList:
            for i in range(len(forwardPrimer)):
                forwardPrimerMatchList.append(forwardPrimer[:i] + forwardPrimer[i + 1:])
                forwardPrimerMatchList.append(forwardPrimer[:i] + '[ATGC]' + forwardPrimer[i + 1:])

        # set the forward / reverse region to find pattern
        fwd_primer_max_len = max([len(fwd) for fwd in forwardPrimerList])
        forwardTrimNum = fwd_primer_max_len + 5

        recogSet = {}

        with open(os.path.join(input_file_name), 'r') as f:
            reader = csv.reader(f)
            header = next(reader)

            if len(header) != 2:
                print('Input file column error')
                raise RuntimeError

            # go to start position
            for i in range(0, line_start):
                next(reader)

            line_count = 0
            accessed_line_num = line_start
            while True:
                try:
                    row = next(reader)
                except StopIteration:
                    break

                line_count += 1
                if line_count > self.chunk_size:
                    break

                sequence = row[0]
                readcount = int(row[1])

                success = False

                for forwardPrimerMatch in forwardPrimerMatchList:
                    forwardMatchResult = re.search(forwardPrimerMatch, sequence[:forwardTrimNum])

                    if forwardMatchResult != None:
                        # newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[0])]  # erase primer region
                        newSequence = sequence
                        success = True
                        recogSet[accessed_line_num] = {'full_NT': newSequence, 'readcount': readcount}
                        accessed_line_num += 1
                        break
                    else:
                        continue

                # try primer recognition process for reverse complement case.
                if success == False:
                    sequence = sHelper.to_reverse_complement(sequence)

                    for forwardPrimerMatch in forwardPrimerMatchList:
                        forwardMatchResult = re.search(forwardPrimerMatch, sequence[:forwardTrimNum])

                        if forwardMatchResult != None:
                            # newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[0])]
                            newSequence = sequence
                            success = True
                            recogSet[accessed_line_num] = {'full_NT': newSequence, 'readcount': readcount}
                            accessed_line_num += 1
                            break
                        else:
                            continue

        self.out_chunk = recogSet

    def primer_recognition_r(self, line_start, input_file_name):
        sHelper = Helper()

        # set the direction of primers into sense strand (plus/plus strand) (fwd, rev)
        reversePrimerList = [sHelper.to_reverse_complement(reversePrimer) for reversePrimer in self.rev_primer_list]

        # add reverse primer sequences containing 1 mismatch (regular expression pattern)
        reversePrimerMatchList = []
        for reversePrimer in reversePrimerList:
            reversePrimerMatchList.append(reversePrimer)
        for reversePrimer in reversePrimerList:
            for i in range(len(reversePrimer)):
                reversePrimerMatchList.append(reversePrimer[:i] + reversePrimer[i + 1:])
                reversePrimerMatchList.append(reversePrimer[:i] + '[ATGC]' + reversePrimer[i + 1:])


        # set the forward / reverse region to find pattern
        rev_primer_max_len = max([len(rev) for rev in reversePrimerList])
        reverseTrimNum = rev_primer_max_len + 14 + 16

        recogSet = {}

        with open(os.path.join(input_file_name), 'r') as f:
            reader = csv.reader(f)
            header = next(reader)

            if len(header) != 2:
                print('Input file column error')
                raise RuntimeError

            # go to start position
            for i in range(0, line_start):
                next(reader)

            line_count = 0
            accessed_line_num = line_start
            while True:
                try:
                    row = next(reader)
                except StopIteration:
                    break

                line_count += 1
                if line_count > self.chunk_size:
                    break

                sequence = row[0]
                readcount = int(row[1])

                success = False

                for reversePrimerMatch in reversePrimerMatchList:
                    reverseMatchResult = re.search(reversePrimerMatch, sequence[-reverseTrimNum:])

                    if reverseMatchResult != None:
                        newSequence = sequence
                        success = True
                        recogSet[accessed_line_num] = {'full_NT': newSequence, 'readcount': readcount}
                        accessed_line_num += 1
                        break
                    else:
                        continue

                # try primer recognition process for reverse complement case.
                if success == False:
                    sequence = sHelper.to_reverse_complement(sequence)

                    for reversePrimerMatch in reversePrimerMatchList:
                        reverseMatchResult = re.search(reversePrimerMatch, sequence[-reverseTrimNum:])

                        if reverseMatchResult != None:
                            newSequence = sequence
                            success = True
                            recogSet[accessed_line_num] = {'full_NT': newSequence, 'readcount': readcount}
                            accessed_line_num += 1
                            break
                        else:
                            continue

        self.out_chunk = recogSet


    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue
