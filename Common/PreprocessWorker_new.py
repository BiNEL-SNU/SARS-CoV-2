#####################################################################################################
# Author    YHLee
# Date      2021-04-26 ~
# Editor
# Revised
# Note      user module for multiprocessing in pre-processing script.
#           cover the case of light chain & RT primers binding head of C genes
#           version 1.1
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
                1a-1) primer recognition -> return dict
                1a-2) UMI merge -> return dict
                1b) sub-clustering -> return dict
                1c) length filtering -> return dict
                1d) consensus sequence extraction -> return error_correction file
        step 2. Isotype annotation (by BLAST)
        step 3. Ig-related Information annotation (by IgBLAST)
        step 4. functional reads filtration
    """

    def __init__(self, work_dir, sample_name_list, infile_suffix='', outfile_suffix='',
                 fwd_primers=None, rev_primers=None, c_ref_file=None, dist_thresh=5, hamming_thresh=1, vote_thresh=0.6,
                 step=None, multiprocess_log_queue=None, chain_type=None, file_type='tsv', debug=False, rev_primer_dist_thresh=2, subtyping=True):
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
        self.c_ref_file = c_ref_file
        self.dist_thresh = dist_thresh
        self.hamming_thresh = hamming_thresh
        self.vote_thresh = vote_thresh
        self.step = step

        self.chain_type = chain_type

        if file_type in ['tsv', 'csv']:
            self.file_type = file_type
        else:
            logging.info('Not supported file type... Must be csv or tsv')
            raise ValueError
        self.delim = ',' if file_type=='csv' else '\t'

        self.debug = debug
        self.rev_primer_dist_thresh = rev_primer_dist_thresh
        self.subtyping = subtyping

        self.initial_RAM_used = psutil.virtual_memory().used / 1024 ** 3


    def run(self, num_process):
        set_num_process = max(2, num_process)
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

            # case separation: Step 1, 2, 3, 4
            if self.step == 1:
                start_time_main = datetime.now()
                logging.info('Started Error Correction for %s', self.sample_name_to_put)

                in_dir = sample_dir
                out_dir = os.path.join(sample_dir, '1_error_correction')
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)
                target_file = os.path.join(in_dir, '%s.csv' % self.sample_name_to_put)
                out_file = os.path.join(out_dir, '%s_%s.%s' % (self.sample_name_to_put, self.outfile_suffix, self.file_type))

                ##################################################################################
                ######## sub-step: 1a-1. primer recognition ######################################
                ##################################################################################
                if True:
                    # 1a-1. primer recognition
                    sub_step = 'a-1'
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

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step 1a-1.", set_num_process)
                    extract_workers = []
                    for i in range(set_num_process):
                        p = ErrorCorrecter(target_queue, chunk_size, sub_step, result_queue, self.log_queue, RAM_used_queue,
                                           input_file=target_file, _dir=in_dir, fwd_primer_list=self.fwd_primers, rev_primer_list=self.rev_primers,
                                           chain_type=self.chain_type)

                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    logging.info('1a-1. Primer recognition was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    primer_recognized_dict = self.merge_chunk_dict_type1(out_result_queue)

                    # temporary file write for the comparison
                    if self.debug:
                        outfile_dir = sample_dir
                        output_file = os.path.join(outfile_dir, '%s_%s_primer_recognized.csv' % (self.sample_name_to_put, self.outfile_suffix))

                        w_header = ['UMI', 'full_NT', 'readcount', 'rev_primer_len']
                        w_data = []
                        for umi in primer_recognized_dict:
                            for seq_dict in primer_recognized_dict[umi]:
                                w_data.append(seq_dict)

                        with open(output_file, 'w') as handle:
                            write_file(handle, w_header, w_data)

                ##################################################################################
                ######## sub-step: 1a-2. UMI merge ###############################################
                ##################################################################################
                if True:
                    # 1a-2. UMI merge
                    sub_step = 'a-2'
                    start_time = datetime.now()

                    # check if primer recognized dict is empty or not
                    if len(primer_recognized_dict) == 0:
                        logging.info('All reads were excluded during primer recognition step... pass following processes')
                        continue

                    sorted_result_dict = self.sort_UMI_with_size(primer_recognized_dict)

                    UMI_list = list(sorted_result_dict.keys())

                    # get chunk size according to num_process and put it to the queue
                    whole_UMI_size = len(UMI_list)
                    chunk_size = int(whole_UMI_size / set_num_process) + 1
                    if chunk_size > whole_UMI_size:
                        chunk_size = whole_UMI_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(whole_UMI_size, chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step 1a-2.", set_num_process)
                    extract_workers = []
                    for i in range(set_num_process):
                        p = ErrorCorrecter(target_queue, chunk_size, sub_step, result_queue, self.log_queue, RAM_used_queue,
                                           input_dict=sorted_result_dict, input_UMI_list=UMI_list)

                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # merge result and update output dictionary if possible
                    out_result_queue = p.get_result_queue()
                    logging.info('# of UMIs extracted :%d' % whole_UMI_size)
                    UMI_survived_hash = self.get_UMI_be_merged_hash(out_result_queue)       # UMI_survived_hash = {UMI1: [UMIi, UMIj], UMI4: [UMIm], ...}
                    UMI_merged_dict = self.get_UMI_merged_dict(sorted_result_dict, UMI_survived_hash)

                    # temporary file write for the comparison
                    if self.debug:
                        outfile_dir = sample_dir
                        output_file = os.path.join(outfile_dir, '%s_%s_UMI_merged.csv' % (self.sample_name_to_put, self.outfile_suffix))

                        w_header = ['mergedUMI', 'UMI', 'full_NT', 'readcount', 'rev_primer_len']
                        w_data = []
                        for merged_umi in UMI_merged_dict:
                            for d in UMI_merged_dict[merged_umi]:
                                d['mergedUMI'] = merged_umi
                                w_data.append(d)

                        with open(output_file, 'w') as handle:
                            write_file(handle, w_header, w_data)

                    logging.info('1a-2. UMI merge was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                    # get RAM_used_queue
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                ##################################################################################
                ######## sub-step: 1b. sub-clustering ############################################
                ##################################################################################
                if True:
                    sub_step = 'b'
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # get chunk size according to num_process and put it to the queue
                    UMI_merged_dict_size = self.get_key_len(UMI_merged_dict)
                    chunk_size = int(UMI_merged_dict_size / set_num_process) + 1
                    if chunk_size > UMI_merged_dict_size :
                        chunk_size = UMI_merged_dict_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(UMI_merged_dict_size , chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step 1b.", set_num_process)
                    extract_workers = []
                    for i in range(set_num_process):
                        p = ErrorCorrecter(target_queue, chunk_size, sub_step, result_queue, self.log_queue, RAM_used_queue,
                                           input_dict=UMI_merged_dict, dist_thresh=self.dist_thresh, rev_primer_dist_thresh=self.rev_primer_dist_thresh)
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    logging.info('1b. Sub-clustering was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    subcluster_dict = self.merge_chunk_dict_type1(out_result_queue)

                    # get RAM_used_queue
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                # temporary ftn for the check
                if self.debug:
                    self.tmp_write(subcluster_dict, sample_dir)

                ##################################################################################
                ######## sub-step: 1c. length_filtering ##########################################
                ##################################################################################
                if True:
                    sub_step = 'c'
                    num_process_for_substep_c = 10
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # check if subclustered dict is empty or not
                    if len(subcluster_dict) == 0:
                        logging.info('All reads were excluded during sub-clustering step... pass following processes')
                        continue

                    # get chunk size according to num_process and put it to the queue
                    subcluster_dict_size = self.get_key_len(subcluster_dict)
                    chunk_size = int(subcluster_dict_size / num_process_for_substep_c) + 1
                    if chunk_size > subcluster_dict_size:
                        chunk_size = subcluster_dict_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(subcluster_dict_size, chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step 1c.", num_process_for_substep_c)
                    extract_workers = []
                    for i in range(num_process_for_substep_c):
                        p = ErrorCorrecter(target_queue, chunk_size, sub_step, result_queue, self.log_queue, RAM_used_queue,
                                           input_dict=subcluster_dict, hamming_thresh=self.hamming_thresh)
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    logging.info('1c. Length filtering was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    len_filtered_dict = self.merge_chunk_dict_type1(out_result_queue)

                    # get RAM_used_queue
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                ##################################################################################
                ######## sub-step: 1d. consensus_sequence_extraction #############################
                ##################################################################################
                if True:
                    sub_step = 'd'
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # check if length filtered dict is empty or not
                    if len(len_filtered_dict) == 0:
                        logging.info('All reads were excluded during length filtering step... pass following processes')
                        continue

                    # get chunk size according to num_process and put it to the queue
                    len_filtered_dict_size = self.get_key_len(len_filtered_dict)
                    chunk_size = int(len_filtered_dict_size / set_num_process) + 1
                    if chunk_size > len_filtered_dict_size:
                        chunk_size = len_filtered_dict_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(len_filtered_dict_size, chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    # start ErrorCorrecter processes (MP)
                    logging.info("Spawning %d ErrorCorrecter processes for step 1d.", set_num_process)
                    extract_workers = []
                    for i in range(set_num_process):
                        p = ErrorCorrecter(target_queue, chunk_size, sub_step, result_queue, self.log_queue, RAM_used_queue,
                                           input_dict=len_filtered_dict, vote_thresh=self.vote_thresh)
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    logging.info('1d. Consensus sequence extraction was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                    # merge result into one dict
                    out_result_queue = p.get_result_queue()
                    error_corrected_dict = self.merge_chunk_dict_with_UMIs(out_result_queue)

                    # get RAM_used_queue
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    # write result files
                    w_header = ['full_NT', 'umi_list', 'readcount']
                    with open(out_file, 'w') as handle:
                        write_file(handle, w_header, list(error_corrected_dict.values()), delim=self.delim)

                logging.info("Finished Error Correction for %s", self.sample_name_to_put)
                logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time_main))


            elif self.step == 2:
                start_time_main = datetime.now()
                logging.info('Started Isotyping for %s', self.sample_name_to_put)

                in_dir = os.path.join(sample_dir, '1_error_correction')
                out_dir = os.path.join(sample_dir, '2_annotation')
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                target_file = os.path.join(in_dir, '%s_%s.%s' % (self.sample_name_to_put, self.infile_suffix, self.file_type))
                out_file = os.path.join(out_dir, '%s_%s.%s' % (self.sample_name_to_put, self.outfile_suffix, self.file_type))

                # check if error-corrected file exists or not
                if not os.path.exists(target_file):
                    logging.info('Error corrected file does not exist... pass following processes')
                    continue

                ##################################################################################
                ######## sub-step: 2a. run blast on C gene reference #############################
                ##################################################################################
                if True:
                    sub_step = 'a'
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    ############################## started MP #########
                    ## Applying the multiprocessing that the BLAST tool supports
                    logging.info("Do blast run (step 2a) of IsotypeExtractor with %d processes", set_num_process)
                    extract_workers = []
                    p = IsotypeExtractor(target_file, sub_step, self.log_queue, RAM_used_queue,
                                         out_dir=out_dir, c_ref_file=self.c_ref_file, num_process=set_num_process, delim=self.delim)
                    p.daemon = True
                    p.start()
                    extract_workers.append(p)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # RAM usage update -> not for sure... because it cannot monitor the inside the tool BLAST...
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info('2a. Run blast on C gene reference was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                ##################################################################################
                ######## sub-step: 2b. Parsing blast result file #################################
                ##################################################################################
                if True:
                    sub_step = 'b'
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # blast_result out directory
                    blast_result_dir = os.path.join(out_dir + 'temp_cblast', 'result')

                    # get chunk size according to num_process and put it to the queue
                    blast_parsing_iterator_size = self.get_file_lines(target_file) - 1      # exclude header
                    chunk_size = int(blast_parsing_iterator_size / set_num_process) + 1
                    if chunk_size > blast_parsing_iterator_size:
                        chunk_size = blast_parsing_iterator_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(blast_parsing_iterator_size, chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    logging.info("Do parse isotype info (step 2b) of IsotypeExtracter with %d processes", set_num_process)
                    extract_workers = []

                    for i in range(set_num_process):
                        p = IsotypeExtractor(target_file, sub_step, self.log_queue, RAM_used_queue,
                                             out_dir=blast_result_dir, in_queue=target_queue, chain_type=self.chain_type,
                                             out_queue=result_queue, chunk_size=chunk_size, delim=self.delim, subtyping=self.subtyping)
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # combine all chunk result, sort following the line order of original input file (error-corrected file), and merge
                    out_result_queue = p.get_result_queue()
                    isotype_added_dict = self.merge_chunk_dict_type2(out_result_queue)

                    # write file
                    w_header = ['full_NT', 'umi_list', 'isotype', 'readcount']
                    with open(out_file, 'w') as handle:
                        write_file(handle, w_header, list(isotype_added_dict.values()), delim=self.delim)

                    # RAM usage update
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info('2b. Parse isotype information was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                logging.info("Finished Isotyping for %s", self.sample_name_to_put)
                logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time_main))


            elif self.step == 3:
                start_time_main = datetime.now()
                logging.info('Started Annotation for %s', self.sample_name_to_put)

                in_dir = os.path.join(sample_dir, '2_annotation')
                out_dir = in_dir
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                target_file = os.path.join(in_dir, '%s_%s.%s' % (self.sample_name_to_put, self.infile_suffix, self.file_type))
                out_file = os.path.join(out_dir, '%s_%s.%s' % (self.sample_name_to_put, self.outfile_suffix, self.file_type))

                # check if C gene annotated file exists or not
                if not os.path.exists(target_file):
                    logging.info('C gene annotated file does not exist... pass following processes')
                    continue

                # count line number of file and continue
                input_file_line = self.get_file_lines(target_file)
                if input_file_line == 1:
                    logging.info('Empty input file. pass the annotation step...')
                    continue

                ##################################################################################
                ######## sub-step: 3a. run IgBLAST  ##############################################
                ##################################################################################
                if True:
                    sub_step = 'a'
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    ############################## started MP #########
                    ## Applying the multiprocessing that the BLAST tool supports
                    logging.info("Do IgBLAST run (step 3a) of InfoAnnotator with %d processes", set_num_process)
                    extract_workers = []
                    p = InfoAnnotator(target_file, sub_step, self.log_queue, RAM_used_queue, out_dir, self.chain_type,
                                      num_process=set_num_process, delim=self.delim)
                    p.daemon = True
                    p.start()
                    extract_workers.append(p)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # RAM usage update -> not for sure... because it cannot monitor the inside the tool BLAST...
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info('3a. Run IgBLAST was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                ##################################################################################
                ######## sub-step: 3b. Parsing igblast result file ###############################
                ##################################################################################
                if True:
                    sub_step = 'b'
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # igblast_result out directory
                    igblast_result_dir = os.path.join(out_dir, 'temp_igblast')

                    # get chunk size according to num_process and put it to the queue
                    igblast_parsing_size = self.get_file_lines(target_file) - 1  # exclude header
                    chunk_size = int(igblast_parsing_size / set_num_process) + 1
                    if chunk_size > igblast_parsing_size:
                        chunk_size = igblast_parsing_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(igblast_parsing_size, chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    logging.info("Do info annotation (step 3b) of InfoAnnotator with %d processes", set_num_process)
                    extract_workers = []

                    for i in range(set_num_process):
                        p = InfoAnnotator(target_file, sub_step, self.log_queue, RAM_used_queue, igblast_result_dir, self.chain_type,
                                          in_queue=target_queue, out_queue=result_queue, chunk_size=chunk_size, delim=self.delim)
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # combine all chunk result, sort following the line order of original input file (error-corrected file), and merge
                    out_result_queue = p.get_result_queue()
                    annotated_dict = self.merge_chunk_dict_type2(out_result_queue)
                    # annotated_dict_origin_order = {k:v for k, v in sorted(annotated_dict.items(), key=lambda x:x[0], reverse=False)}
                    sorted_annotated_dict = {k:v for k, v in sorted(annotated_dict.items(), key=lambda x:x[1]['duplicate_count'], reverse=True)}


                    # write file
                    w_header = ['duplicate_count', 'sequence', 'sequence_aa',
                                'sequence_alignment','germline_alignment', 'v_call', 'd_call', 'j_call', 'c_call',
                                'junction', 'junction_aa','cdr1', 'cdr2', 'cdr3', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa',
                                'rev_comp', 'productive',
                                'v_cigar', 'd_cigar', 'j_cigar', 'v_alignment_length', 'v_alignment_mutation', 'umi_list']
                    if self.chain_type not in [ChainType.HUMAN_HEAVY, ChainType.MOUSE_HEAVY]:
                        w_header.remove('d_call')
                        w_header.remove('d_cigar')
                    with open(out_file, 'w') as handle:
                        write_file(handle, w_header, list(sorted_annotated_dict.values()), delim=self.delim)

                    # RAM usage update
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info('3b. Info annotation was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                logging.info("Finished Annotation for %s", self.sample_name_to_put)
                logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time_main))


            elif self.step == 4:
                start_time_main = datetime.now()
                logging.info('Started Functionality check for %s', self.sample_name_to_put)

                in_dir = os.path.join(sample_dir, '2_annotation')
                out_dir = os.path.join(sample_dir, '3_functionality')
                if not os.path.exists(out_dir):
                    os.mkdir(out_dir)

                target_file = os.path.join(in_dir, '%s_%s.%s' % (self.sample_name_to_put, self.infile_suffix, self.file_type))
                out_file = os.path.join(out_dir, '%s_%s.%s' % (self.sample_name_to_put, self.outfile_suffix, self.file_type))

                # check if annotated file exists or not
                if not os.path.exists(target_file):
                    logging.info('Annotated file does not exist... pass following processes')
                    continue

                # count line number of file and continue
                input_file_line = self.get_file_lines(target_file)
                if input_file_line == 1:
                    logging.info('Empty input file. pass the functional reads filtration step...')
                    continue

                if True:
                    RAM_used_queue = manager.Queue()  # RAM_used queue
                    start_time = datetime.now()

                    # get chunk size according to num_process and put it to the queue
                    target_file_size = self.get_file_lines(target_file) - 1  # exclude header
                    chunk_size = int(target_file_size / set_num_process) + 1
                    if chunk_size > target_file_size:
                        chunk_size = target_file_size
                    elif chunk_size > chunk_size_max:
                        chunk_size = int(chunk_size / (int(chunk_size / chunk_size_max) + 1))
                    chunk_start_pos_list = self.extract_chunk_start_pos(target_file_size, chunk_size)

                    target_queue = manager.Queue()
                    for chunk_range in chunk_start_pos_list:
                        target_queue.put(chunk_range)

                    # make queue to get result data from processes
                    result_queue = manager.Queue()

                    ############################## started MP #########
                    logging.info("Do functional reads extraction (step 4) of functionalExtracter with %d processes", set_num_process)
                    extract_workers = []

                    for i in range(set_num_process):
                        p = functionalExtracter(target_file, self.log_queue, RAM_used_queue,
                                                in_queue=target_queue, out_queue=result_queue, chunk_size=chunk_size, delim=self.delim)
                        p.daemon = True
                        p.start()

                        extract_workers.append(p)
                        target_queue.put_nowait(None)

                    # wait until all extracting works done
                    for p in extract_workers:
                        p.join()
                    ############################## finished MP #########

                    # combine all chunk result, sort following the line order of readcount (duplicate_count), and merge
                    out_result_queue = p.get_result_queue()
                    functional_dict = self.merge_chunk_dict_type2(out_result_queue)
                    sorted_functional_dict = {k: v for k, v in sorted(functional_dict.items(), key=lambda x: int(x[1]['duplicate_count']), reverse=True)}

                    # set sequence_id
                    row_num = 1
                    for k in sorted_functional_dict:
                        sorted_functional_dict[k]['sequence_id'] = '%s-%d' % (self.sample_name_to_put, row_num)
                        row_num += 1

                    # write file
                    w_header = ['sequence_id', 'duplicate_count', 'sequence', 'sequence_aa',
                                'sequence_alignment', 'germline_alignment', 'v_call', 'd_call', 'j_call', 'c_call',
                                'junction', 'junction_aa', 'cdr1', 'cdr2', 'cdr3', 'cdr1_aa', 'cdr2_aa', 'cdr3_aa',
                                'rev_comp', 'productive',
                                'v_cigar', 'd_cigar', 'j_cigar', 'v_alignment_length', 'v_alignment_mutation', 'umi_list']
                    if self.chain_type not in [ChainType.HUMAN_HEAVY, ChainType.MOUSE_HEAVY]:
                        w_header.remove('d_call')
                        w_header.remove('d_cigar')
                    with open(out_file, 'w') as handle:
                        write_file(handle, w_header, list(sorted_functional_dict.values()), delim=self.delim)

                    # RAM usage update
                    out_RAM_used_queue = p.get_RAM_used_queue()
                    while out_RAM_used_queue.qsize() > 0:
                        RAM_used = out_RAM_used_queue.get()
                        if RAM_used > max_RAM_used:
                            max_RAM_used = RAM_used

                    logging.info('4. Functional reads extraction was done.')
                    logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time))

                logging.info("Finished Functionality check for %s", self.sample_name_to_put)
                logging.info('--- %s seconds elapsed ---' % (datetime.now() - start_time_main))


            elif self.step == None:
                print('Need to set "step" parameter... (int type: 1, 2, 3, 4)')
                return ValueError

        logging.info("Step %d for all files are finished" % self.step)
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

    def merge_chunk_dict_type1(self, in_result_queue):
        # dict form
        # {key1: [{}, {}, {}, ...], key2: [{}, {}, ...], ...}

        result_dict = {}
        while in_result_queue.qsize() > 0:
            result_chunk = in_result_queue.get()
            for k in result_chunk:
                if k not in result_dict:
                    result_dict[k] = result_chunk[k]
                else:
                    result_dict[k] += result_chunk[k]

        # check if the result_dict is empty or not
        if len(result_dict) == 0:
            logging.info('Output_chunks are all empty... return null dictionary')
            return result_dict

        ## for the case that key is a tuple -> primer recognized dict
        example_key = list(result_dict.keys())[0]
        if type(example_key) == tuple:
            # a. sort chunk data following the original raw data order
            sorted_result_dict = {k:v for k, v in sorted(result_dict.items(), key=lambda x:x[0][1], reverse=False)}

            # b. merge data with UMI
            tmp_dict = {}
            tmp_order_num = 0
            for umi, start_num in sorted_result_dict:
                if tmp_order_num > start_num:
                    logging.info('Data access with the order of original raw data failed...')
                    return
                if umi not in tmp_dict:
                    tmp_dict[umi] = sorted_result_dict[(umi, start_num)]
                else:
                    tmp_dict[umi] += sorted_result_dict[(umi, start_num)]
            result_dict = tmp_dict.copy()

        return result_dict

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

    # special function to merge consensus sequence chunk dictionaries
    def merge_chunk_dict_with_UMIs(self, in_result_queue, rc_col_name='readcount', umi_list_name='UMI_list'):
        # dict form
        # {key1: {k:v, k:v, k:[]}, key2: {k:v, k:v, k:[]}, ...}

        result_dict = {}
        while in_result_queue.qsize() > 0:
            result_chunk = in_result_queue.get()
            for k in result_chunk:
                if k not in result_dict:
                    result_dict[k] = result_chunk[k]
                else:
                    result_dict[k][rc_col_name] += result_chunk[k][rc_col_name]
                    result_dict[k][umi_list_name] += result_chunk[k][umi_list_name]

        # make umi_list into string
        for k in result_dict:
            result_dict[k]['umi_list'] = '|'.join(list(set(result_dict[k][umi_list_name])))

        return result_dict


    ### functions for UMI merge process
    def sort_UMI_with_size(self, UMI_dict):
        # sort UMISet according to readcount sum -> sortedUMISet is acquired.
        rcUMISet = {}

        for UMI in UMI_dict:
            readcountSum = 0
            for each in UMI_dict[UMI]:
                readcountSum += int(each['readcount'])
            try:
                rcUMISet[str(readcountSum)].append(UMI)
            except KeyError:
                rcUMISet[str(readcountSum)] = [UMI]

        rcList = [int(x) for x in list(set(rcUMISet.keys()))]
        rcList.sort(reverse=True)

        sortedUMISet = {}

        for rc in rcList:
            for rcUMI in rcUMISet[str(rc)]:
                sortedUMISet[rcUMI] = UMI_dict[rcUMI]

        del rcList
        del rcUMISet

        return sortedUMISet

    def get_UMI_be_merged_hash(self, in_result_queue):
        UMI_survived_hash = {}

        # step 1-1. make whole hash -> hash_chunk = {(merged UMI1, chunk_start): [UMIi, UMIj, ...], (merged UMI2, chunk_start): [UMIm, UMIn, ...], ...}
        while in_result_queue.qsize() > 0:
            UMI_merged_hash_chunk = in_result_queue.get()
            UMI_survived_hash.update(UMI_merged_hash_chunk)

        # step 1-2. sort the whole hash following the input (primer_recognized_dict) order
        sorted_UMI_survived_hash = {k[0]:v for k, v in sorted(UMI_survived_hash.items(), key=lambda x:x[0][1], reverse=False)}

        ## step 2 : already-merged UMIs deletion process
        # step 2-1. delete target UMI list if the key UMI already belongs to primary key UMI & delete target UMI in the list if it already merged to primary key UMI.
        checked_UMI_dict = {}
        sorted_UMI_survived_hash_out = deepcopy(sorted_UMI_survived_hash)

        for key_UMI, target_list in sorted_UMI_survived_hash.items():
            if len(target_list) > 0:
                if key_UMI in checked_UMI_dict:
                    sorted_UMI_survived_hash_out[key_UMI] = []
                    continue

                tmp_len = len(target_list)
                tmp_count = 0
                for target_UMI in target_list:
                    if target_UMI in checked_UMI_dict:
                        tmp_count += 1
                        sorted_UMI_survived_hash_out[key_UMI].remove(target_UMI)
                    else:
                        tmp_count += 1
                        checked_UMI_dict[target_UMI] = 0

                if tmp_len > tmp_count:
                    logging.info('Tracking target merged UMIs failed...')
                    return

        logging.info('Delete process 1 finished')
        logging.info('Size of whole merged UMIs: %d' % len(checked_UMI_dict))


        # step 2-2. delete already merged UMIs in the sorted_UMI_survived_hash_out keys
        to_delete_UMI_list = []
        for key_UMI, target_list in sorted_UMI_survived_hash_out.items():
            if len(target_list) > 0:
                to_delete_UMI_list += target_list

        for to_delete_UMI in to_delete_UMI_list:
            del sorted_UMI_survived_hash_out[to_delete_UMI]
        logging.info('Delete process 2 finished')
        logging.info('Final # of UMIs after UMI merge: %d' % len(sorted_UMI_survived_hash_out))

        return sorted_UMI_survived_hash_out

    def get_UMI_merged_dict(self, sorted_UMI_dict, UMI_survived_hash):
        UMI_merged_dict = {}

        for parent_UMI, child_UMI_list in UMI_survived_hash.items():
            UMI_merged_dict[parent_UMI] = sorted_UMI_dict[parent_UMI]
            for child_UMI in child_UMI_list:
                UMI_merged_dict[parent_UMI] += sorted_UMI_dict[child_UMI]

        return UMI_merged_dict


    ### temporary function
    def tmp_write(self, Sub_clustered, sample_dir):
        NewList = []
        for SubCluster in Sub_clustered:
            for Seq in Sub_clustered[SubCluster]:
                NewList.append({'UMI-number': SubCluster, 'seqNum': Seq['seqNum'], 'rcSum': Seq['rcSum'],
                                'readcount': Seq['readcount'], 'full_NT': Seq['full_NT'], 'rev_primer_len': Seq['rev_primer_len']})
        with open(os.path.join(sample_dir, 'SubClusteredsave_%s.tsv' % self.outfile_suffix), 'w') as handle:
            write_file(handle=handle, header=['UMI-number', 'seqNum', 'rcSum', 'readcount', 'full_NT', 'rev_primer_len'], data=NewList, delim='\t')


class ErrorCorrecter(Process):
    """ Do error correction (UMI processing)

    1a-1. primer recognition : file -> dict
    1a-2. UMI merge : dict -> dict
    1b. sub-clustering : dict -> dict
    1c. length filtering : dict -> dict
    1d. consensus sequence extraction : dict -> file
    """

    def __init__(self, in_queue, chunk_size, sub_step, out_queue, log_queue, RAM_used_queue, input_file='', input_dict={}, _dir=None,
                 fwd_primer_list=None, rev_primer_list=None, input_UMI_list=[], dist_thresh=5, hamming_thresh=1, vote_thresh=0.6, rev_primer_dist_thresh=2,
                 chain_type=None):
        Process.__init__(self)

        self.input_file = input_file
        self.in_queue = in_queue
        self.chunk_size = chunk_size
        self.sub_step = sub_step

        self.out_chunk = {}

        self.out_queue = out_queue
        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.input_dict = input_dict
        self.dir = _dir
        self.fwd_primer_list = fwd_primer_list
        self.rev_primer_list = rev_primer_list
        self.input_UMI_list = input_UMI_list

        self.dist_thresh = dist_thresh
        self.hamming_thresh = hamming_thresh
        self.vote_thresh = vote_thresh

        self.rev_primer_dist_thresh = rev_primer_dist_thresh

        self.chain_type = chain_type


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

            if self.sub_step == 'a-1':
                self.primer_recognition(chunk_start_pos, self.input_file)

            elif self.sub_step == 'a-2':
                self.UMI_merge(chunk_start_pos)

            elif self.sub_step == 'b':
                self.sub_cluster_UMI_data(chunk_start_pos)

            elif self.sub_step == 'c':
                self.different_length_byebye(chunk_start_pos)

            elif self.sub_step == 'd':
                self.extract_consensus_sequence(chunk_start_pos)

            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

            self.in_queue.task_done()

            # save chunk data into result queue
            self.out_queue.put(self.out_chunk)

        # logging.info("Finished run")


    # worker functions
    def primer_recognition(self, line_start, input_file_name):
        sHelper = Helper()

        # set the direction of primers into sense strand (plus/plus strand) (fwd, rev)
        forwardPrimerList = self.fwd_primer_list
        reversePrimerList = [sHelper.to_reverse_complement(reversePrimer) for reversePrimer in self.rev_primer_list]

        # add reverse primer sequences containing 1 mismatch (regular expression pattern)
        reversePrimerMatchList = []
        for reversePrimer in reversePrimerList:
            reversePrimerMatchList.append(reversePrimer)
        for reversePrimer in reversePrimerList:
            for i in range(len(reversePrimer)):
                reversePrimerMatchList.append(reversePrimer[:i] + reversePrimer[i + 1:])
                reversePrimerMatchList.append(reversePrimer[:i] + '[ATGC]' + reversePrimer[i + 1:])

        # add forward primer sequences containing 1 mismatch
        forwardPrimerMatchList = []
        for forwardPrimer in forwardPrimerList:
            forwardPrimerMatchList.append(forwardPrimer)
        for forwardPrimer in forwardPrimerList:
            for i in range(len(forwardPrimer)):
                forwardPrimerMatchList.append(forwardPrimer[:i] + forwardPrimer[i + 1:])
                forwardPrimerMatchList.append(forwardPrimer[:i] + '[ATGC]' + forwardPrimer[i + 1:])

        if self.chain_type == ChainType.MOUSE_HEAVY_ATC_MOUSE:
            forwardPrimerMatchList = []
            for forwardPrimer in forwardPrimerList:
                forwardPrimerMatchList.append(forwardPrimer[-15:])

        # set the forward / reverse region to find pattern
        fwd_primer_max_len = max([len(fwd) for fwd in forwardPrimerList])
        rev_primer_max_len = max([len(rev) for rev in reversePrimerList])
        forwardTrimNum = fwd_primer_max_len + 5

        if self.chain_type == ChainType.HUMAN_DELTA:
            reverseTrimNum = rev_primer_max_len + 14 + 16
        else:
            reverseTrimNum = rev_primer_max_len + 14 + 5

        UMISet = {}

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
                        if self.chain_type == ChainType.HUMAN_DELTA:
                            UMI = sequence[-(reverseTrimNum - reverseMatchResult.span()[1] - 11):]
                        else:
                            UMI = sequence[-(reverseTrimNum - reverseMatchResult.span()[1]):]

                        len_recog_rev_primer = reverseMatchResult.span()[1] - reverseMatchResult.span()[0]

                        if len(UMI) in [14, 15, 13]:
                            for forwardPrimerMatch in forwardPrimerMatchList:
                                forwardMatchResult = re.search(forwardPrimerMatch, sequence[:forwardTrimNum])

                                if forwardMatchResult != None:
                                    # # erase fwd primer region, not the rev primer region
                                    # newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[0])]
                                    newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[1])]
                                    success = True
                                    try:
                                        # UMISet[UMI].append({'UMI': UMI, 'full_NT': newSequence, 'readcount': readcount})
                                        UMISet[(UMI, line_start)].append({'UMI': UMI, 'full_NT': newSequence,
                                                                          'readcount': readcount, 'rev_primer_len': len_recog_rev_primer})
                                    except KeyError:
                                        # UMISet[UMI] = [{'UMI': UMI, 'full_NT': newSequence, 'readcount': readcount}]
                                        UMISet[(UMI, line_start)] = [{'UMI': UMI, 'full_NT': newSequence,
                                                                      'readcount': readcount, 'rev_primer_len': len_recog_rev_primer}]
                                    break
                                else:
                                    continue
                            break

                # try primer recognition process for reverse complement case.
                if success == False:
                    sequence = sHelper.to_reverse_complement(sequence)

                    for reversePrimerMatch in reversePrimerMatchList:
                        reverseMatchResult = re.search(reversePrimerMatch, sequence[-reverseTrimNum:])

                        if reverseMatchResult != None:
                            if self.chain_type == ChainType.HUMAN_DELTA:
                                UMI = sequence[-(reverseTrimNum - reverseMatchResult.span()[1] - 11):]
                            else:
                                UMI = sequence[-(reverseTrimNum - reverseMatchResult.span()[1]):]

                            len_recog_rev_primer = reverseMatchResult.span()[1] - reverseMatchResult.span()[0]

                            if len(UMI) in [14, 15, 13]:
                                for forwardPrimerMatch in forwardPrimerMatchList:
                                    forwardMatchResult = re.search(forwardPrimerMatch, sequence[:forwardTrimNum])

                                    if forwardMatchResult != None:
                                        # # erase fwd primer region, not the rev primer region
                                        # newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[0])]
                                        newSequence = sequence[forwardMatchResult.span()[1]:-(reverseTrimNum - reverseMatchResult.span()[1])]
                                        success = True
                                        try:
                                            # UMISet[UMI].append({'UMI': UMI, 'full_NT': newSequence, 'readcount': readcount})
                                            UMISet[(UMI, line_start)].append({'UMI': UMI, 'full_NT': newSequence,
                                                                              'readcount': readcount, 'rev_primer_len': len_recog_rev_primer})
                                        except KeyError:
                                            # UMISet[UMI] = [{'UMI': UMI, 'full_NT': newSequence, 'readcount': readcount}]
                                            UMISet[(UMI, line_start)] = [{'UMI': UMI, 'full_NT': newSequence,
                                                                          'readcount': readcount, 'rev_primer_len': len_recog_rev_primer}]
                                        break
                                    else:
                                        continue
                                break

        self.out_chunk = UMISet


    def UMI_merge(self, key_start):
        sorted_UMI_dict = deepcopy(self.input_dict)
        UMI_list = deepcopy(self.input_UMI_list)

        UMI_survived_hash_withOrder = {}

        key_start_count = 0
        accessed_key_count = 0

        for targetUMI in UMI_list:
            del sorted_UMI_dict[targetUMI]

            # go to start position
            if key_start_count < key_start:
                key_start_count += 1
                continue

            # stop working if the number of accessed key exceed the chunk_size
            accessed_key_count += 1
            if accessed_key_count > self.chunk_size:
                break

            # UMI_survived_hash[targetUMI] = []
            UMI_survived_hash_withOrder[(targetUMI, key_start)] = []

            # make target UMI hash
            if targetUMI == None:
                logging.info('Wrong target UMI...')
                raise ValueError
            target_UMI_hash_dict = {}
            woNTSet = {'A': ['T', 'G', 'C'], 'T': ['A', 'G', 'C'], 'G': ['A', 'T', 'C'], 'C': ['A', 'T', 'G']}
            for i in range(len(targetUMI)):
                target_UMI_hash_dict[targetUMI[:i] + targetUMI[i + 1:]] = 0         # for 1bp deleted cases (?)
                for woNT in woNTSet[targetUMI[i]]:
                    target_UMI_hash_dict[targetUMI[:i] + woNT + targetUMI[i + 1:]] = 0
            target_UMI_hash = list(target_UMI_hash_dict.keys())

            # do UMI merge: 1 mismatch is allowed and merged.
            for UMI_hash in target_UMI_hash:
                try:
                    sorted_UMI_dict[UMI_hash]
                    UMI_survived_hash_withOrder[(targetUMI, key_start)].append(UMI_hash)
                except KeyError:
                    pass

        self.out_chunk = UMI_survived_hash_withOrder


    def sub_cluster_UMI_data(self, key_start):
        mergedUMISet = self.input_dict
        ClusteredUMISet = {}

        key_start_count = 0
        accessed_key_count = 0
        for mergedUMI in mergedUMISet:
            # go to start position
            if key_start_count < key_start:
                key_start_count += 1
                continue

            # stop working if the number of accessed key exceed the chunk_size
            accessed_key_count += 1
            if accessed_key_count > self.chunk_size:
                break

            UMI = mergedUMI

            # sorting sequences within the UMI group by readcount sum
            readcountSet = {}
            for mergedUMInum in mergedUMISet[mergedUMI]:
                temSeq = {'full_NT': mergedUMInum['full_NT'], 'readcount': int(mergedUMInum['readcount']),
                          'rev_primer_len': int(mergedUMInum['rev_primer_len'])}
                try:
                    readcountSet[str(mergedUMInum['readcount'])].append(temSeq)
                except KeyError:
                    readcountSet[str(mergedUMInum['readcount'])] = [temSeq]

            readcountList = []
            for readcount in readcountSet:
                readcountList.append(int(readcount))

            readcountList = list(set(readcountList))
            readcountList.sort(reverse=True)

            seqList = []
            for readcount in readcountList:
                seqList += readcountSet[str(readcount)]

            # start sub-clustering sequences (full_NT) belonging to the UMI group
            subClusterNum = 1
            while (len(seqList) > 0):
                targetSeq = seqList.pop(0)
                sbjct_trim_3_prime = int(targetSeq['rev_primer_len'])
                sbjct_rev_primer = targetSeq['full_NT'][-sbjct_trim_3_prime:]

                resultSeqList = [targetSeq]
                readcountSum = targetSeq['readcount']

                popIndexList = []

                for i, seq in enumerate(seqList):
                    # reverse primer similarity check
                    query_trim_3_prime = int(seq['rev_primer_len'])
                    query_rev_primer = seq['full_NT'][-query_trim_3_prime:]
                    dist_of_rev_primers = Helper.levenshtein_distance(sbjct_rev_primer, query_rev_primer)
                    if dist_of_rev_primers > self.rev_primer_dist_thresh:
                        continue

                    # sequence similarity check
                    distFromTarget = Helper.levenshtein_distance(targetSeq['full_NT'][:-sbjct_trim_3_prime], seq['full_NT'][:-query_trim_3_prime])
                    if distFromTarget <= self.dist_thresh:
                        resultSeqList.append(seq)
                        readcountSum += seq['readcount']
                        popIndexList.append(i)

                popIndexList.sort(reverse=True)

                # erase clustered sequences from seqList
                for popIndex in popIndexList:
                    seqList.pop(popIndex)

                if readcountSum == 2 and len(resultSeqList) == 2:
                    continue
                elif readcountSum > 1:
                    for resultSeq in resultSeqList:
                        try:
                            ClusteredUMISet[UMI + '-' + str(subClusterNum)].append({'seqNum': len(resultSeqList), 'rcSum': readcountSum,
                                                                                    'readcount': resultSeq['readcount'], 'full_NT': resultSeq['full_NT'],
                                                                                    'rev_primer_len': resultSeq['rev_primer_len']})
                        except KeyError:
                            ClusteredUMISet[UMI + '-' + str(subClusterNum)] = [{'seqNum': len(resultSeqList), 'rcSum': readcountSum,
                                                                                'readcount': resultSeq['readcount'], 'full_NT': resultSeq['full_NT'],
                                                                                'rev_primer_len': resultSeq['rev_primer_len']}]

                    subClusterNum += 1

        self.out_chunk = ClusteredUMISet


    def different_length_byebye(self, key_start):
        ClusteredUMISet = self.input_dict

        DiffLengthBye = {}
        shorterhamming = []
        DiffLengthNoBye = {}
        howmany = 0
        whatdiff = {}

        key_start_count = 0
        accessed_key_count = 0
        for ClusteredUMI in ClusteredUMISet:
            # go to start position
            if key_start_count < key_start:
                key_start_count += 1
                continue

            # stop working if the number of accessed key exceed the chunk_size
            accessed_key_count += 1
            if accessed_key_count > self.chunk_size:
                break

            lengthandnum = {}

            for Seq in ClusteredUMISet[ClusteredUMI]:
                try:
                    lengthandnum[len(Seq['full_NT'])] = lengthandnum[len(Seq['full_NT'])] + Seq['readcount']
                except KeyError:
                    lengthandnum[len(Seq['full_NT'])] = Seq['readcount']

            max_num = 0
            dominantlength = 0
            # targetList_arranged = sorted(targetList, key=lambda x: x['readcount'], reverse=True)

            lengthandnumm = {k: v for k, v in (sorted(lengthandnum.items(), reverse=True, key=lambda item: item[1]))}
            if len(lengthandnumm) >= 2:
                if list(lengthandnumm.values())[0] == list(lengthandnumm.values())[1]:
                    if len(lengthandnumm) == 2:
                        DiffLengthAlmostBye = {}
                        howmany += 1
                        fr_trimmed = {}
                        forward_trimmed = {}
                        reverse_trimmed = {}
                        small_length = {}
                        max_length = 0
                        min_length = 0

                        for q in lengthandnumm:
                            if int(q) >= max_length:
                                max_length = int(q)
                        for q in lengthandnumm:
                            if int(q) != max_length:
                                min_length = int(q)
                        for Seq in ClusteredUMISet[ClusteredUMI]:
                            if len(Seq['full_NT']) == max_length:
                                try:
                                    fr_trimmed[ClusteredUMI].append(Seq)
                                except KeyError:
                                    fr_trimmed[ClusteredUMI] = [Seq]

                            else:
                                try:
                                    small_length[ClusteredUMI].append(Seq)
                                except KeyError:
                                    small_length[ClusteredUMI] = [Seq]

                        for Seq in fr_trimmed[ClusteredUMI]:
                            try:
                                forward_trimmed[ClusteredUMI].append(
                                    {'rcSum': Seq['rcSum'], 'seqNum': Seq['rcSum'], 'readcount': Seq['readcount'], 'full_NT': Seq['full_NT'][max_length - min_length:]})
                            except KeyError:
                                forward_trimmed[ClusteredUMI] = [
                                    {'rcSum': Seq['rcSum'], 'seqNum': Seq['rcSum'], 'readcount': Seq['readcount'], 'full_NT': Seq['full_NT'][max_length - min_length:]}]

                        for Seq in fr_trimmed[ClusteredUMI]:
                            try:
                                reverse_trimmed[ClusteredUMI].append(
                                    {'rcSum': Seq['rcSum'], 'seqNum': Seq['rcSum'], 'readcount': Seq['readcount'],
                                     'full_NT': Seq['full_NT'][:min_length - max_length]})
                            except KeyError:
                                reverse_trimmed[ClusteredUMI] = [
                                    {'rcSum': Seq['rcSum'], 'seqNum': Seq['rcSum'], 'readcount': Seq['readcount'],
                                     'full_NT': Seq['full_NT'][:min_length - max_length]}]

                        forward_hamming = []
                        reverse_hamming = []
                        for Seq in forward_trimmed[ClusteredUMI]:
                            hamOneFAllShortList = []
                            for small_seq in small_length[ClusteredUMI]:
                                hamOneFAllShortList.append(Helper.hamming_distance(small_seq['full_NT'], Seq['full_NT']))
                            min_ham = hamOneFAllShortList[0]
                            for hamOneFAllShort in hamOneFAllShortList:
                                if hamOneFAllShort <= min_ham:
                                    min_ham = hamOneFAllShort
                            forward_hamming.append(min_ham)
                        for Seq in reverse_trimmed[ClusteredUMI]:
                            hamOneRAllShortList = []
                            for small_seq in small_length[ClusteredUMI]:
                                hamOneRAllShortList.append(Helper.hamming_distance(small_seq['full_NT'], Seq['full_NT']))
                            min_ham = hamOneRAllShortList[0]
                            for hamOneRAllShort in hamOneRAllShortList:
                                if hamOneRAllShort <= min_ham:
                                    min_ham = hamOneRAllShort
                            reverse_hamming.append(min_ham)
                        for a in range(len(forward_hamming)):
                            # print(len(forward_hamming))
                            if forward_hamming[a] <= reverse_hamming[a] and forward_hamming[a] <= self.hamming_thresh:
                                try:
                                    DiffLengthAlmostBye[ClusteredUMI].append(forward_trimmed[ClusteredUMI][a])
                                    # print(DiffLengthBye[ClusteredUMI])
                                except KeyError:
                                    DiffLengthAlmostBye[ClusteredUMI] = small_length[ClusteredUMI]
                                    DiffLengthAlmostBye[ClusteredUMI].append(forward_trimmed[ClusteredUMI][a])
                            elif reverse_hamming[a] <= forward_hamming[a] and reverse_hamming[a] <= self.hamming_thresh:
                                try:
                                    DiffLengthAlmostBye[ClusteredUMI].append(reverse_trimmed[ClusteredUMI][a])
                                except KeyError:
                                    DiffLengthAlmostBye[ClusteredUMI] = small_length[ClusteredUMI]
                                    DiffLengthAlmostBye[ClusteredUMI].append(reverse_trimmed[ClusteredUMI][a])
                        try:
                            if len(DiffLengthAlmostBye[ClusteredUMI]) > 0.6 * len(ClusteredUMISet[ClusteredUMI]):
                                DiffLengthBye[ClusteredUMI] = DiffLengthAlmostBye[ClusteredUMI]
                        except KeyError:
                            w = 0

                else:
                    DiffAlmostSet = {}
                    for kk in lengthandnumm:
                        if lengthandnumm[kk] >= max_num:
                            max_num = lengthandnumm[kk]
                            dominantlength = kk

                    for Seqq in ClusteredUMISet[ClusteredUMI]:
                        if len(Seqq['full_NT']) == dominantlength:
                            try:
                                DiffAlmostSet[ClusteredUMI].append(Seqq)
                            except KeyError:
                                DiffAlmostSet[ClusteredUMI] = [Seqq]
                    try:
                        if len(DiffAlmostSet[ClusteredUMI]) > 0.6 * len(ClusteredUMISet[ClusteredUMI]):
                            DiffLengthBye[ClusteredUMI] = DiffAlmostSet[ClusteredUMI]
                    except KeyError:
                        qqq = 0
            else:
                DiffLengthBye[ClusteredUMI] = ClusteredUMISet[ClusteredUMI]

        for UMI in DiffLengthBye:
            readcountSum = 0
            for each in DiffLengthBye[UMI]:
                readcountSum += int(each['readcount'])
            for each in DiffLengthBye[UMI]:
                each['rcSum'] = readcountSum

        self.out_chunk = DiffLengthBye


    def extract_consensus_sequence(self, key_start):
        ClusteredUMISet = self.input_dict

        resultSeqSet = {}
        key_start_count = 0
        accessed_key_count = 0
        for ClusteredUMI in ClusteredUMISet:
            # go to start position
            if key_start_count < key_start:
                key_start_count += 1
                continue

            # stop working if the number of accessed key exceed the chunk_size
            accessed_key_count += 1
            if accessed_key_count > self.chunk_size:
                break

            seqNum = int(ClusteredUMISet[ClusteredUMI][0]['seqNum'])
            rcSum = int(ClusteredUMISet[ClusteredUMI][0]['rcSum'])

            if rcSum == 1:
                continue
            if seqNum == 2 and rcSum == 2:
                continue

            seqList = ClusteredUMISet[ClusteredUMI]

            seqScoreList = []

            # get nucleotide counts (considering readcount of the reads) by the position. must have the same length within sub-cluster.
            for i, seq in enumerate(seqList):
                if i == 0:
                    for eachSeq in seq['full_NT']:
                        tempSet = {'A': 0, 'G': 0, 'T': 0, 'C': 0, '-': 0}
                        tempSet[eachSeq] += int(seq['readcount'])
                        seqScoreList.append(tempSet)
                else:
                    for j, eachSeq in enumerate(seq['full_NT']):
                        seqScoreList[j][eachSeq] += int(seq['readcount'])

            # get consensus sequence if maxVote nt occupies over the vote_thresh.
            resultSeq = ''
            for seqScore in seqScoreList:
                maxVote = max(seqScore.values())
                if float(maxVote) / rcSum > self.vote_thresh:
                    maxVoteNT = list(seqScore.keys())[list(seqScore.values()).index(maxVote)]
                    resultSeq += maxVoteNT
                else:
                    resultSeq += 'X'

            try:
                resultSeq.index('X')
            except ValueError:
                try:
                    # resultSeqSet[resultSeq.replace('-', '')] += 1
                    resultSeqSet[resultSeq.replace('-', '')]['readcount'] += 1
                    resultSeqSet[resultSeq.replace('-', '')]['UMI_list'].append(ClusteredUMI.split('-')[0])

                except KeyError:
                    # resultSeqSet[resultSeq.replace('-', '')] = 1
                    resultSeqSet[resultSeq.replace('-', '')] = {'full_NT': resultSeq.replace('-', ''), 'readcount': 1, 'UMI_list': [ClusteredUMI.split('-')[0]]}

        self.out_chunk = resultSeqSet


    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue


class IsotypeExtractor(Process):
    """ Do c gene alignment and isotype annotation

    2a. run blast : file -> file
    2b. parse blast run result and get isotype info : file -> dict
    """

    def __init__(self, input_file, sub_step, log_queue, RAM_used_queue, out_dir,
                 c_ref_file=None, in_queue=None, out_queue=None, num_process=1, chunk_size=None, delim='\t', chain_type='', subtyping=True):
        Process.__init__(self)

        self.input_file = input_file
        self.sub_step = sub_step

        self.out_chunk = {}

        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.out_dir = out_dir
        self.c_ref_file = c_ref_file
        self.in_queue = in_queue
        self.out_queue = out_queue

        self.num_process = num_process
        self.chunk_size = chunk_size

        self.delim = delim

        self.chain_type = chain_type
        self.subtyping = subtyping


    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        # logging.info("Started run")
        if self.sub_step == 'a':
            self.run_blast(self.input_file)

            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

        elif self.sub_step == 'b':
            ## load data into RAM: error-corrected files for loading sequence & readcount info
            targetHeader, targetData = read_file_new(self.input_file, delim=self.delim)

            while True:
                # Get the work from the queue and expand the tuple
                chunk_start_pos = self.in_queue.get()
                if chunk_start_pos is None:
                    self.in_queue.task_done()
                    break

                # extract isotype hash table
                self.parse_isotype(chunk_start_pos, targetData)
                self.in_queue.task_done()

                # RAM usage check
                tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
                self.RAM_used_queue.put(tmp_RAM_used)

                # save chunk data into result queue
                self.out_queue.put(self.out_chunk)

        # logging.info("Finished run")


    # worker functions
    def run_blast(self, input_file, **kwargs):
        chRefHeder, chRefList = read_file_new(in_file=self.c_ref_file, delim=',')
        targetSetHeader, targetSeqList = read_file_new(input_file, delim=self.delim)

        # directory setting
        if not self.out_dir:
            output_dir = '/'.join(os.path.dirname(input_file).split('/')[:-1] + ['2_annotation'])
        else:
            output_dir = self.out_dir
        blast_dir = output_dir + 'temp_cblast'
        query_dir = os.path.join(blast_dir, 'query')
        db_dir = os.path.join(blast_dir, 'db')
        result_dir = os.path.join(blast_dir, 'result')

        # make directories
        t_dirs = [blast_dir, query_dir, db_dir, result_dir]
        for t_dir in t_dirs:
            if not os.path.exists(t_dir):
                os.mkdir(t_dir)

        # set parameters
        if True:
            if 'evalue' in kwargs:
                evalue = kwargs['evalue']
            else:
                evalue = 0.1 ** 1

            if 'gapopen' in kwargs:
                gapopen = kwargs['gapopen']
            else:
                gapopen = 1

            if 'gapextend' in kwargs:
                gapextend = kwargs['gapextend']
            else:
                gapextend = 2

            if 'word_size' in kwargs:
                word_size = kwargs['word_size']
            else:
                word_size = 10

            if 'num_alignments' in kwargs:
                num_alignments = kwargs['num_alignments']
            else:
                num_alignments = 1

            query_list = [x['full_NT'] for x in targetSeqList]
            db_dict = {x[chRefHeder[0]]:x[chRefHeder[1]] for x in chRefList}
            result_file_name = 'blasted_UMI.txt'

        result_file = os.path.join(result_dir, result_file_name)

        with open(os.path.join(query_dir, 'query_UMI.fasta'), 'w') as fasta_writer:
            for i, e_seq in enumerate(query_list):
                fasta_writer.write('>' + str(i) + '\n')
                fasta_writer.write(e_seq + '\n')

        with open(os.path.join(db_dir, 'db_UMI.fasta'), 'w') as fasta_writer:
            for name, e_seq in db_dict.items():
                fasta_writer.write('>' + name + '\n')
                fasta_writer.write(e_seq + '\n')

        # construct db
        format_cmd = '%s -in %s -dbtype nucl -input_type fasta -out %s' % (
            '/Tools/ncbi-blast-2.7.1+/bin/makeblastdb', os.path.join(db_dir, 'db_UMI.fasta'), os.path.join(db_dir, 'db_UMI'))
        os.system(format_cmd)

        # run blastn
        cline = NcbiblastnCommandline(num_threads=self.num_process, query=os.path.join(query_dir, 'query_UMI.fasta'), db=os.path.join(db_dir, 'db_UMI'),
                                      evalue=evalue, out=result_file, gapopen=gapopen, gapextend=gapextend, word_size=word_size,
                                      num_descriptions=num_alignments, num_alignments=num_alignments)
        format_exe = '/Tools/ncbi-blast-2.7.1+/bin/blastn'
        os.system(format_exe + str(cline)[len('blastn'):])

    # make user defined blast parser for the time efficiency -> do not use NCBIStandalone module
    def parse_isotype(self, key_start, input_data_list):
        result_dict = {}
        blast_result_file = os.path.join(self.out_dir, 'blasted_UMI.txt')

        with open(blast_result_file, 'r') as handle:
            key_start_count = 0
            accessed_key_count = 0

            while True:
                info_accessed = 0
                line = handle.readline()
                if line == '':
                    break
                if line[:len('Query=')] == 'Query=':
                    # go to start position
                    if key_start_count < key_start:
                        key_start_count += 1
                        continue

                    # stop working if the number of accessed key exceed the chunk_size
                    accessed_key_count += 1
                    if accessed_key_count > self.chunk_size:
                        break

                    access_num = int(line.split()[-1].strip())
                    seq = input_data_list[access_num]['full_NT']
                    rc = input_data_list[access_num]['readcount']
                    umi_list = input_data_list[access_num]['umi_list']

                    # get query length
                    handle.readline()
                    line = handle.readline()
                    if line[:len('Length=')] == 'Length=':
                        query_len = int(line.split('=')[-1].strip())

                    query_start = 0
                    query_end = 0
                    sbjct_start = 0

                    # user defined blast parser (specific purpose: extract required info for isotyping)
                    while True:
                        line = handle.readline()
                        if 'No hits found' in line:
                            break
                        elif line[:len('>')] == '>':
                            sbjct_name = line.split()[-1].strip()
                            info_accessed += 1
                        elif line.strip()[:len('Identities')] == 'Identities':
                            identities = line.split(',')[0].split('=')[-1].strip().split()[0]
                            matched, whole = identities.split('/')
                            info_accessed += 2
                        elif line.strip()[:len('Strand')] == 'Strand':
                            direction = line.split('=')[-1].strip()
                            if direction != 'Plus/Plus':
                                break
                            info_accessed += 1
                        elif line[:len('Query')] == 'Query':
                            query_end = int(line.split()[-1].strip())
                            if query_start > 0:
                                continue
                            else:
                                query_start = int(line.split()[1].strip())
                                info_accessed += 1
                        elif line[:len('Sbjct')] == 'Sbjct':
                            if sbjct_start > 0:
                                continue
                            sbjct_start = int(line.split()[1].strip())
                            info_accessed += 1

                        # end of alignment
                        elif line[:len('Effective search space used:')] == 'Effective search space used:':
                            if query_end > 0:
                                info_accessed += 1
                            break
                        # if info_accessed == 6:
                        #     break

                    if info_accessed < 7:
                        continue

                    # pass not well aligned cases
                    numMismatches = int(whole) - int(matched)
                    # if int(whole) < 45 or int(whole) > 71 or int(numMismatches) > 5 or int(sbjct_start) > 10 or int(query_start) > 400:
                    #     continue

                    if int(numMismatches) > 5 or int(sbjct_start) > 10 or (int(query_len) - int(query_end)) > 10:
                        continue

                    # get isotype information and save
                    vdj_seq = seq[:query_start - sbjct_start]
                    if self.chain_type in [ChainType.HUMAN_HEAVY, ChainType.MOUSE_HEAVY]:
                        if self.subtyping:
                            isotype = sbjct_name.split('*')[0][3:]
                        else:
                            isotype = sbjct_name.split('*')[0][3]
                    else:
                        isotype = sbjct_name.split('*')[0][:4]
                    new_key = '|'.join([vdj_seq, isotype])
                    try:
                        result_dict[new_key]['readcount'] += rc
                    except KeyError:
                        result_dict[new_key] = {'full_NT': vdj_seq, 'umi_list': umi_list, 'isotype': isotype, 'readcount': rc}

            self.out_chunk = result_dict

    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue


class InfoAnnotator(Process):
    """ Do IgBLAST run and region annotation

    3a. run igblast : file -> file
    3b. parse igblast run result and get region info : file -> dict
    """

    def __init__(self, input_file, sub_step, log_queue, RAM_used_queue, out_dir, chain_type,
                 in_queue=None, out_queue=None, num_process=1, chunk_size=None, delim='\t'):
        Process.__init__(self)

        self.input_file = input_file
        self.sub_step = sub_step

        self.out_chunk = {}

        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.out_dir = out_dir
        self.chain_type = chain_type
        self.in_queue = in_queue
        self.out_queue = out_queue

        self.num_process = num_process
        self.chunk_size = chunk_size

        self.delim = delim


    def run(self):
        if self.log_queue:
            logger = logging.getLogger()
            if len(logger.handlers) == 0:
                h = QueueHandler(self.log_queue)
                logger.addHandler(h)
                logger.setLevel(logging.DEBUG)

        # logging.info("Started run")
        if self.sub_step == 'a':
            self.run_igblast()

            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

        elif self.sub_step == 'b':
            ## load data into RAM: error-corrected files for loading sequence & readcount info
            targetHeader, targetData = read_file_new(self.input_file, delim=self.delim)

            # # some parameters...
            # self.domain_system = 'imgt'

            while True:
                # Get the work from the queue and expand the tuple
                chunk_start_pos = self.in_queue.get()
                if chunk_start_pos is None:
                    self.in_queue.task_done()
                    break

                # extract isotype hash table
                self.annotate_info(chunk_start_pos, targetData)
                self.in_queue.task_done()

                # RAM usage check
                tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
                self.RAM_used_queue.put(tmp_RAM_used)

                # save chunk data into result queue
                self.out_queue.put(self.out_chunk)

        # logging.info("Finished run")


    # worker functions → not use Helper class in SeqUtil.py. separate run_igblast and igblast_parser
    def run_igblast(self):
        ### parameter setting
        targetHeader, targetList = read_file_new(self.input_file, delim=self.delim)
        sequence_list = [d['full_NT'] for d in targetList]

        blast_dir = os.path.join(os.path.join(self.out_dir, 'temp_igblast'))
        if not os.path.exists(blast_dir):
            os.mkdir(blast_dir)

        igblast_name = os.path.basename(self.input_file)[:-len('.tsv')]
        file_blast = os.path.join(blast_dir, igblast_name + '_igblasted_UMI.txt')

        if self.chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT,
                               ChainType.CHICKEN_HEAVY, ChainType.CHICKEN_LIGHT,
                               ChainType.RAT_HEAVY, ChainType.RAT_LIGHT,
                               ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA,
                               ChainType.MOUSE_BALBC_HEAVY, ChainType.MOUSE_BALBC_LIGHT,
                               ChainType.MOUSE_C57BL6_HEAVY, ChainType.MOUSE_C57BL6_LIGHT,
                               ChainType.MOUSE_HEAVY, ChainType.MOUSE_LIGHT,
                               ChainType.ALPACA_HEAVY]:
            seq_type = 'Ig'
            # domain_system = 'kabat'
            domain_system = 'imgt'
        elif self.chain_type in [ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_DELTA, ChainType.HUMAN_GAMMA,
                                 ChainType.MOUSE_BETA, ChainType.MOUSE_ALPHA, ChainType.MOUSE_GAMMA, ChainType.MOUSE_DELTA]:
            seq_type = 'TCR'
            domain_system = 'imgt'
        else:
            raise RuntimeError('chain type error')

        # check sequence type
        is_nucl = is_nucleotide_or_protein(sequence_list[0])
        if is_nucl == None or is_nucl == False:
            raise RuntimeError("parameter sequence_list is not a nucleotide sequence list")

        # make temporal query fasta file
        tp_query = tempfile.NamedTemporaryFile('wt', suffix='.fasta', delete=False)
        for i, seq in enumerate(sequence_list):
            tp_query.write(">%d\n" % i)
            tp_query.write(seq + '\n')
        tp_query.close()

        # run igblast
        run_igblast_new(chain_type=self.chain_type, query=tp_query.name, out=file_blast,
                        seq_type=seq_type, domain_system=domain_system, num_threads=self.num_process)

    def annotate_info(self, key_start, input_data_list):
        sHelper = Helper()
        igblast_name = os.path.basename(self.input_file)[:-len('.tsv')]
        igblast_result_file = os.path.join(self.out_dir, '%s_igblasted_UMI.txt' % igblast_name)

        ### parse igblast result txt file.
        igblast_parsed_dict = {}
        # igblast_parsed_list = []
        with open(igblast_result_file, 'r') as igb_handle:
            key_start_count = 0
            accessed_key_count = 0

            # for the case that domain system is "imgt"
            # include BCR & TCR
            if self.chain_type in [ChainType.HUMAN_HEAVY, ChainType.HUMAN_LIGHT,
                                   ChainType.RAT_HEAVY, ChainType.RAT_LIGHT,
                                   ChainType.RABBIT_HEAVY, ChainType.RABBIT_KAPPA,
                                   ChainType.MOUSE_BALBC_HEAVY, ChainType.MOUSE_BALBC_LIGHT,
                                   ChainType.MOUSE_C57BL6_HEAVY, ChainType.MOUSE_C57BL6_LIGHT,
                                   ChainType.ALPACA_HEAVY,
                                   ChainType.MOUSE_HEAVY, ChainType.MOUSE_LIGHT,
                                   ChainType.HUMAN_BETA, ChainType.HUMAN_ALPHA, ChainType.HUMAN_DELTA, ChainType.HUMAN_GAMMA,
                                   ChainType.MOUSE_BETA, ChainType.MOUSE_ALPHA, ChainType.MOUSE_DELTA, ChainType.MOUSE_GAMMA]:
                # for the case that domain system is "imgt"
                while True:
                    line = igb_handle.readline()
                    if not line:
                        break

                    if line[:len('Query= ')] == 'Query= ':
                        # go to start position
                        if key_start_count < key_start:
                            key_start_count += 1
                            continue

                        # stop working if the number of accessed key exceed the chunk_size
                        accessed_key_count += 1
                        if accessed_key_count > self.chunk_size:
                            break

                        row_num = int(line[len('Query= '):].strip())

                        junction = ''
                        junction_aa = ''
                        junction_start = 0
                        junction_end = 0
                        CDR3_NT = ''
                        sequence_alignment = ''
                        seq_aa = ''

                        result = {'query': line[len('Query= '):-1], 'hit': False, 'v_call': 'N/A', 'j_call': 'N/A',
                                  'cdr3': 'N/A', 'cdr3_aa': 'X',
                                  'junction': 'N/A', 'junction_aa': 'X', 'rev_comp': 'F', 'productive': 'F',
                                  'v_cigar': '', 'j_cigar': '', 'v_alignment_length': 'N/A', 'v_alignment_mutation': 'N/A',
                                  'sequence_alignment': 'N/A', 'germline_alignment': 'N/A'}
                        if self.chain_type in [ChainType.HUMAN_HEAVY, ChainType.RABBIT_HEAVY,
                                               ChainType.MOUSE_C57BL6_HEAVY, ChainType.MOUSE_BALBC_HEAVY, ChainType.MOUSE_HEAVY,
                                               ChainType.HUMAN_BETA, ChainType.HUMAN_DELTA,
                                               ChainType.MOUSE_BETA, ChainType.MOUSE_DELTA]:
                            result['d_call'] = 'N/A'
                            result['d_cigar'] = ''

                        while True:
                            line = igb_handle.readline()

                            if 'No hits found' in line:
                                result['hit'] = False
                                igblast_parsed_dict[row_num] = result
                                break

                            elif 'Note that your query represents the minus strand' in line:
                                result['rev_comp'] = 'T'

                            elif line[:len('Sub-region sequence details')] == 'Sub-region sequence details':
                                next_line_sp = igb_handle.readline().split('\t')
                                if next_line_sp[0] == 'CDR3':
                                    result['hit'] = True
                                    CDR3_NT = next_line_sp[1]
                                    result['cdr3'] = CDR3_NT
                                    result['cdr3_aa'] = sHelper.translate(result['cdr3'])
                                    result['cdr3_start'] = int(next_line_sp[3])
                                    junction_start = int(next_line_sp[3]) - 4
                                    if junction_start != 0:
                                        result['junction_start'] = junction_start
                                    result['cdr3_end'] = int(next_line_sp[4])
                                    junction_end = int(next_line_sp[4]) + 3
                                    if junction_end != 0:
                                        result['junction_end'] = junction_end
                                    result['fr4_from'] = int(next_line_sp[4])

                            elif line[:len('Alignments')] == 'Alignments':
                                next_line_sp = igb_handle.readline()
                                next_line_sp = igb_handle.readline()
                                next_line_sp = igb_handle.readline().split(' ')
                                next_line_trimmed = []
                                for l in range(len(next_line_sp)):
                                    if next_line_sp[l] != '':
                                        next_line_trimmed.append(next_line_sp[l])
                                query_line = []
                                query_start_where = []
                                query_finish_where = []
                                V_germline_mismatch = []
                                D_germline_mismatch = []
                                J_germline_mismatch = []
                                V_mis_ratio = []
                                D_mis_ratio = []
                                J_mis_ratio = []
                                v_germ_start_where = []
                                v_germ_finish_where = []
                                d_germ_start_where = []
                                d_germ_finish_where = []
                                j_germ_start_where = []
                                j_germ_finish_where = []
                                sequence_a = []
                                dog = True
                                while dog:
                                    next_line_trimmed = []
                                    before_line_sp = next_line_sp
                                    next_line_sp = igb_handle.readline().split(' ')
                                    for l in range(len(next_line_sp)):
                                        if next_line_sp[l] != '':
                                            next_line_trimmed.append(next_line_sp[l])

                                    before_line_trimmed = []
                                    for l in range(len(before_line_sp)):
                                        if before_line_sp[l] != '':
                                            before_line_trimmed.append(before_line_sp[l])

                                    findQu = re.compile('Query_')
                                    findQue = findQu.search(next_line_trimmed[0])

                                    findper = re.compile('%')
                                    if len(next_line_trimmed) >= 2:
                                        findperc = findper.search(next_line_trimmed[1])

                                    if findQue:
                                        query_start_where.append(next_line_trimmed[1])
                                        query_line.append(next_line_trimmed[2])
                                        query_finish_where.append(next_line_trimmed[3].replace('\n', ''))
                                        seq_a = ''
                                        for l in range(len(before_line_trimmed)):
                                            seq_a = seq_a + before_line_trimmed[l]
                                        sequence_a.append(seq_a.replace('\n', ''))

                                    elif next_line_trimmed[0] == 'V' and findperc:
                                        if len(next_line_trimmed) >= 7:
                                            if len(V_mis_ratio) == 0:
                                                V_mis_ratio.append(next_line_trimmed[1:3])
                                            v_germ_start_where.append(next_line_trimmed[4])
                                            V_germline_mismatch.append(next_line_trimmed[5])
                                            v_germ_finish_where.append(next_line_trimmed[6].replace('\n', ''))

                                    elif next_line_trimmed[0] == 'D' and findperc:
                                        if len(next_line_trimmed) >= 7:
                                            if len(D_mis_ratio) == 0:
                                                D_mis_ratio.append(next_line_trimmed[1:3])
                                            d_germ_start_where.append(next_line_trimmed[4])
                                            D_germline_mismatch.append(next_line_trimmed[5])
                                            d_germ_finish_where.append(next_line_trimmed[6].replace('\n', ''))

                                    elif next_line_trimmed[0] == 'J' and findperc:
                                        if len(next_line_trimmed) >= 7:
                                            if len(J_mis_ratio) == 0:
                                                J_mis_ratio.append(next_line_trimmed[1:3])
                                            j_germ_start_where.append(next_line_trimmed[4])
                                            J_germline_mismatch.append(next_line_trimmed[5])
                                            j_germ_finish_where.append(next_line_trimmed[6].replace('\n', ''))

                                    elif next_line_trimmed[0] == 'Lambda':
                                        dog = False
                                seq_aa = ''
                                for k in range(len(sequence_a)):
                                    seq_aa = seq_aa + sequence_a[k]
                                result['sequence_aa'] = seq_aa
                                result['query_line'] = query_line
                                for a in range(len(query_line)):
                                    sequence_alignment = sequence_alignment + query_line[a]
                                result['sequence_alignment'] = sequence_alignment
                                seq_without = ''
                                if len(sequence_alignment) != 0:
                                    for z in range(len(sequence_alignment)):
                                        if sequence_alignment[z] != '-':
                                            seq_without = seq_without + sequence_alignment[z]
                                if len(seq_without) != 0 and junction_start != 0:
                                    junction = seq_without[junction_start - int(query_start_where[0]) + 1:junction_end - int(query_start_where[0]) + 1]
                                if CDR3_NT != '':
                                    result['junction'] = junction
                                    result['junction_aa'] = sHelper.translate(junction)
                                else:
                                    result['junction'] = ''
                                    result['junction_aa'] = ''
                                v_finish_from_zero = 0
                                if len(v_germ_start_where) != 0:
                                    v_finish_from_zero = int(v_germ_finish_where[len(v_germ_finish_where) - 1]) - int(v_germ_start_where[0])
                                v_cigar_mismatch = []
                                v_cigar_mismatches = []
                                j_cigar_mismatches = []
                                d_cigar_mismatches = []
                                if (len(sequence_alignment) >= v_finish_from_zero + 1) and v_finish_from_zero != 0:
                                    for a in range(v_finish_from_zero + 1):
                                        if sequence_alignment[a] == '-':
                                            v_cigar_mismatch.append([a, '0'])
                                if len(v_germ_start_where) != 0:
                                    for a in range(len(v_germ_start_where)):
                                        for i in range(int(v_germ_finish_where[a]) - int(v_germ_start_where[a]) + 1):
                                            if V_germline_mismatch[a][i] != '.':
                                                v_cigar_mismatch.append([i + int(v_germ_start_where[a]) - int(v_germ_start_where[0]), V_germline_mismatch[a][i]])
                                v_cigar_mismatch.sort(key=lambda x: x[0])
                                for kk in range(len(v_cigar_mismatch)):
                                    isit = 0
                                    for ll in range(len(v_cigar_mismatches)):
                                        if v_cigar_mismatch[kk][0] == v_cigar_mismatches[ll][0]:
                                            isit += 1
                                    if isit == 0:
                                        v_cigar_mismatches.append(v_cigar_mismatch[kk])

                                d_finish_from_zero = 0
                                if len(d_germ_start_where) != 0:
                                    d_finish_from_zero = int(d_germ_finish_where[len(d_germ_finish_where) - 1]) - int(
                                        d_germ_start_where[0])

                                d_cigar_mismatch = []
                                d_line_pos = []
                                if (len(d_germ_start_where)) != 0:
                                    for k in range(90):
                                        if D_germline_mismatch[0][k] != '-':
                                            d_line_po = k
                                            d_line_pos.append(k)
                                            break

                                seq_before_d_num = 0
                                if d_finish_from_zero != 0 and len(v_germ_finish_where) != 0 and len(sequence_alignment) != 0:
                                    for q in range(len(v_germ_start_where) - 1):
                                        # seq_before_d_num : 90, 180, ...
                                        seq_before_d_num = seq_before_d_num + 90
                                    for a in range(len(d_cigar_mismatch)):
                                        if sequence_alignment[seq_before_d_num + d_line_pos[0] + 90 * a] == '-':
                                            d_cigar_mismatch.append([seq_before_d_num + d_line_pos[0] + 90 * a, '0'])
                                if len(d_germ_start_where) == 1:
                                    for i in range(int(d_germ_finish_where[0]) - int(d_germ_start_where[0])):
                                        if D_germline_mismatch[0][d_line_pos[0] + i] != '.':
                                            d_cigar_mismatch.append(
                                                [i + seq_before_d_num + d_line_pos[0], D_germline_mismatch[0][d_line_pos[0] + i]])
                                elif len(d_germ_start_where) >= 2:
                                    for a in range(90 - d_line_pos[0]):
                                        if D_germline_mismatch[0][d_line_pos[0] + a] != '.':
                                            d_cigar_mismatch.append([seq_before_d_num + d_line_pos[0] + a, D_germline_mismatch[0][d_line_pos[0] + a]])
                                    for a in range(int(d_germ_finish_where[1]) - int(d_germ_start_where[0]) + 1 - 90 + d_line_pos[0]):
                                        if D_germline_mismatch[1][a] != '.':
                                            d_cigar_mismatch.append([seq_before_d_num + 90 + a, D_germline_mismatch[1][a]])
                                d_cigar_mismatch.sort(key=lambda x: x[0])
                                for kk in range(len(d_cigar_mismatch)):
                                    isit = 0
                                    for ll in range(len(d_cigar_mismatches)):
                                        if d_cigar_mismatch[kk][0] == d_cigar_mismatches[ll][0]:
                                            isit += 1
                                    if isit == 0:
                                        d_cigar_mismatches.append(d_cigar_mismatch[kk])

                                j_finish_from_zero = 0
                                if len(j_germ_start_where) != 0:
                                    j_finish_from_zero = int(j_germ_finish_where[len(j_germ_finish_where) - 1]) - int(
                                        j_germ_start_where[0])

                                j_cigar_mismatch = []
                                j_line_pos = []
                                if (len(j_germ_start_where)) != 0:
                                    for k in range(90):
                                        if J_germline_mismatch[0][k] != '-':
                                            j_line_po = k
                                            j_line_pos.append(j_line_po)
                                            break

                                seq_before_j_num = 0
                                if j_finish_from_zero != 0 and len(v_germ_finish_where) != 0:
                                    if (len(d_germ_finish_where)) <= 1:
                                        for q in range(len(v_germ_start_where) - 1):
                                            # seq_before_j_num : 90, 180, ...
                                            seq_before_j_num = seq_before_j_num + 90
                                    elif (len(d_germ_finish_where)) == 2:
                                        for q in range(len(v_germ_start_where)):
                                            # seq_before_j_num : 90, 180, ...
                                            seq_before_j_num = seq_before_j_num + 90
                                    for a in range(len(j_cigar_mismatch)):
                                        if sequence_alignment[seq_before_j_num + j_line_pos[0] + 90 * a] == '-':
                                            j_cigar_mismatch.append([seq_before_j_num + j_line_pos[0] + 90 * a, '0'])
                                if len(j_germ_start_where) == 1:
                                    for i in range(int(j_germ_finish_where[0]) - int(j_germ_start_where[0])):
                                        if J_germline_mismatch[0][j_line_pos[0] + i] != '.':
                                            j_cigar_mismatch.append(
                                                [i + seq_before_j_num + j_line_pos[0],
                                                 J_germline_mismatch[0][j_line_pos[0] + i]])
                                elif len(j_germ_start_where) >= 2:
                                    for a in range(90 - j_line_pos[0]):
                                        if J_germline_mismatch[0][j_line_pos[0] + a] != '.':
                                            j_cigar_mismatch.append([seq_before_j_num + j_line_pos[0] + a,
                                                                     J_germline_mismatch[0][j_line_pos[0] + a]])
                                    for pp in range(len(j_germ_start_where) - 1):
                                        for a in range(int(j_germ_finish_where[pp + 1]) - int(j_germ_start_where[pp + 1]) + 1):
                                            if J_germline_mismatch[pp + 1][a] != '.':
                                                j_cigar_mismatch.append([seq_before_j_num + 90 * (pp + 1) + a, J_germline_mismatch[1][a]])
                                j_cigar_mismatch.sort(key=lambda x: x[0])

                                for kk in range(len(j_cigar_mismatch)):
                                    isit = 0
                                    for ll in range(len(j_cigar_mismatches)):
                                        if j_cigar_mismatch[kk][0] == j_cigar_mismatches[ll][0]:
                                            isit += 1
                                    if isit == 0:
                                        j_cigar_mismatches.append(j_cigar_mismatch[kk])

                                germline_alignment = sequence_alignment
                                germline_alignment_one = ''
                                germline_alignment_two = ''
                                germline_alignment_three = ''
                                if len(j_cigar_mismatches) != 0 and len(germline_alignment) != 0:
                                    for p in range(len(j_cigar_mismatches)):
                                        germline_alignment = germline_alignment[: j_cigar_mismatches[p][0]] + j_cigar_mismatches[p][1] + germline_alignment[
                                                                                                                                         j_cigar_mismatches[p][0] + 1:]
                                if len(d_cigar_mismatches) != 0 and len(germline_alignment) != 0:
                                    for p in range(len(d_cigar_mismatches)):
                                        germline_alignment = germline_alignment[: d_cigar_mismatches[p][0]] + d_cigar_mismatches[p][1] + germline_alignment[
                                                                                                                                         d_cigar_mismatches[p][0] + 1:]
                                if len(j_cigar_mismatch) != 0 and len(germline_alignment) != 0:
                                    for p in range(len(j_cigar_mismatch)):
                                        germline_alignment = germline_alignment[: j_cigar_mismatch[p][0]] + j_cigar_mismatch[p][1] + germline_alignment[j_cigar_mismatch[p][0] + 1:]
                                the_string = ''
                                second_string = ''
                                if (len(d_germ_start_where)) == 0 and (len(j_germ_finish_where)) != 0 and len(v_germ_finish_where) != 0:
                                    seq_before_j_num = 0
                                    if j_finish_from_zero != 0 and len(v_germ_finish_where) != 0:
                                        if (len(d_germ_finish_where)) <= 1:
                                            for q in range(len(v_germ_start_where) - 1):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + 90
                                        elif (len(d_germ_finish_where)) == 2:
                                            for q in range(len(v_germ_start_where)):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + 90
                                    the_string = ''
                                    for u in range(seq_before_j_num + j_line_pos[0] - int(v_germ_finish_where[(len(v_germ_finish_where) - 1)]) + int(v_germ_start_where[0])):
                                        the_string = the_string + '-'
                                    if len(germline_alignment) != 0:
                                        germline_alignment = germline_alignment[
                                                             :int(v_germ_finish_where[len(v_germ_finish_where) - 1]) - int(v_germ_start_where[0])] + the_string + germline_alignment[
                                                                                                                                                                  seq_before_j_num +
                                                                                                                                                                  j_line_pos[0]:]
                                elif (len(d_germ_start_where)) != 0 and len(j_germ_finish_where) != 0 and len(v_germ_finish_where) != 0:
                                    seq_before_d_num = 0
                                    for q in range(len(v_germ_start_where) - 1):
                                        # seq_before_d_num : 90, 180, ...
                                        seq_before_d_num = seq_before_d_num + int(v_germ_finish_where[0]) - int(v_germ_start_where[0]) + 1

                                    seq_before_j_num = 0
                                    if j_finish_from_zero != 0 and len(v_germ_finish_where) != 0:
                                        if (len(d_germ_finish_where)) <= 1:
                                            for q in range(len(v_germ_start_where) - 1):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + int(v_germ_finish_where[0]) - int(
                                                    v_germ_start_where[0]) + 1
                                        elif (len(d_germ_finish_where)) == 2:
                                            for q in range(len(v_germ_start_where)):
                                                # seq_before_j_num : 90, 180, ...
                                                seq_before_j_num = seq_before_j_num + int(v_germ_finish_where[0]) - int(
                                                    v_germ_start_where[0]) + 1

                                    the_string = ''
                                    second_string = ''
                                    for u in range(seq_before_d_num + d_line_pos[0] - int(v_germ_finish_where[(len(v_germ_finish_where) - 1)]) + int(v_germ_start_where[0]) - 1):
                                        the_string = the_string + '-'
                                    for u in range(seq_before_j_num + j_line_pos[0] - seq_before_d_num - d_line_pos[0] - int(d_germ_finish_where[len(d_germ_finish_where) - 1]) + int(
                                            d_germ_start_where[0]) - 1):
                                        second_string = second_string + '-'
                                    germline_alignment_one = germline_alignment[:int(v_germ_finish_where[len(v_germ_finish_where) - 1]) - int(v_germ_start_where[0]) + 1]
                                    germline_alignment_two = germline_alignment[seq_before_d_num + d_line_pos[0]: seq_before_d_num + d_line_pos[0] + int(
                                        d_germ_finish_where[len(d_germ_finish_where) - 1]) - int(d_germ_start_where[0]) + 1]
                                    germline_alignment_three = germline_alignment[seq_before_j_num + j_line_pos[0]:]
                                    germline_alignment = germline_alignment_one + the_string + germline_alignment_two + second_string + germline_alignment_three

                                result['the_string'] = the_string
                                result['second_string'] = second_string
                                result['germline_alignment'] = germline_alignment
                                result['germline_alignment_1'] = germline_alignment_one
                                result['germline_alignment_2'] = germline_alignment_two
                                result['germline_alignment_3'] = germline_alignment_three
                                result['j_line_pos'] = j_line_pos

                                v_cig = ''
                                cur_pos_is = -1
                                type_mis = 0
                                v_cig_num = 0
                                v_cigar = ''
                                if (len(v_cigar_mismatches) >= 1):
                                    cur_pos_is = -1
                                    for e in range(len(v_cigar_mismatches)):
                                        if v_cigar_mismatches[e][1] == '0':
                                            v_cig_num = int(v_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if v_cig_num == -1:
                                                v_cig = v_cig
                                            elif v_cig_num != 0:
                                                v_cig = v_cig + str(v_cig_num) + '=' + str(1) + 'D'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            elif type_mis == 1:
                                                v_cig = v_cig[:len(v_cig) - 2] + str(int(v_cig[len(v_cig) - 2]) + 1) + 'D'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            else:
                                                v_cig = v_cig + str(1) + 'D'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            type_mis = 1
                                        elif v_cigar_mismatches[e][1] == '-':
                                            v_cig_num = int(v_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if v_cig_num == -1:
                                                v_cig = v_cig
                                            elif v_cig_num != 0:
                                                v_cig = v_cig + str(v_cig_num) + '=' + str(1) + 'I'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            elif type_mis == 2:
                                                v_cig = v_cig[:len(v_cig) - 2] + str(int(v_cig[len(v_cig) - 2]) + 1) + 'I'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            else:
                                                v_cig = v_cig + str(1) + 'I'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            type_mis = 2
                                        else:
                                            v_cig_num = int(v_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if v_cig_num == -1:
                                                v_cig = v_cig
                                            elif v_cig_num != 0:
                                                v_cig = v_cig + str(v_cig_num) + '=' + str(1) + 'X'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            elif type_mis == 3:
                                                v_cig = v_cig[:len(v_cig) - 2] + str(int(v_cig[len(v_cig) - 2]) + 1) + 'X'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            else:
                                                v_cig = v_cig + str(1) + 'X'
                                                cur_pos_is = int(v_cigar_mismatches[e][0])
                                            type_mis = 3
                                if len(v_germ_finish_where) != 0:
                                    if cur_pos_is != int(v_germ_finish_where[(len(v_germ_finish_where) - 1)]) - int(v_germ_start_where[0]) + 1:
                                        v_cig_num = int(v_germ_finish_where[(len(v_germ_finish_where) - 1)]) - int(v_germ_start_where[0]) - cur_pos_is
                                        v_cig = v_cig + str(v_cig_num) + '='
                                v_cigar = v_cig

                                d_cig = ''
                                type_mis = 0
                                d_cig_num = 0
                                d_cigar = ''
                                if len(d_germ_finish_where) != 0 and len(v_germ_finish_where) != 0:
                                    cur_pos_is = seq_before_d_num + d_line_pos[0] - 1
                                    for e in range(len(d_cigar_mismatches)):
                                        if d_cigar_mismatches[e][1] == '0':
                                            d_cig_num = int(d_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if d_cig_num == -1:
                                                d_cig = d_cig
                                            elif d_cig_num != 0:
                                                d_cig = d_cig + str(d_cig_num) + '=' + str(1) + 'D'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            elif type_mis == 1:
                                                d_cig = d_cig[:len(d_cig) - 2] + str(int(d_cig[len(d_cig) - 2]) + 1) + 'D'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            else:
                                                d_cig = d_cig + str(1) + 'D'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            type_mis = 1
                                        elif d_cigar_mismatches[e][1] == '-':
                                            d_cig_num = int(d_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if d_cig_num == -1:
                                                d_cig = d_cig
                                            elif d_cig_num != 0:
                                                d_cig = d_cig + str(d_cig_num) + '=' + str(1) + 'I'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            elif type_mis == 2:
                                                d_cig = d_cig[:len(d_cig) - 2] + str(int(d_cig[len(d_cig) - 2]) + 1) + 'I'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            else:
                                                d_cig = d_cig + str(1) + 'I'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            type_mis = 2
                                        else:
                                            d_cig_num = int(d_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if d_cig_num == -1:
                                                d_cig = d_cig
                                            elif d_cig_num != 0:
                                                d_cig = d_cig + str(d_cig_num) + '=' + str(1) + 'X'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            elif type_mis == 3:
                                                d_cig = d_cig[:len(d_cig) - 2] + str(int(d_cig[len(d_cig) - 2]) + 1) + 'X'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            else:
                                                d_cig = d_cig + str(1) + 'X'
                                                cur_pos_is = int(d_cigar_mismatches[e][0])
                                            type_mis = 3
                                    if cur_pos_is != int(d_germ_finish_where[(len(d_germ_finish_where) - 1)]) - int(
                                            d_germ_start_where[0]) + seq_before_d_num + d_line_pos[0] + 1:
                                        d_cig_num = int(d_germ_finish_where[(len(d_germ_finish_where) - 1)]) - int(
                                            d_germ_start_where[0]) - cur_pos_is + seq_before_d_num + d_line_pos[0]
                                        d_cig = d_cig + str(d_cig_num) + '='
                                    d_cigar = d_cig
                                elif len(d_germ_finish_where) == 0:
                                    d_cigar = ''

                                j_cig = ''
                                type_mis = 0
                                j_cig_num = 0
                                j_cigar = ''
                                if (len(j_line_pos)) != 0:
                                    cur_pos_is = seq_before_j_num + j_line_pos[0] - 1
                                    for e in range(len(j_cigar_mismatches)):
                                        if j_cigar_mismatches[e][1] == '0':
                                            j_cig_num = int(j_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if j_cig_num == -1:
                                                j_cig = j_cig
                                            elif j_cig_num != 0:
                                                j_cig = j_cig + str(j_cig_num) + '=' + str(1) + 'D'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            elif type_mis == 1:
                                                j_cig = j_cig[:len(j_cig) - 2] + str(int(j_cig[len(j_cig) - 2]) + 1) + 'D'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            else:
                                                j_cig = j_cig + str(1) + 'D'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            type_mis = 1
                                        elif j_cigar_mismatches[e][1] == '-':
                                            j_cig_num = int(j_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if j_cig_num == -1:
                                                j_cig = j_cig
                                            elif j_cig_num != 0:
                                                j_cig = j_cig + str(j_cig_num) + '=' + str(1) + 'I'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            elif type_mis == 2:
                                                j_cig = j_cig[:len(j_cig) - 2] + str(int(j_cig[len(j_cig) - 2]) + 1) + 'I'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            else:
                                                j_cig = j_cig + str(1) + 'I'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            type_mis = 2
                                        else:
                                            j_cig_num = int(j_cigar_mismatches[e][0]) - cur_pos_is - 1
                                            if j_cig_num == -1:
                                                j_cig = j_cig
                                            elif j_cig_num != 0:
                                                j_cig = j_cig + str(j_cig_num) + '=' + str(1) + 'X'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            elif type_mis == 3:
                                                j_cig = j_cig[:len(j_cig) - 2] + str(int(j_cig[len(j_cig) - 2]) + 1) + 'X'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            else:
                                                j_cig = j_cig + str(1) + 'X'
                                                cur_pos_is = int(j_cigar_mismatches[e][0])
                                            type_mis = 3
                                    if cur_pos_is != int(j_germ_finish_where[(len(j_germ_finish_where) - 1)]) - int(
                                            j_germ_start_where[0]) + seq_before_j_num + j_line_pos[0] + 1:
                                        j_cig_num = int(j_germ_finish_where[(len(j_germ_finish_where) - 1)]) - int(
                                            j_germ_start_where[0]) - cur_pos_is + seq_before_j_num + j_line_pos[0]
                                        j_cig = j_cig + str(j_cig_num) + '='
                                    j_cigar = j_cig

                                result['v_cigar'] = v_cigar
                                result['d_cigar'] = d_cigar
                                result['j_cigar'] = j_cigar
                                result['V_mis_ratio'] = V_mis_ratio
                                result['D_mis_ratio'] = D_mis_ratio
                                result['J_mis_ratio'] = J_mis_ratio
                                result['V_germline_mismatch'] = V_germline_mismatch
                                result['D_germline_mismatch'] = D_germline_mismatch
                                result['J_germline_mismatch'] = J_germline_mismatch
                                result['v_germ_start_where'] = v_germ_start_where
                                result['v_germ_finish_where'] = v_germ_finish_where
                                result['d_germ_start_where'] = d_germ_start_where
                                result['d_germ_finish_where'] = d_germ_finish_where
                                result['j_germ_start_where'] = j_germ_start_where
                                result['j_germ_finish_where'] = j_germ_finish_where
                                result['v_cigar_mismatches'] = v_cigar_mismatches
                                result['d_cigar_mismatches'] = d_cigar_mismatches
                                result['j_cigar_mismatches'] = j_cigar_mismatches

                            elif line[:len('Alignment summary')] == 'Alignment summary':
                                totalDistanceFromVGene = 0
                                totalAlignedLength = 0

                                while True:
                                    next_line_sp = igb_handle.readline().split('\t')
                                    region = next_line_sp[0].replace('-IMGT', '').lower()
                                    if region == 'total':
                                        break
                                    if region == 'FR1':
                                        fr1_start = int(next_line_sp[1]) - 1
                                        fr1_end = int(next_line_sp[2])
                                        if sHelper.RepresentsInt(next_line_sp[6]):
                                            fr1_gap = int(next_line_sp[6])
                                        else:
                                            fr1_gap = 0
                                        while ((fr1_end - fr1_start + fr1_gap) % 3 != 0):
                                            fr1_start += 1
                                        result['fr1_from'] = fr1_start
                                        result['fr1_to'] = fr1_end
                                        result['fr1_length'] = int(next_line_sp[3])
                                        result['fr1_matches'] = int(next_line_sp[4])
                                        result['fr1_mismatches'] = int(next_line_sp[5])
                                        result['fr1_gaps'] = int(next_line_sp[6])
                                        totalAlignedLength += result['fr1_length']
                                        totalDistanceFromVGene += (result['fr1_mismatches'] + result['fr1_gaps'])
                                    else:
                                        try:
                                            result[region + '_from'] = int(next_line_sp[1]) - 1
                                            result[region + '_to'] = int(next_line_sp[2])
                                            result[region + '_length'] = int(next_line_sp[3])
                                            result[region + '_matches'] = int(next_line_sp[4])
                                            result[region + '_mismatches'] = int(next_line_sp[5])
                                            result[region + '_gaps'] = int(next_line_sp[6])
                                            totalAlignedLength += result[region + '_length']
                                            totalDistanceFromVGene += (result[region + '_mismatches'] + result[region + '_gaps'])
                                        except ValueError:
                                            continue

                                result['v_alignment_length'] = totalAlignedLength
                                result['v_alignment_mutation'] = totalDistanceFromVGene

                            elif line[:len('V-(D)-J rearrangement summary')] == 'V-(D)-J rearrangement summary':
                                next_line_sp = igb_handle.readline().split('\t')
                                if self.chain_type in [ChainType.HUMAN_HEAVY, ChainType.RABBIT_HEAVY,
                                                       ChainType.MOUSE_C57BL6_HEAVY, ChainType.MOUSE_BALBC_HEAVY, ChainType.MOUSE_HEAVY,
                                                       ChainType.HUMAN_BETA, ChainType.HUMAN_DELTA,
                                                       ChainType.MOUSE_BETA, ChainType.MOUSE_DELTA]:
                                    result['v_call'] = next_line_sp[0]
                                    result['d_call'] = next_line_sp[1]
                                    result['j_call'] = next_line_sp[2]
                                    result['vj_frame'] = next_line_sp[5].split('-')[0]
                                    result['productive'] = next_line_sp[6]

                                elif self.chain_type in [ChainType.HUMAN_LIGHT, ChainType.RABBIT_KAPPA,
                                                         ChainType.MOUSE_C57BL6_LIGHT, ChainType.MOUSE_BALBC_LIGHT, ChainType.MOUSE_LIGHT,
                                                         ChainType.HUMAN_ALPHA, ChainType.HUMAN_GAMMA,
                                                         ChainType.MOUSE_ALPHA, ChainType.MOUSE_GAMMA]:
                                    result['v_call'] = next_line_sp[0]
                                    result['j_call'] = next_line_sp[1]
                                    result['vj_frame'] = next_line_sp[4].split('-')[0]
                                    result['productive'] = next_line_sp[6]

                            # End of one alignment result
                            elif line[:len('Effective search space used:')] == 'Effective search space used:':
                                igblast_parsed_dict[row_num] = result
                                # igblast_parsed_list.append(result)
                                break


        ### combine annotation information with input_data_list (isotyped input file)
        # get range_from and range_to
        range_from = key_start
        range_to = max(igblast_parsed_dict.keys())

        # extract target sequence list
        target_list = input_data_list[range_from:(range_to+1)]

        # get output dictionary
        result_dict = {}
        w_row_num = key_start
        for target, igblast_parsed_data in zip(target_list, igblast_parsed_dict.values()):
            tmp_dict = {'sequence': target['full_NT'], 'duplicate_count': target['readcount'], 'c_call': target['isotype'], 'umi_list': target['umi_list']}
            igblast_parsed_data.update(tmp_dict)

            ## check not annotated cases and fill it with default info (N/A or blank)
            if 'sequence_aa' not in igblast_parsed_data:
                igblast_parsed_data['sequence_aa'] = ''
            if igblast_parsed_data['productive'] == 'Yes':
                igblast_parsed_data['productive'] = 'T'
            else:
                igblast_parsed_data['productive'] = 'F'


            # add cdr1, 2
            sequence = target['full_NT']
            if igblast_parsed_data['rev_comp'] == 'T':
                sequence = sHelper.to_reverse_complement(sequence)
            if 'fr4_from' in igblast_parsed_data:
                igblast_parsed_data['fr4'] = sequence[igblast_parsed_data['fr4_from']:igblast_parsed_data['fr4_from'] + 33]
            for r in ['fr1', 'fr2', 'fr3', 'cdr1', 'cdr2']:
                if ('%s_from' % r) in igblast_parsed_data and ('%s_to' % r) in igblast_parsed_data:
                    igblast_parsed_data[r] = sequence[igblast_parsed_data['%s_from' % r]:igblast_parsed_data['%s_to' % r]]
                    igblast_parsed_data[r + '_aa'] = sHelper.translate(igblast_parsed_data[r])
                else:
                    igblast_parsed_data[r] = 'N/A'
                    igblast_parsed_data[r + '_aa'] = 'X'

            result_dict[w_row_num] = igblast_parsed_data
            w_row_num += 1

        self.out_chunk = result_dict


    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue


class functionalExtracter(Process):
    """ Do extraction of functional reads

        """

    def __init__(self, input_file, log_queue, RAM_used_queue,
                 in_queue=None, out_queue=None, chunk_size=None, delim='\t'):
        Process.__init__(self)

        self.input_file = input_file

        self.out_chunk = {}

        self.log_queue = log_queue
        self.RAM_used_queue = RAM_used_queue

        self.in_queue = in_queue
        self.out_queue = out_queue

        self.chunk_size = chunk_size

        self.delim = delim

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

            # extract functional reads
            self.extract_functional(chunk_start_pos)
            self.in_queue.task_done()

            # RAM usage check
            tmp_RAM_used = psutil.virtual_memory().used / 1024 ** 3
            self.RAM_used_queue.put(tmp_RAM_used)

            # save chunk data into result queue
            self.out_queue.put(self.out_chunk)

        # logging.info("Finished run")


    # work function
    def extract_functional(self, key_start):
        # get output dictionary
        result_dict = {}
        w_row_num = key_start

        with open(self.input_file, 'r') as igb_handle:
            key_start_count = 0
            accessed_key_count = 0

            header = igb_handle.readline().strip().split(self.delim)        # extract header line
            col_seq_aa = header.index('sequence_aa')
            col_v_call = header.index('v_call')
            col_j_call = header.index('j_call')

            col_cdr1_aa = header.index('cdr1_aa')
            col_cdr2_aa = header.index('cdr2_aa')
            col_cdr3_aa = header.index('cdr3_aa')

            col_rev_comp = header.index('rev_comp')

            while True:
                line = igb_handle.readline()
                if not line:
                    break

                # go to start position
                if key_start_count < key_start:
                    key_start_count += 1
                    continue

                # stop working if the number of accessed key exceed the chunk_size
                accessed_key_count += 1
                if accessed_key_count > self.chunk_size:
                    break

                data_list = line.strip().split(self.delim)
                seq_aa = data_list[col_seq_aa]
                v_call = data_list[col_v_call]
                j_call = data_list[col_j_call]

                cdr1_aa = data_list[col_cdr1_aa]
                cdr2_aa = data_list[col_cdr2_aa]
                cdr3_aa = data_list[col_cdr3_aa]

                rev_comp = data_list[col_rev_comp]


                # full aa stop codon filtering
                try:
                    seq_aa.index('*')
                    continue
                except ValueError:
                    pass

                # full aa frame filtering
                try:
                    seq_aa.index('X')
                    continue
                except ValueError:
                    pass

                # gene annotation filtering
                if v_call == 'N/A' or j_call == 'N/A':
                    continue

                # cdr1,2,3 annotation flitering
                if cdr1_aa == 'N/A' or cdr2_aa == 'N/A' or cdr3_aa == 'N/A':
                    continue

                # cdr1,2,3 stop codon filtering
                try:
                    cdr1_aa.index('*')
                    continue
                except ValueError:
                    pass
                try:
                    cdr2_aa.index('*')
                    continue
                except ValueError:
                    pass
                try:
                    cdr3_aa.index('*')
                    continue
                except ValueError:
                    pass

                # cdr1,2,3 frame filtering
                try:
                    cdr1_aa.index('X')
                    continue
                except ValueError:
                    pass
                try:
                    cdr2_aa.index('X')
                    continue
                except ValueError:
                    pass
                try:
                    cdr3_aa.index('X')
                    continue
                except ValueError:
                    pass

                # rev_comp filtering: must be 'F'
                if rev_comp == 'T':
                    continue

                result_dict[w_row_num] = {k:v for k, v in zip(header, data_list)}
                w_row_num += 1

        self.out_chunk = result_dict


    def get_result_queue(self):
        return self.out_queue

    def get_RAM_used_queue(self):
        return self.RAM_used_queue