#####################################################################################################
# Author    HJ Han, JS Noh
# Date      2021-08-03 ~
# Editor    YHLee
# Revised
# Note      Extract distance files for NN plot
#
#####################################################################################################

import os
import logging
import glob
from datetime import datetime
from multiprocessing import Manager, Process
from itertools import combinations
import numpy as np

# add user module finding path
import sys
sys.path.append('/home/yonghee91/IP_python')

from Common.FileUtil import read_file_new, read_file_fast, write_file
from Common.SeqUtil import Helper


### old functions: calculate levenshtein distance
class EdgeExtracterLD(Process):
    def __init__(self, output_dir: str, distance_threshold: int, col_id: str, col_seq: str, map_length_data: str,
                 in_queue):
        Process.__init__(self)
        self.in_queue = in_queue

        self.output_dir = output_dir
        self.dist_cutoff = distance_threshold
        self.col_id = col_id
        self.col_seq = col_seq
        self.map_length_data = map_length_data

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            item = self.in_queue.get()
            if item is None:
                self.in_queue.task_done()
                break

            self.do_work(item)
            self.in_queue.task_done()

    def do_work(self, length: int):
        cnt_compare, cnt_edge = 0, 0
        temp_start = datetime.now()

        data_l = self.map_length_data[length]
        line_edges = []
        # brute-force
        for i, d_from in enumerate(data_l):
            for d_to in data_l[i + 1:]:
                dist = Helper.levenshtein_distance(d_from[self.col_seq], d_to[self.col_seq])
                cnt_compare += 1
                if dist <= self.dist_cutoff:
                    if d_from[self.col_id] < d_to[self.col_id]:
                        line_edges.append('%s\t%s\t%d\n' % (d_from[self.col_id], d_to[self.col_id], dist))
                    else:
                        line_edges.append('%s\t%s\t%d\n' % (d_to[self.col_id], d_from[self.col_id], dist))
                    cnt_edge += 1

        for length_to in range(length + 3, length + self.dist_cutoff + 1, 3):
            if length_to in self.map_length_data:
                for d_to in self.map_length_data[length_to]:
                    dist = Helper.levenshtein_distance(d_from[self.col_seq], d_to[self.col_seq])
                    cnt_compare += 1
                    if dist <= self.dist_cutoff:
                        if d_from[self.col_id] < d_to[self.col_id]:
                            line_edges.append('%s\t%s\t%d\n' % (d_from[self.col_id], d_to[self.col_id], dist))
                        else:
                            line_edges.append('%s\t%s\t%d\n' % (d_to[self.col_id], d_from[self.col_id], dist))
                        cnt_edge += 1

        print("%s do_work ---------- takes %s, length %d, %d seqs, %d comps, %d edges (Brute-force)" % (
            datetime.now(), datetime.now() - temp_start, length, len(data_l), cnt_compare, cnt_edge))

        if len(line_edges) > 0:
            with open(os.path.join(self.output_dir, 'edges.length-%d.tsv' % length), 'w') as f:
                f.writelines(line_edges)

def extract_edge_all_levenstein_under_thres(file_in_nodes: str, file_out_edges: str, col_id: str, col_seq: str,
                                            threshold: int = 30, num_process: int = 4):
    logging.info('Started extract_edge_levenstein_all_under_thres: %s' % os.path.basename(file_in_nodes))
    time_start = datetime.now()

    temp_start = datetime.now()
    header, data = read_file_new(file_in_nodes, delim='\t')
    logging.info("---------- takes %s : read_file(file_in_nodes)" % (datetime.now() - temp_start))

    temp_start = datetime.now()
    map_length_data = {}
    for d in data:
        if len(d[col_seq]) not in map_length_data:
            map_length_data[len(d[col_seq])] = []
        map_length_data[len(d[col_seq])].append(d)
    logging.info("---------- takes %s : make map_length_data" % (datetime.now() - temp_start))

    dir_edges = os.path.join(os.path.dirname(file_out_edges), os.path.basename(file_out_edges).replace('.tsv', ''))
    if not os.path.exists(dir_edges):
        os.makedirs(dir_edges)

    # Make another worker thread to send query
    manager = Manager()
    in_queue = manager.Queue()
    worker_list = []
    logging.info("Spawning %d EdgeExtracterLD processes.", num_process)
    for i in range(num_process):
        worker = EdgeExtracterLD(dir_edges, threshold, col_id, col_seq, map_length_data, in_queue)
        worker.daemon = True
        worker.start()

        worker_list.append(worker)

    jobs = list(map_length_data.items())
    jobs.sort(key=lambda kv: len(kv[1]), reverse=True)

    for length, data_l in jobs:
        in_queue.put_nowait(length)

    # enqueue stopper
    for i in range(len(worker_list)):
        in_queue.put_nowait(None)

    for worker in worker_list:
        worker.join()

    logging.info('Start merging')

    edges_files = glob.glob(os.path.join(dir_edges, '*.tsv'))
    with open(file_out_edges, 'w') as fw:
        fw.write('from\tto\tdist\n')
        for fr in edges_files:
            with open(fr, 'r') as hr:
                for line in hr:
                    fw.write(line)

    logging.info('Finished extract_edge_levenstein_all_under_thres: %s' % os.path.basename(file_in_nodes))
    logging.info('------------- takes total %s' % (datetime.now() - time_start))


### new functions: calculated hamming distance on each VJ-len groups
class EdgeExtractorHD(Process):
    def __init__(self, distance_threshold: int, col_id: str, col_seq: str, grouped_data: dict, in_queue, out_queue):
        Process.__init__(self)

        self.in_queue = in_queue
        self.out_queue = out_queue
        self.dist_cutoff = distance_threshold
        self.col_id = col_id
        self.col_seq = col_seq
        self.grouped_data = grouped_data
        self.out_chunk = []

    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            group_id = self.in_queue.get()
            if group_id is None:
                self.in_queue.task_done()
                break

            self.extract_dist(group_id)
            self.in_queue.task_done()

            self.out_queue.put(self.out_chunk)

    def extract_dist(self, group):
        # cnt_compare, cnt_edge = 0, 0

        target_group_data = self.grouped_data[group]
        line_edges = []

        # version 1. in-house code, using two combined loops
        if False:
            for i in range(0, len(target_group_data)-1):
                d_from = target_group_data[i]
                id_from = d_from[self.col_id]
                seq_from = d_from[self.col_seq]

                for j in range((i + 1), len(target_group_data)):
                    d_to = target_group_data[j]
                    id_to = d_to[self.col_id]
                    seq_to = d_to[self.col_seq]
                    dist = Helper.hamming_distance(seq_from, seq_to)
                    if dist > self.dist_cutoff:
                        continue

                    if id_from < id_to:
                        tmp_dict = {'from': id_from, 'to': id_to, 'dist': dist}
                    else:
                        tmp_dict = {'from': id_to, 'to': id_from, 'dist': dist}
                    line_edges.append(tmp_dict)

        # version 2. using combinations function of itertools module
        if True:
            data_list = [(d[self.col_id], d[self.col_seq]) for d in target_group_data]
            for compare in combinations(data_list, 2):
                seq_from = compare[0][1]
                seq_to = compare[1][1]
                id_from = compare[0][0]
                id_to = compare[1][0]
                dist = Helper.hamming_distance(seq_from, seq_to)
                if dist > self.dist_cutoff:
                    continue

                if id_from < id_to:
                    tmp_dict = {'from': id_from, 'to': id_to, 'dist': dist, 'region_len': max(len(seq_from), len(seq_to))}
                else:
                    tmp_dict = {'from': id_to, 'to': id_from, 'dist': dist, 'region_len': max(len(seq_from), len(seq_to))}
                line_edges.append(tmp_dict)

        self.out_chunk = line_edges

    def get_result_queue(self):
        return self.out_queue


def extract_edge_all_hamming_under_thres(file_in_nodes: str, file_out_edges: str, col_id: str, col_seq: str,
                                         threshold: int = 10, num_process: int = 5, vj_grouping: bool = True, delim='\t', **kwargs):
    logging.info('Script started on %s' % os.path.basename(file_in_nodes))
    time_start = datetime.now()
    header, data = read_file_fast(file_in_nodes, delim=delim)
    logging.info("--- %s seconds elapsed for reading file" % (datetime.now() - time_start))

    temp_start = datetime.now()
    grouped_data = {}

    # for the case of calculating distances within VJ-len group
    if vj_grouping:
        # get v/j gene column names
        if 'col_v_gene' in kwargs:
            col_v_gene = kwargs['col_v_gene']
        else:
            col_v_gene = 'v_call'
        if 'col_j_gene' in kwargs:
            col_j_gene = kwargs['col_j_gene']
        else:
            col_j_gene = 'j_call'

        # test if the v/j gene column names are correct
        if col_v_gene not in header and col_j_gene not in header:
            print('Wrong column name of v gene & j gene...')
            print('  SOLUTION ==> extract_edge_all_hamming_under_thres(..., col_v_gene="column name of V gene", col_j_gene="column name of J gene")')
            raise KeyError
        elif col_v_gene not in header:
            print('Wrong column name of v gene...')
            print('  SOLUTION ==> extract_edge_all_hamming_under_thres(..., col_v_gene="column name of V gene")')
            raise KeyError
        elif col_j_gene not in header:
            print('Wrong column name of j gene...')
            print('  SOLUTION ==> extract_edge_all_hamming_under_thres(..., col_j_gene="column name of J gene")')
            raise KeyError

        # extract v, j gene and do grouping
        for d in data:
            v_gene = d[col_v_gene].split('*')[0]
            j_gene = d[col_j_gene].split('*')[0]
            seq = d[col_seq]
            group = '|'.join([v_gene, j_gene, str(len(seq))])
            tmp_dict = {col_id: d[col_id], col_seq: seq}
            if group not in grouped_data:
                grouped_data[group] = [tmp_dict]
            else:
                grouped_data[group].append(tmp_dict)

    else:
        # do grouping using sequence length
        for d in data:
            seq = d[col_seq]
            group = str(len(seq))
            tmp_dict = {col_id: d[col_id], col_seq: seq}
            if group not in grouped_data:
                grouped_data[group] = [tmp_dict]
            else:
                grouped_data[group].append(tmp_dict)

    # for group, d_list in grouped_data.items():
    #     print(group)
    #     print(d_list)

    sorted_grouped_data = {k:v for k, v in sorted(grouped_data.items(), key=lambda x: len(x[1]), reverse=True)}

    logging.info("--- %s seconds elapsed for making grouped_data (type=dict)" % (datetime.now() - temp_start))


    # Make another worker thread to send query and get results
    temp_start = datetime.now()
    manager = Manager()
    in_queue = manager.Queue()
    result_queue = manager.Queue()
    for group in sorted_grouped_data:
        in_queue.put_nowait(group)
    worker_list = []

    logging.info("Spawning %d EdgeExtractorHD processes.", num_process)
    for i in range(num_process):
        worker = EdgeExtractorHD(distance_threshold=threshold, col_id=col_id, col_seq=col_seq, grouped_data=grouped_data, in_queue=in_queue, out_queue=result_queue)
        worker.daemon = True
        worker.start()
        worker_list.append(worker)

    # in_queue stopper
    for i in range(len(worker_list)):
        in_queue.put_nowait(None)

    # wait until all worker processes are done
    for worker in worker_list:
        worker.join()

    logging.info("--- %s seconds elapsed for calculating distance values" % (datetime.now() - temp_start))

    # load all results and make edge file
    temp_start = datetime.now()
    out_result_queue = worker.get_result_queue()
    edge_data = []
    while out_result_queue.qsize() > 0:
        result_chunk = out_result_queue.get()
        edge_data += result_chunk

    w_header = ['from', 'to', 'dist', 'region_len']
    with open(file_out_edges, 'w') as handle:
        write_file(handle, w_header, edge_data, delim=delim)

    logging.info('--- %s seconds elapsed for writing edge files' % (datetime.now() - temp_start))

    logging.info('Script finished on %s' % os.path.basename(file_in_nodes))
    logging.info('%s seconds elapsed for whole process' % (datetime.now() - time_start))



### functions for sharing analysis
### new functions: calculated hamming distance on each VJ-len groups
class EdgeExtractor_btw_lin_v1(Process):
    def __init__(self, out_dir: str, dist_thres: float, region: str, in_queue, delim: str='\t'):
        Process.__init__(self)
        self.out_dir = out_dir
        self.in_queue = in_queue
        self.dist_cutoff = dist_thres
        self.region = region
        self.delim = delim
        if self.delim == '\t':
            self.filetype = 'tsv'
        elif self.delim == ',':
            self.filetype = 'csv'


    def run(self):
        while True:
            # Get the work from the queue and expand the tuple
            compare = self.in_queue.get()
            if compare is None:
                self.in_queue.task_done()
                break

            self.get_dist_btw_reps(compare)
            self.in_queue.task_done()


    def get_dist_btw_reps(self, compare):
        ref_rep_file = compare[0]
        target_rep_file = compare[1]

        refSample = os.path.splitext(os.path.basename(ref_rep_file))[0]
        targetSample = os.path.splitext(os.path.basename(target_rep_file))[0]
        writeFilename = refSample + '_vs_' + targetSample + '.%s' % self.filetype

        refHeaer, refList = read_file_fast(ref_rep_file, delim=self.delim)
        targetHeader, targetList = read_file_fast(target_rep_file, delim=self.delim)

        refLinSet = {}
        targetLinSet = {}

        for ref in refList:
            try:
                refLinSet[ref['lineage_id']].append(ref)
            except KeyError:
                refLinSet[ref['lineage_id']] = [ref]

        for target in targetList:
            try:
                targetLinSet[target['lineage_id']].append(target)
            except KeyError:
                targetLinSet[target['lineage_id']] = [target]

        refSet = {}
        targetSet = {}

        for refLin in refLinSet:
            tempList = refLinSet[refLin]
            temp = tempList[0]
            refKey = '|'.join([temp['v_call'].split('*')[0], temp['j_call'].split('*')[0], str(len(temp[self.region]))])
            try:
                refSet[refKey][refLin] = refLinSet[refLin]
            except KeyError:
                refSet[refKey] = {}
                refSet[refKey][refLin] = refLinSet[refLin]

        for targetLin in targetLinSet:
            tempList = targetLinSet[targetLin]
            temp = tempList[0]
            targetKey = '|'.join([temp['v_call'].split('*')[0], temp['j_call'].split('*')[0], str(len(temp[self.region]))])

            try:
                targetSet[targetKey][targetLin] = targetLinSet[targetLin]
            except KeyError:
                targetSet[targetKey] = {}
                targetSet[targetKey][targetLin] = targetLinSet[targetLin]

        writeList = []

        for refKey in refSet:
            try:
                tempTargetLinSet = targetSet[refKey]
            except KeyError:
                continue
            tempRefLinSet = refSet[refKey]
            for tempRefLin in tempRefLinSet:
                tempRefLinList = tempRefLinSet[tempRefLin]
                for tempTargetLin in tempTargetLinSet:
                    tempTargetLinList = tempTargetLinSet[tempTargetLin]
                    linDist = self.cal_dist_btw_two_lin(ref_lin=tempRefLinList, target_lin=tempTargetLinList, region=self.region)
                    normLinDist = float(linDist)/len(tempRefLinList[0][self.region])
                    if normLinDist < self.dist_cutoff:
                        writeList.append({'from':tempRefLin, 'to':tempTargetLin, 'dist':normLinDist})

        with open(os.path.join(self.out_dir, writeFilename), 'w') as handle:
            write_file(handle, ['from', 'to', 'dist'], writeList, delim=self.delim)



    def cal_dist_btw_two_lin(self, ref_lin, target_lin, region):
        distList = []
        for ref in ref_lin:
            for target in target_lin:
                distList.append(Helper.hamming_distance(ref[region], target[region]))
        return (np.mean(distList))


def cal_lin_dist_btw_two_reps(rep_list: list, out_dir: str, region: str, dist_thresh: float, num_threads: int, delim: str='\t'):
    rep_combs = list(combinations(rep_list, 2))
    batch_size = num_threads
    batch_num = 1

    logging.info('The number of pairwise comparison: %d' % len(rep_combs))

    while True:
        if len(rep_combs) == 0:
            break

        batch_targets = rep_combs[:batch_size]
        rep_combs = rep_combs[batch_size:]
        logging.info('Batch number %d started...' % batch_num)

        manager = Manager()
        in_queue = manager.Queue()
        for job in batch_targets:
            in_queue.put(job)
        worker_list = []

        logging.info("Spawning %d EdgeExtractor_btw_lin_v1 processes.", num_threads)
        for i in range(num_threads):
            p = EdgeExtractor_btw_lin_v1(out_dir=out_dir, dist_thres=dist_thresh, region=region, in_queue=in_queue, delim=delim)
            p.daemon = True
            p.start()
            worker_list.append(p)

        # in_queue stopper
        for i in range(len(worker_list)):
            in_queue.put_nowait(None)

        # wait until all extracting works done
        for p in worker_list:
            p.join()

        logging.info('Batch number %d finished...' % batch_num)
        batch_num += 1



