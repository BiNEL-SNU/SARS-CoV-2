# 1) CoV-AbDab mapping to the BCR HC libraries

import os
from nbformat import write
import pandas as pd
import numpy as np
import re
import glob
from shutil import copyfile
from xlsxwriter.workbook import Workbook
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import matplotlib as mpl
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, FixedFormatter
import sys
sys.path.append('/home/jwonchoi20/ip_code/')
from Common.SeqUtil import Helper
from Common.ProcessUtil import run_clustal_omega
from Common.FileUtil import read_file, write_file, read_file_fast, read_fasta
from Common.Enum import ChainType
from Common.EdgeExtractor import extract_edge_all_hamming_under_thres

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.1f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def extract_nn_distance_from_edge_files(vertex_file: str, edge_file: str, out_file: str, seq_id_colname: str, region_type: str):
    vertexHeader, vertexList = read_file_fast(vertex_file, delim='\t')
    edgeHeader, edgeList = read_file_fast(edge_file, delim='\t')

    vertexEdgeSet = {}

    for vertex in vertexList:
        vertexEdgeSet[vertex[seq_id_colname]] = {'sequence_length': len(vertex[region_type]),
                                                 'edges': []}

    for edge in edgeList:
        try:
            vertexEdgeSet[edge['from']]['edges'].append(edge['dist'])
            vertexEdgeSet[edge['to']]['edges'].append(edge['dist'])
        except KeyError:
            print('Vertex name error occured.')
            raise RuntimeError

    nnWriteList = []
    for vertexID in vertexEdgeSet:
        tempSet = vertexEdgeSet[vertexID]
        if len(tempSet['edges']) == 0:
            continue

        nnWriteList.append({'sequence_id': vertexID,
                            'nn_distance': min(tempSet['edges']),
                            'normalized_nn_distance': float(min(tempSet['edges'])) / tempSet['sequence_length'],
                            'sequence_length': tempSet['sequence_length']})

    with open(out_file, 'w') as handle:
        write_file(handle=handle,
                    header=['sequence_id', 'nn_distance', 'normalized_nn_distance', 'sequence_length'],
                    data=nnWriteList, delim='\t')

def extract_nn_distance_from_edge_files_btw_diff_stages(vertex_file: str, edge_file: str, seq_id_colname: str, region_type: str,
                                                        write_dir: str, sample_type: str, edge_type: str, ref_stage_list: list):
    vertexHeader, vertexList = read_file_fast(vertex_file, delim='\t')
    edgeHeader, edgeList = read_file_fast(edge_file, delim='\t')

    vertexEdgeSet = {}
    nnWriteSet = {}

    for ref_stage in ref_stage_list:
        vertexEdgeSet['_'.join(ref_stage)] = {}
        nnWriteSet['_'.join(ref_stage)] = []

    for ref_stage in ref_stage_list:
        for vertex in vertexList:
            vertexEdgeSet['_'.join(ref_stage)][vertex[seq_id_colname]] = {'sequence_length': len(vertex[region_type]),
                                                                          'edges': []}

    for edge in edgeList:
        from_stage = edge['from'].split('-')[0]
        to_stage = edge['to'].split('-')[0]

        stage_list = list(set([from_stage, to_stage]))
        if len(stage_list) == 1:
            continue

        try:
            stage_list.index('pre')
            stage_list.sort(reverse=True)
        except ValueError:
            stage_list.sort()


        try:
            vertexEdgeSet['_'.join(stage_list)]
        except KeyError:
            continue

        try:
            vertexEdgeSet['_'.join(stage_list)][edge['from']]['edges'].append(edge['dist'])
            vertexEdgeSet['_'.join(stage_list)][edge['to']]['edges'].append(edge['dist'])
        except KeyError:
            print('Vertex name error occured.')
            raise RuntimeError


    for stage, stageEdgeSet in vertexEdgeSet.items():
        for vertexID in stageEdgeSet:
            tempSet = stageEdgeSet[vertexID]
            if len(tempSet['edges']) == 0:
                continue

            nnWriteSet[stage].append({'sequence_id': vertexID,
                                      'nn_distance': min(tempSet['edges']),
                                      'normalized_nn_distance': float(min(tempSet['edges'])) / tempSet['sequence_length'],
                                      'sequence_length': tempSet['sequence_length']})

    for stage, nnWriteList in nnWriteSet.items():
        out_filename = '_'.join([sample_type, stage, edge_type, 'nnDist']) + '.tsv'
        out_file = os.path.join(write_dir, out_filename)
        with open(out_file, 'w') as handle:
            write_file(handle=handle,
                        header= ['sequence_id', 'nn_distance', 'normalized_nn_distance', 'sequence_length'],
                        data= nnWriteList, delim='\t')
            
exist_list = [i + 1 for i in range(55)]

pre_rem_list = [2, 3, 46]
first_rem_list = [51]
second_rem_list = []
third_rem_list = []
fourth_rem_list = [2, 12, 20, 42, 46]
fifth_rem_list = [2, 8, 12, 18, 19, 20, 25, 26,
                  31, 37, 42, 46]

pre_exist_list = exist_list.copy()
first_exist_list = exist_list.copy()
second_exist_list = exist_list.copy()
third_exist_list = exist_list.copy()
fourth_exist_list = exist_list.copy()
fifth_exist_list = exist_list.copy()

for pre_rem in pre_rem_list:
    pre_exist_list.remove(pre_rem)
for first_rem in first_rem_list:
    first_exist_list.remove(first_rem)
for second_rem in second_rem_list:
    second_exist_list.remove(second_rem)
for third_rem in third_rem_list:
    third_exist_list.remove(third_rem)
for fourth_rem in fourth_rem_list:
    fourth_exist_list.remove(fourth_rem)
for fifth_rem in fifth_rem_list:
    fifth_exist_list.remove(fifth_rem)

exist_set = {'pre':pre_exist_list,
             '1st':first_exist_list,
             '2nd':second_exist_list,
             '3rd':third_exist_list,
             '4th':fourth_exist_list,
             '5th':fifth_exist_list}

colors = sns.color_palette('colorblind')

WORK_DIR = '/home/Backup1/SARS-CoV-2_vaccination'
REP_DIR = os.path.join(WORK_DIR, 'NovaSeq_data')
PATIENT_REP_OLD_DIR = os.path.join(WORK_DIR, '220225_NGS_data_patient')
PATIENT_REP_DIR = os.path.join(WORK_DIR, '221015_biopanning_NGS_jw', '221015_repertoire_sum')
MAPPING_DIR = os.path.join(WORK_DIR, '221016_CoV-AbDab_mapping')
if not os.path.exists(MAPPING_DIR):
    os.mkdir(MAPPING_DIR)
TSNE_DIR = os.path.join(MAPPING_DIR, 'tSNE')
NAIVE_DIR = os.path.join(MAPPING_DIR, 'naive_effect_quantification')

### Mapping in clonotype level
if False:
    helper = Helper()

    target_thresh = 0.0

    cov_ver = '031022'
    cov_dir = os.path.join(WORK_DIR, 'CoV_AbDab', cov_ver)
    cov_case = ''
    cov_filename = 'CoV-AbDab_' + cov_ver + '_refTyped_v2%s.tsv' % cov_case

    cov_header, cov_data = read_file_fast(os.path.join(cov_dir, cov_filename), delim='\t')

    ## set_num = 0 or 1 or 2
    # 0 : whole repertoires
    # 1 : For subjects with complete set of pre~5th repertoires with > 50,000 reads
    # 2 : For subjects with complete set of pre~5th repertoires with > 10,000 reads
    case_table = {0: 'whole', 1: 'set_1', 2: 'set_2', 3: 'set_2nd+5th'}
    set_num = 2
    if set_num == 0:
        patient_list = [str(i) for i in range(1, 56)]
    if set_num == 1:
        patient_list = [4, 5, 6, 9, 10, 13, 14, 15, 16, 17,
                        22, 23, 24, 27, 28, 29, 30, 32, 33, 35,
                        36, 38, 39, 40, 41, 43, 44, 45, 48, 49,
                        50, 52, 53, 54]
        patient_list = [str(i) for i in patient_list]
    elif set_num == 2:
        patient_list = [1, 4, 5, 6, 7, 9, 10, 11, 13, 14,
                        15, 16, 17, 21, 22, 23, 24, 27, 28, 29,
                        30, 32, 33, 34, 35, 36, 38, 39, 40, 41,
                        43, 44, 45, 47, 48, 49, 50, 52, 53, 54, 55]
        patient_list = [str(i) for i in patient_list]        
    elif set_num == 3:
        patient_list = [3,4,5,6,7,9,10,13,14,15,16,17,21,22,23,24,27,28,29,30,32,33,35,36,38,39,40,41,43,44,45,48,49,50,51,52,53,54,55]
        patient_list = [str(i) for i in patient_list]
    temp_write_dir = os.path.join(MAPPING_DIR, case_table[set_num])
    write_dir = os.path.join(temp_write_dir, 'th' + str(target_thresh) + 'not_neu_to')
    for t_dir in [temp_write_dir, write_dir]:
        os.makedirs(t_dir, exist_ok=True)
    
    print('Mapping to CoV' + str(cov_ver) + ' th' + str(target_thresh) + ' on set %d' % set_num)

    cov_set = {}
    if target_thresh == 0.0:
        for d in cov_data:
            each_key = '|'.join([d['Heavy V Gene'], d['Heavy J Gene'], d['CDRH3']])
            try:
                cov_set[each_key][0].append(d['ref_name'])
                cov_set[each_key][1].append(d['bind_ref_type'])
                cov_set[each_key][2].append(d['neu_ref_type'])
                cov_set[each_key][3].append(d['Not Neutralising Vs'])
            except KeyError:
                cov_set[each_key] = [[],[],[]]
                cov_set[each_key][0].append(d['ref_name'])
                cov_set[each_key][1].append(d['bind_ref_type'])
                cov_set[each_key][2].append(d['neu_ref_type'])
                cov_set[each_key][3].append(d['Not Neutralising Vs'])

        for _, temp_list in cov_set.items():
            cov_set[_][0] = '&'.join(temp_list[0])
            cov_set[_][1] = '&'.join(temp_list[1])
            cov_set[_][2] = '&'.join(temp_list[2])
            cov_set[_][3] = '&'.join(temp_list[3])

        for patient in patient_list:
            print('CoV-AbDab_mapping to ' + patient)
            patient_dir = os.path.join(PATIENT_REP_DIR, patient)
            patient_filename = patient + '_new.tsv'

            write_list = []
            write_filename = 's' + patient + '_to_CoV' + cov_ver + '.tsv'
            
            with open(os.path.join(patient_dir, patient_filename), 'r') as handle:
                header = handle.readline().strip().split('\t')
                col_dict = {h:header.index(h) for h in header}
                while True:
                    line = handle.readline()
                    if not line:
                        break
                    d = line.strip().split('\t')
                    target_cdr3 = d[col_dict['cdr3_aa']]
                    target_v = d[col_dict['v_call']].split('*')[0]
                    target_j = d[col_dict['j_call']].split('*')[0]
                    target_key = '|'.join([target_v, target_j, target_cdr3])
                    if target_key in cov_set:
                        data_dict = {h:d[col_dict[h]] for h in header}
                        data_dict['subject'] = patient
                        data_dict['mapped_to'] = cov_set[target_key][0]
                        data_dict['mapped_bind'] = cov_set[target_key][1]
                        data_dict['mapped_neu'] = cov_set[target_key][2]
                        data_dict['not_neu_to'] = cov_set[target_key][3]
                        write_list.append(data_dict)

            w_header = ['subject'] + header + ['mapped_to', 'mapped_bind', 'mapped_neu', 'not_neu_to']
            with open(os.path.join(write_dir, write_filename), 'w') as handle:
                write_file(handle, w_header, write_list, delim='\t')

    else:
        for d in cov_data:
            each_vjlen = '|'.join([d['Heavy V Gene'], d['Heavy J Gene'], str(len(d['CDRH3']))])
            if each_vjlen not in cov_set:
                cov_set[each_vjlen] = {d['CDRH3']: [[d['ref_name']],[d['bind_ref_type']],[d['neu_ref_type']],[d['Not Neutralising Vs']]]}
            else:
                if d['CDRH3'] not in cov_set[each_vjlen]:
                    cov_set[each_vjlen][d['CDRH3']] = [[d['ref_name']],[d['bind_ref_type']],[d['neu_ref_type']],[d['Not Neutralising Vs']]]
                else:
                    cov_set[each_vjlen][d['CDRH3']][0].append(d['ref_name'])
                    cov_set[each_vjlen][d['CDRH3']][1].append(d['bind_ref_type'])
                    cov_set[each_vjlen][d['CDRH3']][2].append(d['neu_ref_type'])
                    cov_set[each_vjlen][d['CDRH3']][3].append(d['Not Neutralising Vs'])

        for vjlen, cdrh3_dict in cov_set.items():
            for cdrh3, temp_list in cdrh3_dict.items():
                cov_set[vjlen][cdrh3][0] = '&'.join(temp_list[0])
                cov_set[vjlen][cdrh3][1] = '&'.join(temp_list[1])
                cov_set[vjlen][cdrh3][2] = '&'.join(temp_list[2])
                cov_set[vjlen][cdrh3][3] = '&'.join(temp_list[3])
        for patient in patient_list:
            print('CoV-AbDab_mapping to ' + patient)
            patient_dir = os.path.join(PATIENT_REP_DIR, patient)
            patient_filename = patient + '_new.tsv'

            write_list = []
            write_filename = 's' + patient + '_to_CoV' + cov_ver + cov_case + '.tsv'
            
            with open(os.path.join(patient_dir, patient_filename), 'r') as handle:
                header = handle.readline().strip().split('\t')
                col_dict = {h:header.index(h) for h in header}
                while True:
                    line = handle.readline()
                    if not line:
                        break
                    d = line.strip().split('\t')
                    target_cdr3 = d[col_dict['cdr3_aa']]
                    target_v = d[col_dict['v_call']].split('*')[0]
                    target_j = d[col_dict['j_call']].split('*')[0]
                    target_key = '|'.join([target_v, target_j, str(len(target_cdr3))])
                    if target_key in cov_set:
                        cov_cdrh3_dict = cov_set[target_key]
                        for cov_cdrh3 in cov_cdrh3_dict:
                            dist = helper.hamming_distance(cov_cdrh3, target_cdr3) / len(target_cdr3)
                            if dist <= target_thresh:
                                data_dict = {h:d[col_dict[h]] for h in header}
                                data_dict['subject'] = patient
                                data_dict['CoV_cdr3_aa'] = cov_cdrh3
                                data_dict['mapped_to'] = cov_cdrh3_dict[cov_cdrh3][0]
                                data_dict['mapped_bind'] = cov_cdrh3_dict[cov_cdrh3][1]
                                data_dict['mapped_neu'] = cov_cdrh3_dict[cov_cdrh3][2]
                                data_dict['not_neu_to'] = cov_cdrh3_dict[cov_cdrh3][3]                                
                                write_list.append(data_dict)

            w_header = ['subject'] + header + ['CoV_cdr3_aa', 'mapped_to', 'mapped_bind', 'mapped_neu', 'not_neu_to']
            with open(os.path.join(write_dir, write_filename), 'w') as handle:
                write_file(handle, w_header, write_list, delim='\t')

# merge the separated mapped sequences
if False:
    set_list = ['set_2']
    target_thresh = 0.0
    cov_ver = '031022'
    cov_case = ''
    write_filename = 'Lineage_mapping_to_CoV' + cov_ver + '_th' + str(target_thresh) + '_by_clonotype_v2_not_neu_to.tsv'
    # write_filename = 'V3J_mapping_to_CoV' + cov_ver + '_th' + str(target_thresh) + cov_case + '.tsv'
    
    for set_case in set_list:
        case_dir = os.path.join(MAPPING_DIR, set_case)
        target_dir = os.path.join(case_dir, 'th' + str(target_thresh) + 'not_neu_to')
        write_list = []

        for idx, temp_file in enumerate(glob.glob(os.path.join(target_dir, '*CoV%s%s.tsv' % (cov_ver, cov_case)))):
            temp_df = pd.read_csv(temp_file, sep='\t')
            write_list.append(temp_df)
        write_df = pd.concat(write_list)

        print('\n')
        print(str(idx+1) + ' subjcts were merged.')

        write_df.to_csv(os.path.join(case_dir, write_filename),
                        header=True,
                        index=False,
                        sep='\t')