# bio-panning binder mapping to the BCR HC libraries

import csv
import os
import pandas as pd
import numpy as np
import re
from Common.SeqUtil import Helper
from Common.FileUtil import write_file

WORK_DIR = '/home/Backup1/SARS-CoV-2_vaccination/221015_biopanning_NGS_jw/'
REP_DIR = os.path.join(WORK_DIR, '221015_repertoire_sum')
PATIENT_REP_DIR = os.path.join(WORK_DIR, '221015_repertoire_sum')
MAPPING_DIR = os.path.join(WORK_DIR, '220226_CoV-AbDab_mapping')
TSNE_DIR = os.path.join(MAPPING_DIR, 'tSNE')
NAIVE_DIR = os.path.join(MAPPING_DIR, 'naive_effect_quantification')
BIOPAN_DIR = os.path.join(WORK_DIR, '221016_mapping')
OMICRON_DIR = os.path.join(WORK_DIR, '221016_mapping')

if True:
    patient_list = [str(i + 1) for i in range(55)]

    seq_list = ['IGHV1-69_IGHJ4_ARVRGYSGYGASGYFDN',
                'IGHV5-51_IGHJ4_ATTYHYDTDGPYGEFYY',
                'IGHV3-30_IGHJ3_ARTGSGWTDAFDI',
                'IGHV1-69_IGHJ3_ARVHGYSGYGANDAFDI',
                'IGHV1-69_IGHJ4_ARAEDHGTYYSDSSGYHFDY',
                'IGHV1-69_IGHJ4_AREPGILGYCSSTSCYID',
                'IGHV1-69_IGHJ5_AREVGYSGFGASPKFDP',
                'IGHV1-69_IGHJ5_AREIGYSGSGSAKYFDP',
                'IGHV3-53_IGHJ6_ARDLMEAGGMDV']
    seq_search = {k:None for k in seq_list}

    target_thresh = 0.0

    # Find perfect match
    seq_write_set = {}
    for seq in seq_list:
        seq_write_set[seq] = []

    # helper = Helper()

    for patient in patient_list:
        print(patient)
        patient_dir = os.path.join(PATIENT_REP_DIR, patient)
        target_filename = patient + '_new.tsv'

        with open(os.path.join(patient_dir, target_filename), 'r') as f:
            header = f.readline().strip().split('\t')
            col_dict = {}
            for h in header:
                col_dict[h] = header.index(h)
            while True:
                line = f.readline()
                if not line:
                    break
                data = line.strip().split('\t')
                vgene = data[col_dict['v_call']].split('*')[0]
                jgene = data[col_dict['j_call']].split('*')[0]
                cdr3aa = data[col_dict['cdr3_aa']]
                v3j = vgene + '_' + jgene + '_' + cdr3aa
                if v3j in seq_search:
                    line_dict = {h:data[col_dict[h]] for h in header}
                    seq_write_set[v3j].append(line_dict)
    
    for seq, write_list in seq_write_set.items():
        write_filename = seq + '_th' + str(target_thresh) + 'jw3_1418.tsv'
        if len(write_list) == 0:
            continue
        with open(os.path.join(BIOPAN_DIR, write_filename), 'w') as handle:
            write_file(handle, header, write_list, delim='\t')