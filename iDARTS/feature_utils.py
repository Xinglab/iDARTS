# -*- coding: UTF-8 -*-

"""
iDARTS - feature_utils
"""

import string
import os
import yaml
import pyBigWig
from pysam import FastaFile
from pkg_resources import resource_filename
import numpy as np
import regex as re
from collections import defaultdict
from skimage.util.shape import view_as_windows
from collections import Counter
import pandas as pd
from subprocess import *
from tqdm import tqdm
from keras.models import model_from_json

from . import config

REVERSE_COMPLIMENT_TAB = string.maketrans("ACTG", "TGAC")
DATA_CONFIG = yaml.safe_load(open(config.iDARTS_DIRECTORY + 'resource.yaml', 'r'))
DATA_DIR = DATA_CONFIG['data_dir']
HG19SEQ_FASTA_OBJ = FastaFile(os.path.join(DATA_DIR, 'hg19.fa'))
HG19_ALU_REPEAT = resource_filename('iDARTS.resources', 'hg19_alu.bed')
BURGE_CHASIN_ESE_ESS = resource_filename('iDARTS.resources.features', 'ESS_ESE.info.txt')
ISE_ISS_MOTIF = resource_filename('iDARTS.resources.features', 'ISS_ISE.info.txt')
RNA_BINDING_MOTIF = resource_filename('iDARTS.resources.features', 'RNA_Binding_Motif.txt')
RNA_BINDING_PSSM = resource_filename('iDARTS.resources.features', 'RBP_PSSMs')
ESE2Matrix_INFO = resource_filename('iDARTS.resources.features', 'ESE2Matrix.info')

def makedirs(_dir):
    try:
        os.stat(_dir)
    except:
        os.makedirs(_dir)

def remove_file(fn):
    cmd = 'rm -rf ' + fn
    call(cmd, shell = True)

def read_file(filename):
    with open(filename) as fp:
        List = [x.strip() for x in fp if len(x.strip()) > 0]
        return List

def get_file_num_lines(filename):
    fp = open(filename)
    count = 0
    for line in fp:
        count += 1
    fp.close()
    return count

def one_hot_encode(seq):
    map = np.asarray([[0, 0, 0, 0],
                      [1, 0, 0, 0],
                      [0, 1, 0, 0],
                      [0, 0, 1, 0],
                      [0, 0, 0, 1]])
    seq = seq.upper().replace('A', '\x01').replace('C', '\x02')
    seq = seq.replace('G', '\x03').replace('T', '\x04').replace('N', '\x00')
    return map[np.frombuffer(seq, np.int8) % 5]

def seq_fetch(chrom, start, end, strand, mutate_pos_list = [], mutate_base_list = []):
    start, end = int(start), int(end)
    seq = HG19SEQ_FASTA_OBJ.fetch(chrom, start, end).upper()
    if len(mutate_pos_list) > 0:
        seq = list(seq)
        for mutate_pos, mutate_base in zip(mutate_pos_list, mutate_base_list):
            mutate_pos = int(mutate_pos)
            if mutate_pos >= start and mutate_pos < end:
                seq[mutate_pos - start] = mutate_base
        seq = ''.join(seq)
    if strand == '-':
        seq = seq.translate(REVERSE_COMPLIMENT_TAB)[::-1]
    return seq

def load_splice_feature_list(splice_feature_fn):
    splice_feature_list = []
    fp = open(splice_feature_fn)
    fp.readline()
    for line in fp:
        sp = line.strip().split('\t')
        splice_feature_list.append(sp[1])
    return splice_feature_list

def load_alu_counts(alu_code_fn):
    dict_alu_counts = {}
    fp = open(alu_code_fn)
    for line in fp:
        if len(line.strip()) > 0:
            sp = line.strip().split('\t')
            key = '_'.join(sp[0:4])
            dict_alu_counts[key] = int(sp[-1])
    fp.close()
    return dict_alu_counts

def load_rna_binding_motif_pssm():
    dict_rna_binding_pssm = {}
    path_list = os.listdir(RNA_BINDING_PSSM)
    for path in path_list:
        if 'human' in path:
            pssm = []
            for line in read_file(RNA_BINDING_PSSM + '/' + path)[1:]:
                sp = line.strip().split('\t')[1:]
                sp = [float(x) for x in sp]
                pssm.append(sp)
            sq = path.strip().split('_')
            dict_rna_binding_pssm[sq[0] + '.' + sq[1] + '.pssm'] = np.array(pssm)
    return dict_rna_binding_pssm

def get_pssm_score(seq_matrix, dict_rna_binding_pssm):
    seq_len = len(seq_matrix)
    idx = {}
    rank_weight = {}
    if seq_len >= 6:
        idx[6] = view_as_windows(seq_matrix, (6, 4))
        idx_len = len(idx[6])
        rank_weight[6] = np.logspace(-idx_len + 1, 0, idx_len, base = np.e)[::-1]
    if seq_len >= 7:
        idx[7] = view_as_windows(seq_matrix, (7, 4))
        idx_len = len(idx[7])
        rank_weight[7] = np.logspace(-idx_len + 1, 0, idx_len, base = np.e)[::-1]
    if seq_len >= 8:
        idx[8] = view_as_windows(seq_matrix, (8, 4))
        idx_len = len(idx[8])
        rank_weight[8] = np.logspace(-idx_len + 1, 0, idx_len, base = np.e)[::-1]
    score_list = []
    for rbp in sorted(dict_rna_binding_pssm.keys()):
        pssm = np.array(dict_rna_binding_pssm[rbp])
        pssm_len = len(pssm)
        if len(pssm) > seq_len:
            score = 0
        else:
            vectorized_score = pssm * idx[pssm_len]
            scores = vectorized_score.sum(axis = (2,3)).reshape(-1,)
            scores = np.sort(scores)[::-1]
            score = np.sum(scores * rank_weight[pssm_len])
        score_list.append(score)
    return score_list

def translatable(seq):  # translatable without stop codon all three frames tested
    seq1_3, seq2_3, seq3_3 = re.findall('...', seq), re.findall('...', seq[1:]), re.findall('...', seq[2:])
    if ('TAG' in seq1_3 or 'TAA' in seq1_3 or 'TGA' in seq1_3) and (
                        'TAG' in seq2_3 or 'TAA' in seq2_3 or 'TGA' in seq2_3) and (
                        'TAG' in seq3_3 or 'TAA' in seq3_3 or 'TGA' in seq3_3):
        return 0
    else:
        return 1

def load_ESE_matrix_PSSM():
    matrix_list = read_file(ESE2Matrix_INFO)
    dict_ESE_matrix = {}
    dict_ESE_matrix['SF2_ASF'] = np.array([x.split()[1:] for x in matrix_list[2:6]], dtype='float32').T
    dict_ESE_matrix['SC35'] = np.array([x.split()[1:] for x in matrix_list[8:12]], dtype='float32').T
    dict_ESE_matrix['SRp40'] = np.array([x.split()[1:] for x in matrix_list[14:18]], dtype='float32').T
    dict_ESE_matrix['SRp55'] = np.array([x.split()[1:] for x in matrix_list[20:24]], dtype='float32').T
    return dict_ESE_matrix

def get_ESE_matrix_score(seq_matrix, dict_ESE_matrix): # http://krainer01.cshl.edu/tools/ESE2/ESEmatrix.html
    seq_len = len(seq_matrix)
    idx = {}
    rank_weight = {}
    if seq_len >= 6:
        idx[6] = view_as_windows(seq_matrix, (6, 4))
        idx_len = len(idx[6])
        rank_weight[6] = np.logspace(-idx_len + 1, 0, idx_len, base = np.e)[::-1]
    if seq_len >= 7:
        idx[7] = view_as_windows(seq_matrix, (7, 4))
        idx_len = len(idx[7])
        rank_weight[7] = np.logspace(-idx_len + 1, 0, idx_len, base = np.e)[::-1]
    if seq_len >= 8:
        idx[8] = view_as_windows(seq_matrix, (8, 4))
        idx_len = len(idx[8])
        rank_weight[8] = np.logspace(-idx_len + 1, 0, idx_len, base = np.e)[::-1]
    score_list = []
    for rbp in sorted(dict_ESE_matrix.keys()):
        pssm = np.array(dict_ESE_matrix[rbp])
        pssm_len = len(pssm)
        if len(pssm) > seq_len:
            score = 0
        else:
            vectorized_score = pssm * idx[pssm_len]
            scores = vectorized_score.sum(axis = (2,3)).reshape(-1,)
            scores = np.sort(scores)[::-1]
            score = np.sum(scores * rank_weight[pssm_len])
            score = np.mean(scores)
        score_list.append(score)
    return score_list

def load_phatsConsScore():
    chromIndex = ['chr' + str(x) for x in range(1, 23)] + ['chrX', 'chrY']
    dict_chrom_phatsConsScore = {}
    for chrom in chromIndex:
        dict_chrom_phatsConsScore[chrom] = pyBigWig.open(os.path.join(DATA_DIR, '{0}.phastCons46way.bw'.format(chrom)))
    return dict_chrom_phatsConsScore

def load_rna_binding_motif():
    dict_rna_binding_motif = {}
    for line in read_file(RNA_BINDING_MOTIF):
        sp = line.strip().split()
        dict_rna_binding_motif[sp[0]] = sp[-1]
    return dict_rna_binding_motif

def load_esr_matrix():
    file_name = resource_filename('iDARTS.resources.features', 'chasin_ESS_ESE_numbers.txt')
    esr_temp = pd.read_csv(file_name, header=None, names=["wtMotif", "ESR_chasin_score", "ESR_type"], sep = '\t')
    esr_matrix = dict()
    for i in range(0, len(esr_temp)):
        esr_matrix[esr_temp["wtMotif"][i]] = esr_temp["ESR_chasin_score"][i]
    return esr_matrix

def load_rosenberg_matrix():
    A3SS_R1 = resource_filename('iDARTS.resources.features', 'mean_effects_sdpos_A3SS_R1.txt')
    A3SS_R2 = resource_filename('iDARTS.resources.features', 'mean_effects_sdpos_A3SS_R2.txt')
    A5SS_R1 = resource_filename('iDARTS.resources.features', 'mean_effects_sdpos_A5SS_R1.txt')
    A5SS_R2 = resource_filename('iDARTS.resources.features', 'mean_effects_sdpos_A5SS_R2.txt')
    A3SS_R1 = pd.read_csv(A3SS_R1, sep = '\t')
    A3SS_R2 = pd.read_csv(A3SS_R2, sep = '\t')
    A5SS_R1 = pd.read_csv(A5SS_R1, sep = '\t')
    A5SS_R2 = pd.read_csv(A5SS_R2, sep = '\t')
    rosenberg_matrix = dict()
    for i, k in [(A3SS_R1, "A3SS_R1"), (A3SS_R2, "A3SS_R2"), (A5SS_R1, "A5SS_R1"), (A5SS_R2, "A5SS_R2")]:
        rosenberg_matrix[k] = dict()
        for j in range(0, len(i)):
            rosenberg_matrix[k][i["motif"][j]] = i[k + "_mean_effects"][j]
    return rosenberg_matrix

def load_cnn_splice_predictor():
    acceptor_json = resource_filename('iDARTS.resources.model', 'acceptor.json')
    acceptor_weight = resource_filename('iDARTS.resources.model', 'acceptor.h5')
    donor_json = resource_filename('iDARTS.resources.model', 'donor.json')
    donor_weight = resource_filename('iDARTS.resources.model', 'donor.h5')
    acceptor_cnn = model_from_json(open(acceptor_json).read())
    acceptor_cnn.load_weights(acceptor_weight)
    donor_cnn = model_from_json(open(donor_json).read())
    donor_cnn.load_weights(donor_weight)
    return acceptor_cnn, donor_cnn

def rnafold_predict_pu_score(seq, tmp_file, event_type):  # bed format, not include the end_position
    tmp_dir = os.path.dirname(os.path.abspath(tmp_file)) + '/'
    file_name = os.path.basename(tmp_file)
    cwd_path = os.getcwd()
    os.chdir(tmp_dir)
    fw = open(file_name, 'w')
    fw.write('>{}\n'.format('iDARTS_seq_' + event_type))
    fw.write(seq + '\n')
    fw.close()
    seq_len = len(seq)
    cmd = 'RNAplfold -u {} -W {} < {}'.format(1, min(70, seq_len), file_name)
    Popen(cmd, shell = True, stdout=PIPE, stderr=PIPE).communicate()
    lunp_probs = open('iDARTS_seq_' + event_type + '_lunp')
    os.chdir(cwd_path)
    lunp_probs.readline()
    lunp_probs.readline()
    unpair_prob_list = []
    for line in lunp_probs:
        unpair_prob = line.strip().split('\t')[1]
        if 'nan' in unpair_prob:
            unpair_prob = 0
        unpair_prob_list.append(float(unpair_prob))
    return unpair_prob_list
    
def get_rosenberg_mean_score(seq, rosenberg_matrix, rosenberg_region):
    dict_rosenberg_score = defaultdict(lambda:0)
    dict_rosenberg_count = defaultdict(lambda:0)
    for i in xrange(0, len(seq) - 6 + 1):
        sub_seq = seq[i:i+6]
        for region in rosenberg_region:
            if sub_seq in rosenberg_matrix[region]:
                dict_rosenberg_score[region] += rosenberg_matrix[region][sub_seq]
                dict_rosenberg_count[region] += 1
    dict_rosenberg_mean_score = {}
    for region in rosenberg_region:
        if dict_rosenberg_count[region] == 0:
            dict_rosenberg_mean_score[region] = 0.0
        else:
            dict_rosenberg_mean_score[region] = dict_rosenberg_score[region] / dict_rosenberg_count[region]
    return dict_rosenberg_mean_score

def get_esr_mean_score(seq, esr_matrix):
    esr_score, esr_count = 0.0, 0.0
    for i in xrange(0, len(seq) - 6 + 1):
        sub_seq = seq[i:i+6]
        if sub_seq in esr_matrix:
            esr_score += esr_matrix[sub_seq]
            esr_count += 1
    return esr_score / esr_count
