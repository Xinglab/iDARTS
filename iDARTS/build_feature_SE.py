# -*- coding: UTF-8 -*-

"""
iDARTS - build_feature_SE
Implements a splicing feature build module for 
constructing exon skipping features
"""

from subprocess import *
from collections import defaultdict
from pkg_resources import resource_filename
import time
import logging
import tempfile
logger = logging.getLogger('iDARTS.build_feature')

from .feature_utils import *

CHROM_INDEX = 3  # start from 0 in input file
STRAND_INDEX = 4
UPSTREAM_ES_INDEX = 7
UPSTREAM_EE_INDEX = 8
EXON_START_INDEX = 5
EXON_END_INDEX = 6
DOWNSTREAM_ES_INDEX = 9
DOWNSTREAM_EE_INDEX = 10
POS0BASE_INDEX = 12
REF_INDEX = 13
ALT_INDEX = 14
TEMP_FILE_NAME = next(tempfile._get_candidate_names())

def get_region_coordinate(rmats_input_file_line):
    sp = rmats_input_file_line.strip().split('\t')
    _1 = [int(sp[UPSTREAM_ES_INDEX]), int(sp[UPSTREAM_EE_INDEX])]
    _2 = [int(sp[UPSTREAM_EE_INDEX]), int(sp[UPSTREAM_EE_INDEX]) + 300]
    _3 = [int(sp[EXON_START_INDEX]) - 300, int(sp[EXON_START_INDEX])]
    _4 = [int(sp[EXON_START_INDEX]), int(sp[EXON_END_INDEX])]
    _5 = [int(sp[EXON_END_INDEX]), int(sp[EXON_END_INDEX]) + 300]
    _6 = [int(sp[DOWNSTREAM_ES_INDEX]) - 300, int(sp[DOWNSTREAM_ES_INDEX])]
    _7 = [int(sp[DOWNSTREAM_ES_INDEX]), int(sp[DOWNSTREAM_EE_INDEX])]
    if sp[STRAND_INDEX] == '+':
        C1, A, C2 = _1, _4, _7
        I1_5p, I1_3p, I2_5p, I2_3p = _2, _3, _5, _6
    if sp[STRAND_INDEX] == '-':
        C1, A, C2 = _7, _4, _1
        I1_5p, I1_3p, I2_5p, I2_3p = _6, _5, _3, _2
    return C1, I1_5p, I1_3p, A, I2_5p, I2_3p, C2

def fetch_alu_code(args): # the input file (rMATS format)
    input_file = args.input
    event_type = args.event_type
    iDARTStmp_dir = os.path.dirname(os.path.abspath(input_file)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    makedirs(iDARTStmp_dir)
    tmpAluBed = iDARTStmp_dir + 'tmp.{}.alu.bed'.format(event_type)
    alufw = open(tmpAluBed, 'w')
    file_num_lines = get_file_num_lines(input_file)
    fp = open(input_file)
    if args.mutate == 'True':
        file_header = fp.readline().strip().split('\t')
        POS0BASE_INDEX = int(file_header.index('Pos0base'))
        REF_INDEX = int(file_header.index('Ref'))
        ALT_INDEX = int(file_header.index('Alt'))
    else:
        fp.readline()
    for line in fp:
        sp = line.strip().split('\t')
        chrom = sp[CHROM_INDEX]
        strand = sp[STRAND_INDEX]
        C1, I1_5p, I1_3p, A, I2_5p, I2_3p, C2 = get_region_coordinate(line)
        if sp[STRAND_INDEX] == '-':
            I1_3p_long = [sp[EXON_END_INDEX], int(sp[EXON_END_INDEX]) + 2000]
            I2_5p_long = [int(sp[EXON_START_INDEX]) - 2000, sp[EXON_START_INDEX]]
        else:
            I1_3p_long = [int(sp[EXON_START_INDEX]) - 2000, sp[EXON_START_INDEX]]
            I2_5p_long = [sp[EXON_END_INDEX], int(sp[EXON_END_INDEX]) + 2000]
        # Alu repeat writer...
        strInput = '{0}\t{1}\t{2}\tforward\t1\t+\n'
        alufw.write(strInput.format(chrom, I1_3p[0], I1_3p[1]))
        alufw.write(strInput.format(chrom, I2_5p[0], I2_5p[1]))
        alufw.write(strInput.format(chrom, I1_3p_long[0], I1_3p_long[1]))
        alufw.write(strInput.format(chrom, I2_5p_long[0], I2_5p_long[1]))
        strInput = '{0}\t{1}\t{2}\treverse\t1\t-\n'
        alufw.write(strInput.format(chrom, I1_3p[0], I1_3p[1]))
        alufw.write(strInput.format(chrom, I2_5p[0], I2_5p[1]))
        alufw.write(strInput.format(chrom, I1_3p_long[0], I1_3p_long[1]))
        alufw.write(strInput.format(chrom, I2_5p_long[0], I2_5p_long[1]))
    fp.close()
    alufw.close()
    cmd = 'sort {0} | uniq > {1}'.format(tmpAluBed, tmpAluBed + '.tmp')
    call(cmd, shell=True)
    cmd = 'cp {1} {0}'.format(tmpAluBed, tmpAluBed + '.tmp')
    call(cmd, shell=True)
    cmd = 'sort -k1,1 -k2,2n {0} > {1}'.format(tmpAluBed, tmpAluBed + '.tmp')
    call(cmd, shell=True)
    cmd = 'cp {1} {0}'.format(tmpAluBed, tmpAluBed + '.tmp')
    call(cmd, shell=True)
    cmd = 'bedtools intersect -a {0} -b {1} -s -sorted -wa > {2}'.format(tmpAluBed, HG19_ALU_REPEAT,
                                                                         tmpAluBed + '.result')
    call(cmd, shell=True)
    fw = open(iDARTStmp_dir + 'event.{}.alu'.format(event_type), 'w')
    fp = open(tmpAluBed + '.result')
    dict_alu = defaultdict(lambda: 0)
    for line in fp:
        if len(line.strip()) == 0:
            continue
        sp = line.strip().split('\t')
        key = '_'.join([sp[0], sp[1], sp[2], sp[-1]])
        dict_alu[key] += 1
    for key in dict_alu:
        sp = key.split('_')
        fw.write(sp[0] + '\t' + sp[1] + '\t' + sp[2] + '\t' + sp[-1] + '\t' + str(dict_alu[key]) + '\n')
    fp.close()
    fw.close()

def get_alu_feature_score(annotation_list, dict_alu_counts):
    dict_alu_feature_score = {}
    if annotation_list[STRAND_INDEX] == '+':
        I1_3p_plus = [int(annotation_list[EXON_START_INDEX]) - 300, annotation_list[EXON_START_INDEX]]
        I2_5p_plus = [annotation_list[EXON_END_INDEX], int(annotation_list[EXON_END_INDEX]) + 300]
        I1_3p_long_plus = [int(annotation_list[EXON_START_INDEX]) - 2000,
                           annotation_list[EXON_START_INDEX]]
        I2_5p_long_plus = [annotation_list[EXON_END_INDEX], int(annotation_list[EXON_END_INDEX]) + 2000]
    else:
        I1_3p_plus = [annotation_list[EXON_END_INDEX], int(annotation_list[EXON_END_INDEX]) + 300]
        I2_5p_plus = [int(annotation_list[EXON_START_INDEX]) - 300, annotation_list[EXON_START_INDEX]]
        I1_3p_long_plus = [annotation_list[EXON_END_INDEX], int(annotation_list[EXON_END_INDEX]) + 2000]
        I2_5p_long_plus = [int(annotation_list[EXON_START_INDEX]) - 2000,
                           annotation_list[EXON_START_INDEX]]
    chrom = annotation_list[CHROM_INDEX]
    dict_alu_region = {}
    dict_alu_region['I1_3p'] = I1_3p_plus
    dict_alu_region['I2_5p'] = I2_5p_plus
    dict_alu_region['I1_3p_long'] = I1_3p_long_plus
    dict_alu_region['I2_5p_long'] = I2_5p_long_plus
    for strand in ['Plus', 'Minus', 'Total']:
        for region in dict_alu_region:
            region_interval = dict_alu_region[region]
            if strand == 'Plus':
                strand_key = '{0}_{1}_{2}_{3}'.format(chrom, region_interval[0], region_interval[1], '+')
            if strand == 'Minus':
                strand_key = '{0}_{1}_{2}_{3}'.format(chrom, region_interval[0], region_interval[1], '-')
            if strand == 'Total':
                plus_strand_key = '{0}_{1}_{2}_{3}'.format(chrom, region_interval[0], region_interval[1], '+')
                minus_strand_key = '{0}_{1}_{2}_{3}'.format(chrom, region_interval[0], region_interval[1], '-')
                plus_strand_key_score = 0
                minus_strand_key_score = 0
                if plus_strand_key in dict_alu_counts:
                    plus_strand_key_score = dict_alu_counts[plus_strand_key]
                if minus_strand_key in dict_alu_counts:
                    minus_strand_key_score = dict_alu_counts[minus_strand_key]
                dict_alu_feature_score['Alu-' + region + '-' + strand] = plus_strand_key_score + minus_strand_key_score
            else:
                if strand_key in dict_alu_counts:
                    dict_alu_feature_score['Alu-' + region + '-' + strand] = dict_alu_counts[strand_key]
                else:
                    dict_alu_feature_score['Alu-' + region + '-' + strand] = 0
    return dict_alu_feature_score

def get_conservation_score(chrom, start, end, dict_chrom_phatsConsScore):
    cons_score_values = dict_chrom_phatsConsScore[chrom].values(chrom, start, end)
    cons_score_values = [x if x == x else 0 for x in cons_score_values]
    return cons_score_values

def get_regional_conservation_score(annotation_list, dict_chrom_phatsConsScore):
    strand = annotation_list[STRAND_INDEX]
    chrom = annotation_list[CHROM_INDEX]
    if strand == '+':
        I1_5p_index = int(annotation_list[UPSTREAM_EE_INDEX])
        I1_3p_index = int(annotation_list[EXON_START_INDEX])
        I2_5p_index = int(annotation_list[EXON_END_INDEX])
        I2_3p_index = int(annotation_list[DOWNSTREAM_ES_INDEX])
        I1_5p_conservation = get_conservation_score(chrom, I1_5p_index, I1_5p_index + 300, dict_chrom_phatsConsScore)
        I1_3p_conservation = get_conservation_score(chrom, I1_3p_index - 300, I1_3p_index, dict_chrom_phatsConsScore)
        I2_5p_conservation = get_conservation_score(chrom, I2_5p_index, I2_5p_index + 300, dict_chrom_phatsConsScore)
        I2_3p_conservation = get_conservation_score(chrom, I2_3p_index - 300, I2_3p_index, dict_chrom_phatsConsScore)
    if strand == '-':
        I2_3p_index = int(annotation_list[UPSTREAM_EE_INDEX])
        I2_5p_index = int(annotation_list[EXON_START_INDEX])
        I1_3p_index = int(annotation_list[EXON_END_INDEX])
        I1_5p_index = int(annotation_list[DOWNSTREAM_ES_INDEX])
        I1_5p_conservation = get_conservation_score(chrom, I1_5p_index - 300, I1_5p_index, dict_chrom_phatsConsScore)[::-1]
        I1_3p_conservation = get_conservation_score(chrom, I1_3p_index, I1_3p_index + 300, dict_chrom_phatsConsScore)[::-1]
        I2_5p_conservation = get_conservation_score(chrom, I2_5p_index - 300, I2_5p_index, dict_chrom_phatsConsScore)[::-1]
        I2_3p_conservation = get_conservation_score(chrom, I2_3p_index, I2_3p_index + 300, dict_chrom_phatsConsScore)[::-1]
    return I1_5p_conservation, I1_3p_conservation, I2_5p_conservation, I2_3p_conservation

def get_nova_binding_score(annotation_list, dict_seq_len, dict_seq, dict_region2conservation, dict_chrom_phatsConsScore):
    chrom = annotation_list[CHROM_INDEX]
    strand = annotation_list[STRAND_INDEX]
    dict_nova_region = {}
    NISE1_Seq = dict_seq['I1_3p'][-300:-35]
    NISE1_conservation = dict_region2conservation['I1_3p'][-300:-35]
    NISE2_Seq = dict_seq['I2_5p'][20:105]
    NISE2_conservation = dict_region2conservation['I2_5p'][20:105]
    NISE3_Seq = dict_seq['I2_3p'][-180:-75]
    NISE3_conservation = dict_region2conservation['I2_3p'][-180:-75]
    NISS1_Seq = dict_seq['I1_5p'][40:125]
    NISS1_conservation = dict_region2conservation['I1_5p'][40:125]
    UPSTREAM_ES_POS = int(annotation_list[UPSTREAM_ES_INDEX])
    UPSTREAM_EE_POS = int(annotation_list[UPSTREAM_EE_INDEX])
    EXON_START_POS = int(annotation_list[EXON_START_INDEX])
    EXON_END_POS = int(annotation_list[EXON_END_INDEX])
    DOWNSTREAM_ES_POS = int(annotation_list[DOWNSTREAM_ES_INDEX])
    DOWNSTREAM_EE_POS = int(annotation_list[DOWNSTREAM_EE_INDEX])
    if dict_seq_len['C1'] < 60:
        NESE1_Seq = dict_seq['C1'] + dict_seq['I1_5p'][0:5]
        if strand == '+':
            NESE1_conservation = get_conservation_score(chrom, UPSTREAM_ES_POS, UPSTREAM_EE_POS + 5, dict_chrom_phatsConsScore)
        if strand == '-':
            NESE1_conservation = get_conservation_score(chrom, DOWNSTREAM_ES_POS - 5, DOWNSTREAM_EE_POS, dict_chrom_phatsConsScore)[::-1]
    else:
        NESE1_Seq = dict_seq['C1'][-60:] + dict_seq['I1_5p'][0:5]
        if strand == '+':
            NESE1_conservation = get_conservation_score(chrom, UPSTREAM_EE_POS - 60, UPSTREAM_EE_POS + 5, dict_chrom_phatsConsScore)
        if strand == '-':
            NESE1_conservation = get_conservation_score(chrom, DOWNSTREAM_ES_POS - 5, DOWNSTREAM_ES_POS + 60, dict_chrom_phatsConsScore)[::-1]
    if dict_seq_len['A'] < 25:
        NISS2_Seq = dict_seq['I1_3p'][-60:] + dict_seq['A']
        if strand == '+':
            NISS2_conservation = get_conservation_score(chrom, EXON_START_POS - 60, EXON_END_POS, dict_chrom_phatsConsScore)
        if strand == '-':
            NESE1_conservation = get_conservation_score(chrom, EXON_START_POS, EXON_END_POS + 60, dict_chrom_phatsConsScore)[::-1]
            NISS2_conservation = [x if x == x else 0 for x in dict_chrom_phatsConsScore[chrom].values(chrom, int(annotation_list[EXON_START_INDEX]), int(annotation_list[EXON_END_INDEX]) + 60)][::-1]
    else:
        NISS2_Seq = dict_seq['I1_3p'][-60:] + dict_seq['A'][0:25]
        if strand == '+':
            NISS2_conservation = get_conservation_score(chrom, EXON_START_POS - 60, EXON_START_POS + 25, dict_chrom_phatsConsScore)
        if strand == '-':
            NISS2_conservation = get_conservation_score(chrom, EXON_END_POS - 25, EXON_END_POS + 60, dict_chrom_phatsConsScore)[::-1]
    if dict_seq_len['A'] < 65:
        NESS1_Seq = dict_seq['A'] + dict_seq['I2_5p'][0:65 - dict_seq_len['A']]
    else:
        NESS1_Seq = dict_seq['A'][0:65]
    if strand == '+':
        NESS1_conservation = get_conservation_score(chrom, EXON_START_POS, EXON_START_POS + 65, dict_chrom_phatsConsScore)
        NESS1_conservation = [x if x == x else 0 for x in dict_chrom_phatsConsScore[chrom].values(chrom, int(annotation_list[EXON_START_INDEX]), int(annotation_list[EXON_START_INDEX]) + 65)]
    if strand == '-':
        NESS1_conservation = get_conservation_score(chrom, EXON_END_POS - 65, EXON_END_POS, dict_chrom_phatsConsScore)[::-1]
    dict_nova_region['NESE1'] = {}
    dict_nova_region['NESE1']['seq'] = NESE1_Seq
    dict_nova_region['NESE1']['conservation'] = NESE1_conservation

    dict_nova_region['NISE1'] = {}
    dict_nova_region['NISE1']['seq'] = NISE1_Seq
    dict_nova_region['NISE1']['conservation'] = NISE1_conservation

    dict_nova_region['NISE2'] = {}
    dict_nova_region['NISE2']['seq'] = NISE2_Seq
    dict_nova_region['NISE2']['conservation'] = NISE2_conservation

    dict_nova_region['NISE3'] = {}
    dict_nova_region['NISE3']['seq'] = NISE3_Seq
    dict_nova_region['NISE3']['conservation'] = NISE3_conservation

    dict_nova_region['NISS1'] = {}
    dict_nova_region['NISS1']['seq'] = NISS1_Seq
    dict_nova_region['NISS1']['conservation'] = NISS1_conservation

    dict_nova_region['NISS2'] = {}
    dict_nova_region['NISS2']['seq'] = NISS2_Seq
    dict_nova_region['NISS2']['conservation'] = NISS2_conservation

    dict_nova_region['NESS1'] = {}
    dict_nova_region['NESS1']['seq'] = NESS1_Seq
    dict_nova_region['NESS1']['conservation'] = NESS1_conservation
    return dict_nova_region

def get_rnafold_predict_pu_score(dict_seq, tmp_file, event_type):
    dict_rnafold_pu_scores = {}
    if len(dict_seq['A']) >= 50:
        seq = dict_seq['A'][-50:] + dict_seq['I2_5p'][0:250]
    else:
        seq = dict_seq['I1_3p'][-(50 - len(dict_seq['A'])):] + dict_seq['A'] + dict_seq['I2_5p'][0:250]
    unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.1', event_type)
    for i in xrange(3):
        start, end = 70 * i, 70 * (i + 1)
        dict_rnafold_pu_scores['SecStr.max{}_{}.I2_5p'.format(start + 1, end)] = np.max(unpair_prob_list[start + 50:end + 50])
        dict_rnafold_pu_scores['SecStr.avg{}_{}.I2_5p'.format(start + 1, end)] = np.mean(unpair_prob_list[start + 50:end + 50])
    dict_rnafold_pu_scores['SecStr.maxJunc.I2_5p'] = np.max(unpair_prob_list[50:60])
    dict_rnafold_pu_scores['SecStr.avgJunc.I2_5p'] = np.mean(unpair_prob_list[50:60])

    if len(dict_seq['A']) >= 50:
        seq = dict_seq['I1_3p'][-250:] + dict_seq['A'][0:50]
    else:
        seq = dict_seq['I1_3p'][-250:] + dict_seq['A'] + dict_seq['I2_5p'][0:50 - len(dict_seq['A'])]
    unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.2', event_type)
    for i in xrange(3):
        start, end = 70 * i, 70 * (i + 1)
        dict_rnafold_pu_scores['SecStr.max{}_{}.I1_3p'.format(start + 1, end)] = np.max(unpair_prob_list[250 - end:250 - start])
        dict_rnafold_pu_scores['SecStr.avg{}_{}.I1_3p'.format(start + 1, end)] = np.mean(unpair_prob_list[250 - end:250 - start])
    ## junction region for 3' splice site longer (20nt)
    dict_rnafold_pu_scores['SecStr.maxJunc.I1_3p'] = np.max(unpair_prob_list[230:250])
    dict_rnafold_pu_scores['SecStr.avgJunc.I1_3p'] = np.mean(unpair_prob_list[230:250])
    return dict_rnafold_pu_scores

def build_feature(input_file, event_type, output_file, mutate):
    logger.info('building features for {}'.format(event_type))
    iDARTStmp_dir = os.path.dirname(os.path.abspath(input_file)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    splice_feature_fn = resource_filename('iDARTS.resources.features', 'SpliceCode_feature_list_{}.txt'.format(event_type))
    dict_chrom_phatsConsScore = load_phatsConsScore()
    dict_alu_counts = load_alu_counts(iDARTStmp_dir + 'event.{}.alu'.format(event_type))
    dict_ESE_matrix = load_ESE_matrix_PSSM()
    splice_feature_list = load_splice_feature_list(splice_feature_fn)
    acceptor_cnn, donor_cnn = load_cnn_splice_predictor()
    esr_matrix = load_esr_matrix()
    rosenberg_matrix = load_rosenberg_matrix()
    dict_ISE_ISS_index = {}
    dict_ShortSeq_index = {}
    dict_Yeo_cluster_index = {}
    for index, feature in enumerate(splice_feature_list):
        if feature.startswith('YeoClust'):
            dict_Yeo_cluster_index[feature] = index
        if feature.startswith('ISE.') or feature.startswith('ISS.'):
            dict_ISE_ISS_index[feature] = index
        if feature.startswith('ShortSeq.'):
            dict_ShortSeq_index[feature] = index
    ShortSeq_List = [x.split('.')[1] for x in dict_ShortSeq_index]

    file_num_lines = get_file_num_lines(input_file)
    annotation_fp = open(input_file)
    header = annotation_fp.readline().strip()
    fw = open(output_file, 'w')
    fw.write(header + '\t' + '\t'.join(splice_feature_list) + '\n')

    for annotation in tqdm(annotation_fp, total = file_num_lines - 1, desc = 'SE'):
        if len(annotation.strip()) == 0:
            continue
        annotation_list = annotation.strip().split('\t')
        # get sequence for different regions
        C1, I1_5p, I1_3p, A, I2_5p, I2_3p, C2 = get_region_coordinate(annotation.strip())
        chrom = annotation_list[CHROM_INDEX]
        strand = annotation_list[STRAND_INDEX]
        mutate_pos_list, alt_base_list = [], []
        if mutate == 'True':
            mutate_pos_list = annotation_list[POS0BASE_INDEX].split(',')
            alt_base_list = annotation_list[ALT_INDEX].split(',')
        C1_seq = seq_fetch(chrom, C1[0], C1[1], strand, mutate_pos_list, alt_base_list)
        I1_5p_seq = seq_fetch(chrom, I1_5p[0], I1_5p[1], strand, mutate_pos_list, alt_base_list)
        I1_3p_seq = seq_fetch(chrom, I1_3p[0], I1_3p[1], strand, mutate_pos_list, alt_base_list)
        A_seq = seq_fetch(chrom, A[0], A[1], strand, mutate_pos_list, alt_base_list)
        I2_5p_seq = seq_fetch(chrom, I2_5p[0], I2_5p[1], strand, mutate_pos_list, alt_base_list)
        I2_3p_seq = seq_fetch(chrom, I2_3p[0], I2_3p[1], strand, mutate_pos_list, alt_base_list)
        C2_seq = seq_fetch(chrom, C2[0], C2[1], strand, mutate_pos_list, alt_base_list)
        cnn_splice_window_size = 200 # CNN splice predictor window size
        if strand == '+':
            donor_seq_I1 = seq_fetch(chrom, C1[1] - cnn_splice_window_size, C1[1] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            donor_seq_I2 = seq_fetch(chrom, A[1] - cnn_splice_window_size, A[1] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            acceptor_seq_I1 = seq_fetch(chrom, A[0] - 2 - cnn_splice_window_size, A[0] - 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            acceptor_seq_I2 = seq_fetch(chrom, C2[0] - 2 - cnn_splice_window_size, C2[0] - 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
        if strand == '-':
            donor_seq_I1 = seq_fetch(chrom, C1[0] - cnn_splice_window_size, C1[0] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            donor_seq_I2 = seq_fetch(chrom, A[0] - cnn_splice_window_size, A[0] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            acceptor_seq_I1 = seq_fetch(chrom, A[1] + 2 - cnn_splice_window_size, A[1] + 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            acceptor_seq_I2 = seq_fetch(chrom, C1[1] + 2 - cnn_splice_window_size, C1[1] + 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
        A_seq_matrix = one_hot_encode(A_seq)
        C1_seq_matrix = one_hot_encode(C1_seq)
        C2_seq_matrix = one_hot_encode(C2_seq)
        I1_3p_seq_matrix = one_hot_encode(I1_3p_seq)
        I1_5p_seq_matrix = one_hot_encode(I1_5p_seq)
        I2_3p_seq_matrix = one_hot_encode(I2_3p_seq)
        I2_5p_seq_matrix = one_hot_encode(I2_5p_seq)
        dict_seq2matrix = {}
        dict_seq2matrix['C1'] = C1_seq_matrix
        dict_seq2matrix['A'] = A_seq_matrix
        dict_seq2matrix['C2'] = C2_seq_matrix
        dict_seq2matrix['I1_3p'] = I1_3p_seq_matrix
        dict_seq2matrix['I1_5p'] = I1_5p_seq_matrix
        dict_seq2matrix['I2_3p'] = I2_3p_seq_matrix
        dict_seq2matrix['I2_5p'] = I2_5p_seq_matrix
        dict_seq = {}
        dict_seq['C1'] = C1_seq
        dict_seq['A'] = A_seq
        dict_seq['C2'] = C2_seq
        dict_seq['I1_3p'] = I1_3p_seq
        dict_seq['I1_5p'] = I1_5p_seq
        dict_seq['I2_3p'] = I2_3p_seq
        dict_seq['I2_5p'] = I2_5p_seq
        ## get exon length features
        dict_seq_feature2value = {}
        if annotation_list[STRAND_INDEX] == '+':
            length_C1 = int(annotation_list[UPSTREAM_EE_INDEX]) - int(annotation_list[UPSTREAM_ES_INDEX])
            length_A = int(annotation_list[EXON_END_INDEX]) - int(annotation_list[EXON_START_INDEX])
            length_C2 = int(annotation_list[DOWNSTREAM_EE_INDEX]) - int(annotation_list[DOWNSTREAM_ES_INDEX])
            length_I1 = int(annotation_list[EXON_START_INDEX]) - int(annotation_list[UPSTREAM_EE_INDEX])
            length_I2 = int(annotation_list[DOWNSTREAM_ES_INDEX]) - int(annotation_list[EXON_END_INDEX])
        else:
            length_C1 = int(annotation_list[DOWNSTREAM_EE_INDEX]) - int(annotation_list[DOWNSTREAM_ES_INDEX])
            length_A = int(annotation_list[EXON_END_INDEX]) - int(annotation_list[EXON_START_INDEX])
            length_C2 = int(annotation_list[UPSTREAM_EE_INDEX]) - int(annotation_list[UPSTREAM_ES_INDEX])
            length_I1 = int(annotation_list[DOWNSTREAM_ES_INDEX]) - int(annotation_list[EXON_END_INDEX])
            length_I2 = int(annotation_list[EXON_START_INDEX]) - int(annotation_list[UPSTREAM_EE_INDEX])
        dict_seq_len = {}
        dict_seq_len['C1'] = length_C1
        dict_seq_len['A'] = length_A
        dict_seq_len['C2'] = length_C2
        dict_seq_feature2value['LogLen.C1'] = np.log10(length_C1)
        dict_seq_feature2value['LogLen.A'] = np.log10(length_A)
        dict_seq_feature2value['LogLen.C2'] = np.log10(length_C2)
        dict_seq_feature2value['LogLen.I1'] = np.log10(length_I1)
        dict_seq_feature2value['LogLen.I2'] = np.log10(length_I2)
        dict_seq_feature2value['LogLenRatio.A_I1'] = np.log10((length_A + 0.0) / length_I1)
        dict_seq_feature2value['LogLenRatio.A_I2'] = np.log10((length_A + 0.0) / length_I2)
        dict_seq_feature2value['LogLenRatio.I1_I2'] = np.log10((length_I1 + 0.0) / length_I2)
        dict_seq_feature2value['Translatable.C1'] = translatable(C1_seq)
        dict_seq_feature2value['Translatable.C1C2'] = translatable(C1_seq + C2_seq)
        dict_seq_feature2value['Translatable.C1A'] = translatable(C1_seq + A_seq)
        dict_seq_feature2value['Translatable.C1AC2'] = translatable(C1_seq + A_seq + C2_seq)
        if len(A_seq) % 3 == 0:
            dict_seq_feature2value['FrameShift.A'] = 0
        else:
            dict_seq_feature2value['FrameShift.A'] = 1
        ### Conservation score
        I1_5p_conservation, I1_3p_conservation, I2_5p_conservation, I2_3p_conservation = get_regional_conservation_score(annotation_list, dict_chrom_phatsConsScore)
        dict_region2conservation = {}
        dict_region2conservation['I1_5p'] = I1_5p_conservation
        dict_region2conservation['I1_3p'] = I1_3p_conservation
        dict_region2conservation['I2_5p'] = I2_5p_conservation
        dict_region2conservation['I2_3p'] = I2_3p_conservation
        dict_seq_feature2value['Cons.MeanP1_100.I1_3p'] = np.mean(I1_3p_conservation[-100:])
        dict_seq_feature2value['Cons.MeanP1_100.I1_5p'] = np.mean(I1_5p_conservation[0:100])
        dict_seq_feature2value['Cons.MeanP1_100.I2_3p'] = np.mean(I2_3p_conservation[-100:])
        dict_seq_feature2value['Cons.MeanP1_100.I2_5p'] = np.mean(I2_5p_conservation[0:100])
        JunctionAvgP1_100_I1_3p = np.mean(I1_3p_conservation[-2:]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.I1_3p']) + 1e-4)
        JunctionAvgP1_100_I2_5p = np.mean(I2_5p_conservation[0:2]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.I2_5p']) + 1e-4)
        JunctionAvgP1_100_I1_3p = np.sqrt(JunctionAvgP1_100_I1_3p)
        JunctionAvgP1_100_I2_5p = np.sqrt(JunctionAvgP1_100_I2_5p)
        ## get splice AG/GT position AltAGpos AltGTpos
        if 'AG' in I1_3p_seq:
            dict_seq_feature2value['AltAGpos'] = np.log2(300 - I1_3p_seq.rindex('AG')) - 1
        else:
            dict_seq_feature2value['AltAGpos'] = np.log2(300) - 1
        if 'GT' in I2_5p_seq:
            dict_seq_feature2value['AltGTpos'] = np.log2(I2_5p_seq.index('GT') + 2) - 1
        else:
            dict_seq_feature2value['AltGTpos'] = np.log2(300) - 1
        ### get splice site strength by CNN splicing predictor
        donor_score_I1 = round(donor_cnn.predict(np.asarray([one_hot_encode(donor_seq_I1)])).reshape(-1,), 3)
        acceptor_score_I1 = round(acceptor_cnn.predict(np.asarray([one_hot_encode(acceptor_seq_I1)])).reshape(-1,), 3)
        donor_score_I2 = round(donor_cnn.predict(np.asarray([one_hot_encode(donor_seq_I2)])).reshape(-1,), 3)
        acceptor_score_I2 = round(acceptor_cnn.predict(np.asarray([one_hot_encode(acceptor_seq_I2)])).reshape(-1,), 3)

        dict_seq_feature2value['donor.CNN_splice.I1'] = donor_score_I1
        dict_seq_feature2value['acceptor.CNN_splice.I1'] = acceptor_score_I1
        dict_seq_feature2value['donor.CNN_splice.I2'] = donor_score_I2
        dict_seq_feature2value['acceptor.CNN_splice.I2'] = acceptor_score_I2
        ## get acceptor and donor probablity
        dict_seq_feature2value['acceptor.CNN_splice.weighted.JunctionAvgP1_100.I1_3p'] = acceptor_score_I1 * JunctionAvgP1_100_I1_3p
        dict_seq_feature2value['donor.CNN_splice.weighted.JunctionAvgP1_100.I2_5p'] = donor_score_I2 * JunctionAvgP1_100_I2_5p
        # Alu annotation 300bp, for long, 2000bp.
        dict_alu_feature_score = get_alu_feature_score(annotation_list, dict_alu_counts)
        dict_seq_feature2value.update(dict_alu_feature_score)
        # ESS ESE Matrix score
        ESS_ESE_RBP_LIST = sorted(dict_ESE_matrix.keys())
        for region in ['C1', 'A', 'C2']:
            seq_matrix = dict_seq2matrix[region]
            score_list = get_ESE_matrix_score(seq_matrix, dict_ESE_matrix)
            for score, rbp in zip(score_list, ESS_ESE_RBP_LIST):
                dict_seq_feature2value['PSSM.{}.{}'.format(rbp, region)] = score
        # Rosenberg and ESR score
        dict_rosenberg_mean_score = get_rosenberg_mean_score(A_seq, rosenberg_matrix, ['A3SS_R2', 'A5SS_R1'])
        dict_A3SS_R1_I1_3p = get_rosenberg_mean_score(I1_3p_seq, rosenberg_matrix, ['A3SS_R1'])
        dict_A5SS_R2_I2_5p = get_rosenberg_mean_score(I2_5p_seq, rosenberg_matrix, ['A5SS_R2'])
        dict_seq_feature2value['A3SS_R2.A'] = dict_rosenberg_mean_score['A3SS_R2']
        dict_seq_feature2value['A5SS_R1.A'] = dict_rosenberg_mean_score['A5SS_R1']
        dict_seq_feature2value['A5SS_R2.I2_5p'] = dict_A5SS_R2_I2_5p['A5SS_R2']
        dict_seq_feature2value['A3SS_R1.I1_3p'] = dict_A3SS_R1_I1_3p['A3SS_R1']
        dict_seq_feature2value['ESR.A'] = get_esr_mean_score(A_seq, esr_matrix)
        # ISS and ISE counts
        ISS_ISE_Kmers = [x.split('.')[-1] for x in dict_ISE_ISS_index]
        intron_upstream_seq = I1_3p_seq[-100:-15]
        intron_downstream_seq = I2_5p_seq[15:100]
        seq_len = 100 - 15
        motif_len = 6 # ISS ISE kmer length
        split_seq_list = [intron_upstream_seq[i:i+motif_len] for i in xrange(0, seq_len - motif_len + 1)]
        kmer_count_intron_upstream = Counter(split_seq_list)
        split_seq_list = [intron_downstream_seq[i:i+motif_len] for i in xrange(0, seq_len - motif_len + 1)]
        kmer_count_intron_downstream = Counter(split_seq_list)
        for ISE_ISS_feature in dict_ISE_ISS_index:
            if 'upstream' in ISE_ISS_feature:
                kmer = ISE_ISS_feature.split('.')[-1]
                if kmer in kmer_count_intron_upstream:
                    dict_seq_feature2value[ISE_ISS_feature] = kmer_count_intron_upstream[kmer]
                else:
                    dict_seq_feature2value[ISE_ISS_feature] = 0
            if 'downstream' in ISE_ISS_feature:
                kmer = ISE_ISS_feature.split('.')[-1]
                if kmer in kmer_count_intron_downstream:
                    dict_seq_feature2value[ISE_ISS_feature] = kmer_count_intron_downstream[kmer]
                else:
                    dict_seq_feature2value[ISE_ISS_feature] = 0
        for region in ['C1', 'A', 'C2', 'I1_3p', 'I1_5p', 'I2_3p', 'I2_5p']:
            seq = dict_seq[region]
            seq_len = len(seq)
            kmer_1_count = Counter([seq[i:i+1] for i in xrange(seq_len - 1 + 1)])
            kmer_2_count = Counter([seq[i:i+2] for i in xrange(seq_len - 2 + 1)])
            kmer_3_count = Counter([seq[i:i+3] for i in xrange(seq_len - 3 + 1)])
            kmer_count = {}
            kmer_count.update(kmer_1_count)
            kmer_count.update(kmer_2_count)
            kmer_count.update(kmer_3_count)
            for short_seq in ShortSeq_List:
                if 'ShortSeq.{}.{}'.format(short_seq, region) not in dict_ShortSeq_index:
                    continue
                if short_seq in kmer_count:
                    dict_seq_feature2value['ShortSeq.{}.{}'.format(short_seq, region)] = kmer_count[short_seq] / (seq_len - len(short_seq) + 1.0)
                else:
                    dict_seq_feature2value['ShortSeq.{}.{}'.format(short_seq, region)] = 0.0
        ## NOVA annotation
        nova_motif = '[TC]CA[TC]'
        nova_motif_len = 4
        dict_nova_region = get_nova_binding_score(annotation_list, dict_seq_len, dict_seq, dict_region2conservation, dict_chrom_phatsConsScore)
        for region in dict_nova_region:
            seq = dict_nova_region[region]['seq']
            region_conservation_scores = dict_nova_region[region]['conservation']
            pos_list = [m for m in re.finditer(r'{0}'.format(nova_motif), seq, overlapped=True)]
            count = len(pos_list)
            conserve_score = 0.0
            if len(pos_list) > 0:
                for pos in pos_list:
                    pos = pos.start(0)
                    conserve_score += np.mean(region_conservation_scores[pos:pos + nova_motif_len])
            dict_seq_feature2value['Nova.{0}.withCons'.format(region)] = conserve_score
            dict_seq_feature2value['Nova.{0}'.format(region)] = count
        ## Yeo Clust
        dict_seq2kmer_loc = {}
        for kmer_len in [5, 6, 7]:
            dict_seq2kmer_loc[kmer_len] = {}
            for region in ['I1_3p', 'I2_5p']:
                dict_seq2kmer_loc[kmer_len][region] = defaultdict(lambda:[])
                seq = dict_seq[region]
                for i in xrange(0, 300 - kmer_len + 1):
                    dict_seq2kmer_loc[kmer_len][region][seq[i:i+kmer_len]].append(i)
        for yeo_cluster in dict_Yeo_cluster_index:
            sp = yeo_cluster.split('.')
            motif = sp[2]
            region = sp[3]
            kmer_len = len(motif)
            region_conservation_scores = dict_region2conservation[region]
            pos_list = dict_seq2kmer_loc[kmer_len][region][motif]
            conserve_score = 0.0
            if len(pos_list) > 0:
                for pos in pos_list:
                    conserve_score += np.mean(region_conservation_scores[pos:pos+kmer_len])
            dict_seq_feature2value[yeo_cluster] = conserve_score
        ## RNA secondary structure
        dict_rnafold_pu_scores = get_rnafold_predict_pu_score(dict_seq, iDARTStmp_dir + event_type + '.sec', event_type)
        dict_seq_feature2value.update(dict_rnafold_pu_scores)

        out_line_str = ''
        for feature in splice_feature_list:
            out_line_str += str(round(dict_seq_feature2value[feature], 5)) + '\t'
        fw.write(annotation.strip() + '\t' + out_line_str.strip() + '\n')
    fw.close()

def main(args):
    fetch_alu_code(args)
    build_feature(args.input, args.event_type, args.output, args.mutate)
    iDARTStmp_dir = os.path.dirname(os.path.abspath(args.input)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    logger.info('Deleting tmp files')
    remove_file(iDARTStmp_dir)
