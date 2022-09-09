# -*- coding: UTF-8 -*-

"""
iDARTS - build_feature_ASS
Implements a splicing feature build module for 
alternative 5' splice sites events or 
alternative 3' splice sites features
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
LONGEXONSTART = 5
LONGEXONEND = 6
SHORTEXONSTART = 7
SHORTEXONEND = 8
FLANKINGEXONSTART = 9
FLANKINGEXONEND = 10
POS0BASE_INDEX = 12
REF_INDEX = 13
ALT_INDEX = 14
TEMP_FILE_NAME = next(tempfile._get_candidate_names())

def get_region_coordinate(rmats_input_list, event_type):
    sp = rmats_input_list
    C_L = [int(sp[LONGEXONSTART]), int(sp[LONGEXONEND])]
    C_S = [int(sp[SHORTEXONSTART]), int(sp[SHORTEXONEND])]
    C_F = [int(sp[FLANKINGEXONSTART]), int(sp[FLANKINGEXONEND])]
    if sp[STRAND_INDEX] == '-':
        if event_type == 'A3SS':
            C_L_300 = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 300]
            C_S_300 = [int(sp[SHORTEXONEND]), int(sp[SHORTEXONEND]) + 300]
            C_F_300 = [int(sp[FLANKINGEXONSTART]) - 300, int(sp[FLANKINGEXONSTART])]
            C_L_long = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 2000]
            C_S_long = [int(sp[SHORTEXONEND]), int(sp[SHORTEXONEND]) + 2000]
            C_F_long = [int(sp[FLANKINGEXONSTART]) - 2000, int(sp[FLANKINGEXONSTART])]
        else:
            C_L_300 = [int(sp[LONGEXONSTART]) - 300, int(sp[LONGEXONSTART])]
            C_S_300 = [int(sp[SHORTEXONSTART]) - 300, int(sp[SHORTEXONSTART])]
            C_F_300 = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 300]
            C_L_long = [int(sp[LONGEXONSTART]) - 2000, int(sp[LONGEXONSTART])]
            C_S_long = [int(sp[SHORTEXONSTART]) - 2000, int(sp[SHORTEXONSTART])]
            C_F_long = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 2000]
    else:
        if event_type == 'A5SS':
            C_L_300 = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 300]
            C_S_300 = [int(sp[SHORTEXONEND]), int(sp[SHORTEXONEND]) + 300]
            C_F_300 = [int(sp[FLANKINGEXONSTART]) - 300, int(sp[FLANKINGEXONSTART])]
            C_L_long = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 2000]
            C_S_long = [int(sp[SHORTEXONEND]), int(sp[SHORTEXONEND]) + 2000]
            C_F_long = [int(sp[FLANKINGEXONSTART]) - 2000, int(sp[FLANKINGEXONSTART])]
        else:
            C_L_300 = [int(sp[LONGEXONSTART]) - 300, int(sp[LONGEXONSTART])]
            C_S_300 = [int(sp[SHORTEXONSTART]) - 300, int(sp[SHORTEXONSTART])]
            C_F_300 = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 300]
            C_L_long = [int(sp[LONGEXONSTART]) - 2000, int(sp[LONGEXONSTART])]
            C_S_long = [int(sp[SHORTEXONSTART]) - 2000, int(sp[SHORTEXONSTART])]
            C_F_long = [int(sp[LONGEXONEND]), int(sp[LONGEXONEND]) + 2000]
    return C_L, C_S, C_F, C_L_300, C_S_300, C_F_300, C_L_long, C_S_long, C_F_long

def fetch_seq_alu_code(args): # the input file (rMATS format)
    input_file = args.input
    event_type = args.event_type
    iDARTStmp_dir = os.path.dirname(os.path.abspath(input_file)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    makedirs(iDARTStmp_dir)
    tmpAluBed = iDARTStmp_dir + 'tmp.{}.alu.bed'.format(event_type)
    alufw = open(tmpAluBed, 'w')
    fw = open(iDARTStmp_dir + 'event.{}.seq'.format(event_type), 'w')
    fw.write('ID\tC_L\tC_S\tC_F\tC_L_300\tC_S_300\tC_F_300\n')
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
        C_L, C_S, C_F, C_L_300, C_S_300, C_F_300, C_L_long, C_S_long, C_F_long = get_region_coordinate(sp, event_type)
        if args.mutate == 'True':
            mutate_pos_list = sp[POS0BASE_INDEX].split(',')
            ref_base_list = sp[REF_INDEX].split(',')
            alt_base_list = sp[ALT_INDEX].split(',')
            C_L_seq = seq_fetch(chrom, C_L[0], C_L[1], strand, mutate_pos_list, alt_base_list)
            C_S_seq = seq_fetch(chrom, C_S[0], C_S[1], strand, mutate_pos_list, alt_base_list)
            C_F_seq = seq_fetch(chrom, C_F[0], C_F[1], strand, mutate_pos_list, alt_base_list)
            C_L_300_seq = seq_fetch(chrom, C_L_300[0], C_L_300[1], strand, mutate_pos_list, alt_base_list)
            C_S_300_seq = seq_fetch(chrom, C_S_300[0], C_S_300[1], strand, mutate_pos_list, alt_base_list)
            C_F_300_seq = seq_fetch(chrom, C_F_300[0], C_F_300[1], strand, mutate_pos_list, alt_base_list)
        else:
            C_L_seq = seq_fetch(chrom, C_L[0], C_L[1], strand)
            C_S_seq = seq_fetch(chrom, C_S[0], C_S[1], strand)
            C_F_seq = seq_fetch(chrom, C_F[0], C_F[1], strand)
            C_L_300_seq = seq_fetch(chrom, C_L_300[0], C_L_300[1], strand)
            C_S_300_seq = seq_fetch(chrom, C_S_300[0], C_S_300[1], strand)
            C_F_300_seq = seq_fetch(chrom, C_F_300[0], C_F_300[1], strand)
        fw.write(sp[0] + '\t' + '\t'.join([C_L_seq, C_S_seq, C_F_seq, C_L_300_seq, C_S_300_seq, C_F_300_seq]) + '\n')
        # Alu repeat writer...
        strInput = '{0}\t{1}\t{2}\tforward\t1\t+\n'
        alufw.write(strInput.format(sp[CHROM_INDEX], C_L_300[0], C_L_300[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_S_300[0], C_S_300[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_F_300[0], C_F_300[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_L_long[0], C_L_long[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_S_long[0], C_S_long[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_F_long[0], C_F_long[1]))
        strInput = '{0}\t{1}\t{2}\treverse\t1\t-\n'
        alufw.write(strInput.format(sp[CHROM_INDEX], C_L_300[0], C_L_300[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_S_300[0], C_S_300[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_F_300[0], C_F_300[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_L_long[0], C_L_long[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_S_long[0], C_S_long[1]))
        alufw.write(strInput.format(sp[CHROM_INDEX], C_F_long[0], C_F_long[1]))
    fp.close()
    alufw.close()
    fw.close()
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

def get_alu_feature_score(annotation_list, dict_alu_counts, event_type):
    C_L, C_S, C_F, C_L_300, C_S_300, C_F_300, C_L_long, C_S_long, C_F_long = get_region_coordinate(annotation_list, event_type)
    dict_alu_feature_score = {}
    chrom = annotation_list[CHROM_INDEX]
    dict_alu_region = {}
    dict_alu_region['C_S'] = C_S
    dict_alu_region['C_L'] = C_L
    dict_alu_region['C_F'] = C_F
    dict_alu_region['C_S_long'] = C_S_long
    dict_alu_region['C_L_long'] = C_L_long
    dict_alu_region['C_F_long'] = C_F_long
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

def get_conservation_score(dict_chrom_phatsConsScore, chrom, start, end):
    conservation_scores = dict_chrom_phatsConsScore[chrom].values(chrom, start, end)
    return [x if x == x else 0 for x in conservation_scores]

def get_regional_conservation_score(annotation_list, dict_chrom_phatsConsScore, event_type):
    strand = annotation_list[STRAND_INDEX]
    chrom = annotation_list[CHROM_INDEX]
    flanking_exon_start, flanking_exon_end = int(annotation_list[FLANKINGEXONSTART]), int(annotation_list[FLANKINGEXONEND])
    long_exon_start, long_exon_end = int(annotation_list[LONGEXONSTART]), int(annotation_list[LONGEXONEND])
    short_exon_start, short_exon_end = int(annotation_list[SHORTEXONSTART]), int(annotation_list[SHORTEXONEND])
    if strand == '+':
        if event_type == 'A3SS':
            C_F_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, flanking_exon_end, flanking_exon_end + 300)
            C_S_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, short_exon_start - 300, short_exon_start)
            C_L_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, long_exon_start - 300, long_exon_start)
        if event_type == 'A5SS':
            C_F_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, flanking_exon_start - 300, flanking_exon_start)
            C_S_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, short_exon_end, short_exon_end + 300)
            C_L_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, long_exon_end, long_exon_end + 300)
    if strand == '-':
        if event_type == 'A3SS':
            C_F_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, flanking_exon_start - 300, flanking_exon_start)[::-1]
            C_S_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, short_exon_end, short_exon_end + 300)[::-1]
            C_L_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, long_exon_end, long_exon_end + 300)[::-1]
        if event_type == 'A5SS':
            C_F_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, flanking_exon_end, flanking_exon_end + 300)[::-1]
            C_S_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, short_exon_start - 300, short_exon_start)[::-1]
            C_L_300_conservation = get_conservation_score(dict_chrom_phatsConsScore, chrom, long_exon_start - 300, long_exon_start)[::-1]
    return C_F_300_conservation, C_S_300_conservation, C_L_300_conservation

def parse_rnafold_predict_pu_score(unpair_prob_list, junction_start, region, intron_region_direction):
    dict_rnafold_pu_scores = {}
    if intron_region_direction == 'right':
        for i in xrange(3):
            start, end = 70 * i, 70 + 70 * i
            dict_rnafold_pu_scores['SecStr.max{}_{}.{}'.format(start + 1, end, region)] = np.max(unpair_prob_list[start + junction_start:end + junction_start])
            dict_rnafold_pu_scores['SecStr.avg{}_{}.{}'.format(start + 1, end, region)] = np.mean(unpair_prob_list[start + junction_start:end + junction_start])
        dict_rnafold_pu_scores['SecStr.maxJunc.{}'.format(region)] = np.max(unpair_prob_list[junction_start:junction_start + 10])
        dict_rnafold_pu_scores['SecStr.avgJunc.{}'.format(region)] = np.mean(unpair_prob_list[junction_start:junction_start + 10])
    if intron_region_direction == 'left':
        for i in xrange(3):
            start, end = 70 * i, 70 + 70 * i
            dict_rnafold_pu_scores['SecStr.max{}_{}.{}'.format(start + 1, end, region)] = np.max(unpair_prob_list[250 - end:250 - start])
            dict_rnafold_pu_scores['SecStr.avg{}_{}.{}'.format(start + 1, end, region)] = np.mean(unpair_prob_list[250 - end:250 - start])
        dict_rnafold_pu_scores['SecStr.maxJunc.{}'.format(region)] = np.max(unpair_prob_list[junction_start:junction_start + 10])
        dict_rnafold_pu_scores['SecStr.avgJunc.{}'.format(region)] = np.mean(unpair_prob_list[junction_start:junction_start + 10])
    return dict_rnafold_pu_scores

def get_rnafold_predict_pu_score(dict_seq, tmp_file, event_type):
    dict_rnafold_pu_scores = {}
    if event_type == 'A3SS':
        seq = dict_seq['C_F'][-min(50, len(dict_seq['C_F']))] + dict_seq['C_F_300'][0:250]
        unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.1', event_type)
        dict_rnafold_pu_scores.update(parse_rnafold_predict_pu_score(unpair_prob_list, min(50, len(dict_seq['C_F'])), 'C_F', 'right'))

        seq = dict_seq['C_L_300'][-250:] + dict_seq['C_L'][0:min(50, len(dict_seq['C_L']))]
        unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.2', event_type)
        dict_rnafold_pu_scores.update(parse_rnafold_predict_pu_score(unpair_prob_list, 240, 'C_L', 'left'))

        seq = dict_seq['C_S_300'][-250:] + dict_seq['C_S'][0:min(50, len(dict_seq['C_S']))]
        unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.3', event_type)
        dict_rnafold_pu_scores.update(parse_rnafold_predict_pu_score(unpair_prob_list, 240, 'C_S', 'left'))

    if event_type == 'A5SS':
        seq = dict_seq['C_F_300'][-250:] + dict_seq['C_F'][0:min(50, len(dict_seq['C_F']))]
        unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.1', event_type)
        dict_rnafold_pu_scores.update(parse_rnafold_predict_pu_score(unpair_prob_list, 240, 'C_F', 'left'))

        seq = dict_seq['C_L'][-min(50, len(dict_seq['C_L'])):] + dict_seq['C_L_300'][0:250]
        unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.2', event_type)
        dict_rnafold_pu_scores.update(parse_rnafold_predict_pu_score(unpair_prob_list, min(50, len(dict_seq['C_L'])), 'C_L', 'right'))

        seq = dict_seq['C_S'][-min(50, len(dict_seq['C_S'])):] + dict_seq['C_S_300'][0:250]
        unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.3', event_type)
        dict_rnafold_pu_scores.update(parse_rnafold_predict_pu_score(unpair_prob_list, min(50, len(dict_seq['C_S'])), 'C_S', 'right'))
    return dict_rnafold_pu_scores

def build_feature(input_file, event_type, output_file, mutate):
    logger.info('building features for {}'.format(event_type))
    iDARTStmp_dir = os.path.dirname(os.path.abspath(input_file)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    splice_feature_fn = resource_filename('iDARTS.resources.features', 'SpliceCode_feature_list_{}.txt'.format(event_type))
    dict_chrom_phatsConsScore = load_phatsConsScore()
    dict_ESE_matrix = load_ESE_matrix_PSSM()
    dict_rna_binding_pssm = load_rna_binding_motif_pssm()
    dict_rna_binding_motif = load_rna_binding_motif()
    splice_feature_list = load_splice_feature_list(splice_feature_fn)
    acceptor_cnn, donor_cnn = load_cnn_splice_predictor()
    esr_matrix = load_esr_matrix()
    rosenberg_matrix = load_rosenberg_matrix()
    dict_alu_counts = load_alu_counts(iDARTStmp_dir + 'event.{}.alu'.format(event_type))

    dict_ShortSeq_index = {}
    dict_Yeo_cluster_index = {}
    for index, feature in enumerate(splice_feature_list):
        if feature.startswith('YeoClust'):
            dict_Yeo_cluster_index[feature] = index
        if feature.startswith('ShortSeq.'):
            dict_ShortSeq_index[feature] = index
    ShortSeq_List = [x.split('.')[1] for x in dict_ShortSeq_index]

    file_num_lines = get_file_num_lines(input_file)
    annotation_fp = open(input_file)
    header = annotation_fp.readline().strip()
    fw = open(output_file, 'w')
    fw.write(header + '\t' + '\t'.join(splice_feature_list) + '\n')

    for annotation in tqdm(annotation_fp, total = file_num_lines - 1):
        if len(annotation.strip()) == 0:
            continue
        annotation_list = annotation.strip().split('\t')
        C_L, C_S, C_F, C_L_300, C_S_300, C_F_300, C_L_long, C_S_long, C_F_long = get_region_coordinate(annotation_list, event_type)
        chrom = annotation_list[CHROM_INDEX]
        strand = annotation_list[STRAND_INDEX]
        mutate_pos_list, alt_base_list = [], []
        if mutate == 'True':
            mutate_pos_list = annotation_list[POS0BASE_INDEX].split(',')
            alt_base_list = annotation_list[ALT_INDEX].split(',')
        C_L_seq = seq_fetch(chrom, C_L[0], C_L[1], strand, mutate_pos_list, alt_base_list)
        C_S_seq = seq_fetch(chrom, C_S[0], C_S[1], strand, mutate_pos_list, alt_base_list)
        C_F_seq = seq_fetch(chrom, C_F[0], C_F[1], strand, mutate_pos_list, alt_base_list)
        C_L_300_seq = seq_fetch(chrom, C_L_300[0], C_L_300[1], strand, mutate_pos_list, alt_base_list)
        C_S_300_seq = seq_fetch(chrom, C_S_300[0], C_S_300[1], strand, mutate_pos_list, alt_base_list)
        C_F_300_seq = seq_fetch(chrom, C_F_300[0], C_F_300[1], strand, mutate_pos_list, alt_base_list)
        C_L_seq_matrix = one_hot_encode(C_L_seq)
        C_S_seq_matrix = one_hot_encode(C_S_seq)
        C_F_seq_matrix = one_hot_encode(C_F_seq)
        C_L_300_seq_matrix = one_hot_encode(C_L_300_seq)
        C_S_300_seq_matrix = one_hot_encode(C_S_300_seq)
        C_F_300_seq_matrix = one_hot_encode(C_F_300_seq)
        dict_seq2matrix = {}
        dict_seq2matrix['C_L'] = C_L_seq_matrix
        dict_seq2matrix['C_S'] = C_S_seq_matrix
        dict_seq2matrix['C_F'] = C_F_seq_matrix
        dict_seq2matrix['C_L_300'] = C_L_300_seq_matrix
        dict_seq2matrix['C_S_300'] = C_S_300_seq_matrix
        dict_seq2matrix['C_F_300'] = C_F_300_seq_matrix
        dict_seq = {}
        dict_seq['C_L'] = C_L_seq
        dict_seq['C_S'] = C_S_seq
        dict_seq['C_F'] = C_F_seq
        dict_seq['C_L_300'] = C_L_300_seq
        dict_seq['C_S_300'] = C_S_300_seq
        dict_seq['C_F_300'] = C_F_300_seq
        ## get exon length features
        length_C_L = int(annotation_list[LONGEXONEND]) - int(annotation_list[LONGEXONSTART])
        length_C_S = int(annotation_list[SHORTEXONEND]) - int(annotation_list[SHORTEXONSTART])
        length_C_F = int(annotation_list[FLANKINGEXONEND]) - int(annotation_list[FLANKINGEXONSTART])
        if annotation_list[STRAND_INDEX] == '+':
            if event_type == 'A3SS':
                length_I_S = int(annotation_list[SHORTEXONSTART]) - int(annotation_list[FLANKINGEXONEND])
                length_I_L = int(annotation_list[LONGEXONSTART]) - int(annotation_list[FLANKINGEXONEND])
            elif event_type == 'A5SS':
                length_I_S = int(annotation_list[FLANKINGEXONSTART]) - int(annotation_list[SHORTEXONEND])
                length_I_L = int(annotation_list[FLANKINGEXONSTART]) - int(annotation_list[LONGEXONEND])
        else:
            if event_type == 'A5SS':
                length_I_S = int(annotation_list[SHORTEXONSTART]) - int(annotation_list[FLANKINGEXONEND])
                length_I_L = int(annotation_list[LONGEXONSTART]) - int(annotation_list[FLANKINGEXONEND])
            elif event_type == 'A3SS':
                length_I_S = int(annotation_list[FLANKINGEXONSTART]) - int(annotation_list[SHORTEXONEND])
                length_I_L = int(annotation_list[FLANKINGEXONSTART]) - int(annotation_list[LONGEXONEND])
        dict_seq_feature2value = {}
        dict_seq_feature2value['LogLen.C_S'] = np.log10(length_C_S)
        dict_seq_feature2value['LogLen.C_L'] = np.log10(length_C_L)
        dict_seq_feature2value['LogLen.C_F'] = np.log10(length_C_F)
        dict_seq_feature2value['LogLen.I_S'] = np.log10(length_I_S + 1.0)
        dict_seq_feature2value['LogLen.I_L'] = np.log10(length_I_L + 1.0)
        dict_seq_feature2value['LogLenRatio.I_S_C_S'] = np.log10((length_I_S + 1.0) / length_C_S)
        dict_seq_feature2value['LogLenRatio.I_L_C_L'] = np.log10((length_I_L + 1.0) / length_C_L)
        dict_seq_feature2value['Translatable.C_S'] = translatable(C_S_seq)
        dict_seq_feature2value['Translatable.C_L'] = translatable(C_L_seq)
        dict_seq_feature2value['Translatable.C_F'] = translatable(C_F_seq)
        if event_type == 'A3SS':
            dict_seq_feature2value['Translatable.C_F_C_S'] = translatable(C_F_seq + C_S_seq)
            dict_seq_feature2value['Translatable.C_F_C_L'] = translatable(C_F_seq + C_L_seq)
        elif event_type == 'A5SS':
            dict_seq_feature2value['Translatable.C_S_C_F'] = translatable(C_S_seq + C_F_seq)
            dict_seq_feature2value['Translatable.C_L_C_F'] = translatable(C_L_seq + C_F_seq)
        if len(C_S_seq) % 3 == 0:
            dict_seq_feature2value['FrameShift.C_S'] = 0
        else:
            dict_seq_feature2value['FrameShift.C_S'] = 1
        if len(C_L_seq) % 3 == 0:
            dict_seq_feature2value['FrameShift.C_L'] = 0
        else:
            dict_seq_feature2value['FrameShift.C_L'] = 1
        ## get splice AG/GT position AltAGpos AltGTpos
        if event_type == 'A3SS':
            if 'AG' in C_S_300_seq:
                dict_seq_feature2value['AltAGpos.C_S'] = len(C_S_300_seq) - C_S_300_seq.rindex('AG')
            else:
                dict_seq_feature2value['AltAGpos.C_S'] = 300
            if 'AG' in C_L_300_seq:
                dict_seq_feature2value['AltAGpos.C_L'] = len(C_L_300_seq) - C_L_300_seq.rindex('AG')
            else:
                dict_seq_feature2value['AltAGpos.C_L'] = 300
            if 'GT' in C_F_300_seq:
                dict_seq_feature2value['AltGTpos.C_F'] = C_F_300_seq.index('GT')
            else:
                dict_seq_feature2value['AltGTpos.C_F'] = 300
        if event_type == 'A5SS':
            if 'GT' in C_S_300_seq:
                dict_seq_feature2value['AltGTpos.C_S'] = C_S_300_seq.index('GT')
            else:
                dict_seq_feature2value['AltGTpos.C_S'] = 300
            if 'GT' in C_L_300_seq:
                dict_seq_feature2value['AltGTpos.C_L'] = C_L_300_seq.index('GT')
            else:
                dict_seq_feature2value['AltGTpos.C_L'] = 300
            if 'AG' in C_F_300_seq:
                dict_seq_feature2value['AltAGpos.C_F'] = len(C_F_300_seq) - C_F_300_seq.rindex('AG')
            else:
                dict_seq_feature2value['AltAGpos.C_F'] = 300
        ### Conservation score
        C_F_300_conservation, C_S_300_conservation, C_L_300_conservation = get_regional_conservation_score(annotation_list, dict_chrom_phatsConsScore, event_type)
        dict_region2conservation = {}
        dict_region2conservation['C_F_300'] = C_F_300_conservation
        dict_region2conservation['C_S_300'] = C_S_300_conservation
        dict_region2conservation['C_L_300'] = C_L_300_conservation
        if event_type == 'A3SS':
            dict_seq_feature2value['Cons.MeanP1_100.C_F'] = np.mean(C_F_300_conservation[0:100])
            dict_seq_feature2value['Cons.MeanP1_100.C_S'] = np.mean(C_S_300_conservation[-100:])
            dict_seq_feature2value['Cons.MeanP1_100.C_L'] = np.mean(C_L_300_conservation[-100:])
            JunctionAvgP1_100_C_F = np.mean(C_F_300_conservation[0:2]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.C_F']) + 1e-4)
            JunctionAvgP1_100_C_S = np.mean(C_S_300_conservation[-2:]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.C_S']) + 1e-4)
            JunctionAvgP1_100_C_L = np.mean(C_L_300_conservation[-2:]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.C_L']) + 1e-4)
        if event_type == 'A5SS':
            dict_seq_feature2value['Cons.MeanP1_100.C_F'] = np.mean(C_F_300_conservation[-100:])
            dict_seq_feature2value['Cons.MeanP1_100.C_S'] = np.mean(C_S_300_conservation[0:100])
            dict_seq_feature2value['Cons.MeanP1_100.C_L'] = np.mean(C_L_300_conservation[0:100])
            JunctionAvgP1_100_C_F = np.mean(C_F_300_conservation[-2:]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.C_F']) + 1e-4)
            JunctionAvgP1_100_C_S = np.mean(C_S_300_conservation[0:2]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.C_S']) + 1e-4)
            JunctionAvgP1_100_C_L = np.mean(C_L_300_conservation[0:2]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.C_L']) + 1e-4)

        cnn_splice_window_size = 200 # CNN splice predictor window size
        if event_type == 'A3SS':
            if strand == '+':
                donor_C_F = seq_fetch(chrom, C_F[1] - cnn_splice_window_size, C_F[1] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                acceptor_C_S = seq_fetch(chrom, C_S[0] - 2 - cnn_splice_window_size, C_S[0] - 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                acceptor_C_L = seq_fetch(chrom, C_L[0] - 2 - cnn_splice_window_size, C_L[0] - 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            if strand == '-':
                donor_C_F = seq_fetch(chrom, C_F[0] - cnn_splice_window_size, C_F[0] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                acceptor_C_S = seq_fetch(chrom, C_S[1] + 2 - cnn_splice_window_size, C_S[1] + 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                acceptor_C_L = seq_fetch(chrom, C_L[1] + 2 - cnn_splice_window_size, C_L[1] + 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            donor_C_F_score = round(donor_cnn.predict(np.asarray([one_hot_encode(donor_C_F)])).reshape(-1,), 3)
            acceptor_C_S_score = round(acceptor_cnn.predict(np.asarray([one_hot_encode(acceptor_C_S)])).reshape(-1,), 3)
            acceptor_C_L_score = round(acceptor_cnn.predict(np.asarray([one_hot_encode(acceptor_C_L)])).reshape(-1,), 3)
            ## get acceptor and donor probablity
            donor_C_F_score_prob = 2 ** donor_C_F_score / (2 ** donor_C_F_score + 1.0)
            acceptor_C_S_score_prob = 2 ** acceptor_C_S_score / (2 ** acceptor_C_S_score + 1.0)
            acceptor_C_L_score_prob = 2 ** acceptor_C_L_score / (2 ** acceptor_C_L_score + 1.0)
            dict_seq_feature2value['donor.CNN_C_F'] = donor_C_F_score
            dict_seq_feature2value['acceptor.CNN_C_S'] = acceptor_C_S_score
            dict_seq_feature2value['acceptor.CNN_C_L'] = acceptor_C_L_score
            dict_seq_feature2value['donor.weighted.JunctionAvgP1_100_C_F'] = JunctionAvgP1_100_C_F * donor_C_F_score_prob
            dict_seq_feature2value['acceptor.weighted.JunctionAvgP1_100_C_S'] = JunctionAvgP1_100_C_S * acceptor_C_S_score_prob
            dict_seq_feature2value['acceptor.weighted.JunctionAvgP1_100_C_L'] = JunctionAvgP1_100_C_L * acceptor_C_L_score_prob

        if event_type == 'A5SS':
            if strand == '+':
                acceptor_C_F = seq_fetch(chrom, C_F[0] - 2 - cnn_splice_window_size, C_F[0] - 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                donor_C_S = seq_fetch(chrom, C_S[1] - cnn_splice_window_size, C_S[1] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                donor_C_L = seq_fetch(chrom, C_L[1] - cnn_splice_window_size, C_L[1] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            if strand == '-':
                acceptor_C_F = seq_fetch(chrom, C_F[1] + 2 - cnn_splice_window_size, C_F[1] + 2 + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                donor_C_S = seq_fetch(chrom, C_S[0] - cnn_splice_window_size, C_S[0] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
                donor_C_L = seq_fetch(chrom, C_L[0] - cnn_splice_window_size, C_L[0] + cnn_splice_window_size, strand, mutate_pos_list, alt_base_list)
            acceptor_C_F_score = round(acceptor_cnn.predict(np.asarray([one_hot_encode(acceptor_C_F)])).reshape(-1,), 3)
            donor_C_S_score = round(donor_cnn.predict(np.asarray([one_hot_encode(donor_C_S)])).reshape(-1,), 3)
            donor_C_L_score = round(donor_cnn.predict(np.asarray([one_hot_encode(donor_C_L)])).reshape(-1,), 3)
            ## get acceptor and donor probablity
            acceptor_C_F_score_prob = 2 ** acceptor_C_F_score / (2 ** acceptor_C_F_score + 1.0)
            donor_C_S_score_prob = 2 ** donor_C_S_score / (2 ** donor_C_S_score + 1.0)
            donor_C_L_score_prob = 2 ** donor_C_L_score / (2 ** donor_C_L_score + 1.0)
            dict_seq_feature2value['acceptor.CNN_C_F'] = acceptor_C_F_score
            dict_seq_feature2value['donor.CNN_C_S'] = donor_C_S_score
            dict_seq_feature2value['donor.CNN_C_L'] = donor_C_L_score
            dict_seq_feature2value['acceptor.weighted.JunctionAvgP1_100_C_F'] = JunctionAvgP1_100_C_F * acceptor_C_F_score_prob
            dict_seq_feature2value['donor.weighted.JunctionAvgP1_100_C_S'] = JunctionAvgP1_100_C_S * donor_C_S_score_prob
            dict_seq_feature2value['donor.weighted.JunctionAvgP1_100_C_L'] = JunctionAvgP1_100_C_L * donor_C_L_score_prob
        # Alu annotation 300bp, for long, 2000bp.
        dict_alu_feature_score = get_alu_feature_score(annotation_list, dict_alu_counts, event_type)
        dict_seq_feature2value.update(dict_alu_feature_score)
        # ESS ESE Matrix score
        ESS_ESE_RBP_LIST = sorted(dict_ESE_matrix.keys())
        for region in ['C_F', 'C_L', 'C_S']:
            seq_matrix = dict_seq2matrix[region]
            score_list = get_ESE_matrix_score(seq_matrix, dict_ESE_matrix)
            for score, rbp in zip(score_list, ESS_ESE_RBP_LIST):
                dict_seq_feature2value['PSSM.{}.{}'.format(rbp, region)] = score
        # Rosenberg and ESR score
        dict_rosenberg_mean_score_L = get_rosenberg_mean_score(C_L_seq, rosenberg_matrix, ['A3SS_R2', 'A5SS_R1'])
        dict_rosenberg_mean_score_S = get_rosenberg_mean_score(C_S_seq, rosenberg_matrix, ['A3SS_R2', 'A5SS_R1'])
        dict_seq_feature2value['ESR.C_L'] = get_esr_mean_score(C_L_seq, esr_matrix)
        dict_seq_feature2value['ESR.C_S'] = get_esr_mean_score(C_S_seq, esr_matrix)
        dict_seq_feature2value['A3SS_R2.C_L'] = dict_rosenberg_mean_score_L['A3SS_R2']
        dict_seq_feature2value['A5SS_R1.C_L'] = dict_rosenberg_mean_score_L['A5SS_R1']
        dict_seq_feature2value['A3SS_R2.C_S'] = dict_rosenberg_mean_score_S['A3SS_R2']
        dict_seq_feature2value['A5SS_R1.C_S'] = dict_rosenberg_mean_score_S['A5SS_R1']
        # RBP binding motif PSSM profile annotation
        for region in dict_seq2matrix:
            seq_matrix = dict_seq2matrix[region]
            score_list = get_pssm_score(seq_matrix, dict_rna_binding_pssm)
            for index, binding_protein in enumerate(sorted(dict_rna_binding_pssm.keys())):
                dict_seq_feature2value[binding_protein + '.' + region] = score_list[index]
        ## ShortSeq Annotation
        for region in dict_seq:
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
        ## RNA BINDING MOTIF (NON-PSSM VERSION)
        for region in ['C_L', 'C_S', 'C_F']:
            seq = dict_seq[region]
            seq_len = len(seq)
            for binding_protein in dict_rna_binding_motif:  # RBP Binding motif annotation
                re_motif = dict_rna_binding_motif[binding_protein]
                pos_list = [m for m in re.finditer(r'{0}'.format(re_motif), seq, overlapped=True)]
                count = len(pos_list)
                if len(pos_list) > 0:
                    motif_length = pos_list[0].end(0) - pos_list[0].start(0)
                else:
                    motif_length = 0.0
                dict_seq_feature2value[binding_protein + '.' + region] = count / (seq_len - motif_length + 1.0)
        for region in ['C_L_300', 'C_S_300', 'C_L_300']:
            seq = dict_seq[region]
            region_conservation_scores = dict_region2conservation[region]
            for binding_protein in dict_rna_binding_motif:  # RBP Binding motif annotation
                count, conserve_score, motif_length = 0.0, 0.0, 0.0
                re_motif = dict_rna_binding_motif[binding_protein]
                pos_list = [m for m in re.finditer(r'{0}'.format(re_motif), seq, overlapped=True)]
                count = len(pos_list)
                if len(pos_list) > 0:
                    motif_length = pos_list[0].end(0) - pos_list[0].start(0)
                    for pos in pos_list:
                        pos = pos.start(0)
                        conserve_score += np.mean(region_conservation_scores[pos:pos + motif_length])
                dict_seq_feature2value[binding_protein + '.' + region + '.withCons'] = conserve_score 
                dict_seq_feature2value[binding_protein + '.' + region] = count
        ## Yeo Clust
        dict_seq2kmer_loc = {}
        for kmer_len in [5, 6, 7]:
            dict_seq2kmer_loc[kmer_len] = {}
            for region in ['C_L_300', 'C_S_300', 'C_F_300']:
                dict_seq2kmer_loc[kmer_len][region] = defaultdict(lambda:[])
                seq = dict_seq[region]
                for i in xrange(0, 300 - kmer_len + 1):
                    dict_seq2kmer_loc[kmer_len][region][seq[i:i+kmer_len]].append(i)
        for yeo_cluster in dict_Yeo_cluster_index:
            sp = yeo_cluster.split('.')
            motif = sp[2]
            region = sp[3] + '_300'
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
    fetch_seq_alu_code(args)
    build_feature(args.input, args.event_type, args.output, args.mutate)
    iDARTStmp_dir = os.path.dirname(os.path.abspath(args.input)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    logger.info('Deleting tmp files')
    remove_file(iDARTStmp_dir)