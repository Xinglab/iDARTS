# -*- coding: UTF-8 -*-

"""
iDARTS - get_resources
Implements an internal downloading module for 
getting data from internet. Data includes training
data, cis feature files, etc.
This module depends on the url and md5sum stored in ``resources/download.yaml``
"""
from subprocess import *
from collections import defaultdict
from pkg_resources import resource_filename
import time
import logging
logger = logging.getLogger('iDARTS.build_feature')

from feature_utils import *

CHROM_INDEX = 3  # start from 0 in input file
STRAND_INDEX = 4
UPSTREAM_ES_INDEX = 7
UPSTREAM_EE_INDEX = 8
DOWNSTREAM_ES_INDEX = 9
DOWNSTREAM_EE_INDEX = 10
POS0BASE_INDEX = 12
REF_INDEX = 13
ALT_INDEX = 14

def get_region_coordinate(rmats_input_file_line):
    sp = rmats_input_file_line.strip().split('\t')
    _1 = [int(sp[UPSTREAM_ES_INDEX]), int(sp[UPSTREAM_EE_INDEX])]
    _2 = [int(sp[UPSTREAM_EE_INDEX]), int(sp[UPSTREAM_EE_INDEX]) + 300]
    _6 = [int(sp[DOWNSTREAM_ES_INDEX]) - 300, int(sp[DOWNSTREAM_ES_INDEX])]
    _7 = [int(sp[DOWNSTREAM_ES_INDEX]), int(sp[DOWNSTREAM_EE_INDEX])]
    if sp[STRAND_INDEX] == '+':
        C1, C2 = _1, _7
        I1_5p, I2_3p = _2, _6
    if sp[STRAND_INDEX] == '-':
        C1, C2 = _7, _1
        I1_5p, I2_3p = _6, _2
    return C1, I1_5p, I2_3p, C2

def fetch_seq_alu_code(args): # the input file (rMATS format)
    input_file = args.input
    event_type = args.event_type
    iDARTStmp_dir = os.path.dirname(os.path.abspath(input_file)) + '/_iDARTStmp/'
    makedirs(iDARTStmp_dir)
    fw = open(iDARTStmp_dir + 'event.{}.seq'.format(event_type), 'w')
    fw.write('ID\tC1\tI1_5p\tI2_3p\tC2\n')
    file_num_lines = get_file_num_lines(input_file)
    fp = open(input_file)
    fp.readline()
    logger.info('fetching sequences and Alu for splicing events')
    for line in tqdm(fp, total=file_num_lines - 1):
        sp = line.strip().split('\t')
        chrom = sp[CHROM_INDEX]
        strand = sp[STRAND_INDEX]
        C1, I1_5p, I2_3p, C2 = get_region_coordinate(line)
        if args.mutate == 'True':
            mutate_pos_list = sp[POS0BASE_INDEX].split(',')
            ref_base_list = sp[REF_INDEX].split(',')
            alt_base_list = sp[ALT_INDEX].split(',')
            C1_seq = seq_fetch(chrom, C1[0], C1[1], strand, mutate_pos_list, alt_base_list)
            I1_5p_seq = seq_fetch(chrom, I1_5p[0], I1_5p[1], strand, mutate_pos_list, alt_base_list)
            I2_3p_seq = seq_fetch(chrom, I2_3p[0], I2_3p[1], strand, mutate_pos_list, alt_base_list)
            C2_seq = seq_fetch(chrom, C2[0], C2[1], strand, mutate_pos_list, alt_base_list)
        else:
            C1_seq = seq_fetch(chrom, C1[0], C1[1], strand)
            I1_5p_seq = seq_fetch(chrom, I1_5p[0], I1_5p[1], strand)
            I2_3p_seq = seq_fetch(chrom, I2_3p[0], I2_3p[1], strand)
            C2_seq = seq_fetch(chrom, C2[0], C2[1], strand)
        fw.write(sp[0] + '\t' + '\t'.join([C1_seq, I1_5p_seq, I2_3p_seq, C2_seq]) + '\n')
    fp.close()
    fw.close()

def get_regional_conservation_score(annotation_list, dict_chrom_phatsConsScore):
    strand = annotation_list[STRAND_INDEX]
    chrom = annotation_list[CHROM_INDEX]
    if strand == '+':
        I1_5p_index = int(annotation_list[UPSTREAM_EE_INDEX])
        I2_3p_index = int(annotation_list[DOWNSTREAM_ES_INDEX])
        I1_5p_conservation = [x if x == x else 0 for x in dict_chrom_phatsConsScore[chrom].values(chrom, I1_5p_index, I1_5p_index + 300)]
        I2_3p_conservation = [x if x == x else 0 for x in dict_chrom_phatsConsScore[chrom].values(chrom, I2_3p_index - 300, I2_3p_index)]
    if strand == '-':
        I2_3p_index = int(annotation_list[UPSTREAM_EE_INDEX])
        I1_5p_index = int(annotation_list[DOWNSTREAM_ES_INDEX])
        I1_5p_conservation = [x if x == x else 0 for x in dict_chrom_phatsConsScore[chrom].values(chrom, I1_5p_index - 300, I1_5p_index)][::-1]
        I2_3p_conservation = [x if x == x else 0 for x in dict_chrom_phatsConsScore[chrom].values(chrom, I2_3p_index, I2_3p_index + 300)][::-1]
    return I1_5p_conservation, I2_3p_conservation

def tmp():
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

    if len(dict_seq['C1']) >= 50:
        seq = dict_seq['C1'][-50:] + dict_seq['I1_5p'][0:250]
        junction_start = 50
    else:
        seq = dict_seq['C1'] + dict_seq['I1_5p'][0:250]
        junction_start = len(dict_seq['C1'])
    unpair_prob_list = rnafold_predict_pu_score(seq, tmp_file + '.1', event_type)
    for i in xrange(3):
        start, end = 70 * i, 70 * (i + 1)
        dict_rnafold_pu_scores['SecStr.max{}_{}.I1_5p'.format(start + 1, end)] = np.max(unpair_prob_list[start + junction_start: end + junction_start])
        dict_rnafold_pu_scores['SecStr.avg{}_{}.I1_5p'.format(start + 1, end)] = np.mean(unpair_prob_list[start + junction_start: end + junction_start])
    dict_rnafold_pu_scores['SecStr.maxJunc.I1_5p'] = np.max(unpair_prob_list[junction_start:junction_start + 10])
    dict_rnafold_pu_scores['SecStr.avgJunc.I1_5p'] = np.mean(unpair_prob_list[junction_start:junction_start + 10])

    seq = dict_seq['I2_3p'][-250:] + dict_seq['C2'][0:min(50, len(dict_seq['C2']))]
    for i in xrange(3):
        start, end = 70 * i, 70 * (i + 1)
        dict_rnafold_pu_scores['SecStr.max{}_{}.I2_3p'.format(start + 1, end)] = np.max(unpair_prob_list[250 - end:250 - start])
        dict_rnafold_pu_scores['SecStr.avg{}_{}.I2_3p'.format(start + 1, end)] = np.mean(unpair_prob_list[250 - end:250 - start])
    ## junction region for 3' splice site longer (20nt)
    dict_rnafold_pu_scores['SecStr.maxJunc.I2_3p'] = np.max(unpair_prob_list[230:250])
    dict_rnafold_pu_scores['SecStr.avgJunc.I2_3p'] = np.mean(unpair_prob_list[230:250])
    return dict_rnafold_pu_scores

def build_feature(input_file, event_type, output_file):
    logger.info('building features for {}'.format(event_type))
    iDARTStmp_dir = os.path.dirname(os.path.abspath(input_file)) + '/_iDARTStmp/'
    splice_feature_fn = resource_filename('iDARTS.resources.features', 'SpliceCode_feature_list_{}.txt'.format(event_type))
    dict_chrom_phatsConsScore = load_phatsConsScore()
    dict_ESE_matrix = load_ESE_matrix_PSSM()
    dict_rna_binding_pssm = load_rna_binding_motif_pssm()
    dict_rna_binding_motif = load_rna_binding_motif()
    splice_feature_list = load_splice_feature_list(splice_feature_fn)
    donor_maxent_matrix = maxent_fast.load_matrix(5)
    acceptor_maxent_matrix = maxent_fast.load_matrix(3)
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
    seq_fp = open(iDARTStmp_dir + 'event.{}.seq'.format(event_type))
    seq_fp_header = seq_fp.readline().strip()
    fw = open(output_file, 'w')
    fw.write(header + '\t' + '\t'.join(splice_feature_list) + '\n')

    for annotation in tqdm(annotation_fp, total = file_num_lines - 1):
        seq_file_line_list = seq_fp.readline().strip().split('\t')
        if len(annotation.strip()) == 0:
            continue
        annotation_list = annotation.strip().split('\t')
        C1_seq, I1_5p_seq, I2_3p_seq, C2_seq =  seq_file_line_list[1:]
        C1_seq_matrix = np.array([NT_DICT[x] for x in C1_seq])
        C2_seq_matrix = np.array([NT_DICT[x] for x in C2_seq])
        I1_5p_seq_matrix = np.array([NT_DICT[x] for x in I1_5p_seq])
        I2_3p_seq_matrix = np.array([NT_DICT[x] for x in I2_3p_seq])
        dict_seq2matrix = {}
        dict_seq2matrix['C1'] = C1_seq_matrix
        dict_seq2matrix['C2'] = C2_seq_matrix
        dict_seq2matrix['I1_5p'] = I1_5p_seq_matrix
        dict_seq2matrix['I2_3p'] = I2_3p_seq_matrix
        dict_seq = {}
        dict_seq['C1'] = C1_seq
        dict_seq['C2'] = C2_seq
        dict_seq['I1_5p'] = I1_5p_seq
        dict_seq['I2_3p'] = I2_3p_seq
        ## get exon length features
        dict_seq_feature2value = {}
        length_I = int(annotation_list[DOWNSTREAM_ES_INDEX]) - int(annotation_list[UPSTREAM_EE_INDEX])
        if annotation_list[STRAND_INDEX] == '+':
            length_C1 = int(annotation_list[UPSTREAM_EE_INDEX]) - int(annotation_list[UPSTREAM_ES_INDEX])
            length_C2 = int(annotation_list[DOWNSTREAM_EE_INDEX]) - int(annotation_list[DOWNSTREAM_ES_INDEX])
        else:
            length_C1 = int(annotation_list[DOWNSTREAM_EE_INDEX]) - int(annotation_list[DOWNSTREAM_ES_INDEX])
            length_C2 = int(annotation_list[UPSTREAM_EE_INDEX]) - int(annotation_list[UPSTREAM_ES_INDEX])
        dict_seq_len = {}
        dict_seq_len['C1'] = length_C1
        dict_seq_len['C2'] = length_C2
        dict_seq_feature2value['LogLen.C1'] = np.log10(length_C1)
        dict_seq_feature2value['LogLen.C2'] = np.log10(length_C2)
        dict_seq_feature2value['LogLen.I'] = np.log10(length_I)
        dict_seq_feature2value['LogLenRatio.C1_I'] = np.log10((length_C1 + 0.0) / length_I)
        dict_seq_feature2value['LogLenRatio.C2_I'] = np.log10((length_C2 + 0.0) / length_I)
        if length_I % 3 == 0:
            dict_seq_feature2value['FrameShift.I'] = 0
        else:
            dict_seq_feature2value['FrameShift.I'] = 1
        ### Conservation score
        I1_5p_conservation, I2_3p_conservation = get_regional_conservation_score(annotation_list, dict_chrom_phatsConsScore)
        dict_region2conservation = {}
        dict_region2conservation['I1_5p'] = I1_5p_conservation
        dict_region2conservation['I2_3p'] = I2_3p_conservation
        dict_seq_feature2value['Cons.MeanP1_100.I1_5p'] = np.mean(I1_5p_conservation[0:100])
        dict_seq_feature2value['Cons.MeanP1_100.I2_3p'] = np.mean(I2_3p_conservation[-100:])
        JunctionAvgP1_100_I1_5p = np.mean(I1_5p_conservation[0:2]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.I1_5p']) + 1e-4)
        JunctionAvgP1_100_I2_3p = np.mean(I2_3p_conservation[-2:]) / (np.mean(dict_seq_feature2value['Cons.MeanP1_100.I2_3p']) + 1e-4)
        JunctionAvgP1_100_I1_5p = np.sqrt(JunctionAvgP1_100_I1_5p)
        JunctionAvgP1_100_I2_3p = np.sqrt(JunctionAvgP1_100_I2_3p)
        ## get splice AG/GT position AltAGpos AltGTpos
        if 'AG' in I2_3p_seq:
            dict_seq_feature2value['AltAGpos'] = np.log2(300 - I2_3p_seq.rindex('AG')) - 1
        else:
            dict_seq_feature2value['AltAGpos'] = np.log2(300) - 1
        if 'GT' in I1_5p_seq:
            dict_seq_feature2value['AltGTpos'] = np.log2(I1_5p_seq.index('GT') + 2) - 1
        else:
            dict_seq_feature2value['AltGTpos'] = np.log2(300) - 1
        ### get splice site strength by maxent
        # maxent_fast.score5('cagGTAAGT', matrix=donor)
        # maxent_fast.score3('ttccaaacgaacttttgtAGgga', matrix=acceptor)
        donor_seq = C1_seq[-3:] + I1_5p_seq[0:6]
        acceptor_seq = I2_3p_seq[-20:] + C2_seq[0:3]
        if I1_5p_seq[0:2] == 'GT' or I1_5p_seq[0:2] == 'GC' or I1_5p_seq[0:2] == 'AT':
            donor_score = round(maxent_fast.score5(donor_seq, matrix = donor_maxent_matrix), 2)
        else:
            donor_score = -49.17
        if I2_3p_seq[-2:] == 'AG' or I2_3p_seq[-2:] == 'AC':
            acceptor_score = round(maxent_fast.score3(acceptor_seq, matrix = acceptor_maxent_matrix), 2)
        else:
            acceptor_score = -67.35
        dict_seq_feature2value['donor.maxent'] = donor_score
        dict_seq_feature2value['acceptor.maxent'] = acceptor_score
        acceptor_prob = 2 ** acceptor_score / (2 ** acceptor_score + 1.0)
        donor_prob = 2 ** donor_score / (2 ** donor_score + 1.0)
        ## get acceptor and donor probablity
        dict_seq_feature2value['acceptor.maxent.weighted.JunctionAvgP1_100.I2_3p'] = acceptor_prob * JunctionAvgP1_100_I2_3p
        dict_seq_feature2value['donor.maxent.weighted.JunctionAvgP1_100.I1_5p'] = donor_prob * JunctionAvgP1_100_I1_5p
        # ESS ESE Matrix score
        ESS_ESE_RBP_LIST = sorted(dict_ESE_matrix.keys())
        for region in ['C1', 'C2']:
            seq_matrix = dict_seq2matrix[region]
            score_list = get_ESE_matrix_score(seq_matrix, dict_ESE_matrix)
            for score, rbp in zip(score_list, ESS_ESE_RBP_LIST):
                dict_seq_feature2value['PSSM.{}.{}'.format(rbp, region)] = score

        # Rosenberg and ESR score
        dict_seq_feature2value['A5SS_R1.C1'] = get_rosenberg_mean_score(C1_seq, rosenberg_matrix, ['A5SS_R1'])['A5SS_R1']
        dict_seq_feature2value['A3SS_R2.C2'] = get_rosenberg_mean_score(C2_seq, rosenberg_matrix, ['A3SS_R2'])['A3SS_R2']
        dict_seq_feature2value['A3SS_R1.I2_3p'] = get_rosenberg_mean_score(I2_3p_seq, rosenberg_matrix, ['A3SS_R1'])['A3SS_R1']
        dict_seq_feature2value['A5SS_R2.I1_5p'] = get_rosenberg_mean_score(I1_5p_seq, rosenberg_matrix, ['A5SS_R2'])['A5SS_R2']
        dict_seq_feature2value['ESR.C1'] = get_esr_mean_score(C1_seq, esr_matrix)
        dict_seq_feature2value['ESR.C2'] = get_esr_mean_score(C2_seq, esr_matrix)

        # ISS and ISE counts
        ISS_ISE_Kmers = [x.split('.')[-1] for x in dict_ISE_ISS_index]
        intron_upstream_seq = I2_3p_seq[-100:-15]
        intron_downstream_seq = I1_5p_seq[15:100]
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
        # RBP binding motif PSSM profile annotation
        for region in ['C1', 'C2', 'I1_5p', 'I2_3p']:
            seq_matrix = dict_seq2matrix[region]
            score_list = get_pssm_score(seq_matrix, dict_rna_binding_pssm)
            for index, binding_protein in enumerate(sorted(dict_rna_binding_pssm.keys())):
                dict_seq_feature2value[binding_protein + '.' + region] = score_list[index]
        ## ShortSeq Annotation
        for region in ['C1', 'C2', 'I1_5p', 'I2_3p']:
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
        for region in ['C1', 'C2']:
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
        for region in ['I1_5p', 'I2_3p']:
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
            for region in ['I2_3p', 'I1_5p']:
                dict_seq2kmer_loc[kmer_len][region] = defaultdict(lambda:[])
                seq = dict_seq[region]
                for i in xrange(0, 300 - kmer_len + 1):
                    dict_seq2kmer_loc[kmer_len][region][seq[i:i+kmer_len]].append(i)
        for yeo_cluster in dict_Yeo_cluster_index:
            sp = yeo_cluster.split('.')
            motif = sp[2]
            region = sp[-1]
            kmer_len = len(motif)
            region_conservation_scores = dict_region2conservation[region]
            pos_list = dict_seq2kmer_loc[kmer_len][region][motif]
            conserve_score = 0.0
            if len(pos_list) > 0:
                for pos in pos_list:
                    conserve_score += np.mean(region_conservation_scores[pos:pos+kmer_len])
            yeo_name = '.'.join(sp[0:3]) + '.' + region
            dict_seq_feature2value[yeo_name] = conserve_score
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
    build_feature(args.input, args.event_type, args.output)
    iDARTStmp_dir = os.path.dirname(os.path.abspath(args.input)) + '/_iDARTStmp/'
    logger.info('Deleting tmp files')
    remove_file(iDARTStmp_dir)
