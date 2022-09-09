# -*- coding: UTF-8 -*-

"""
iDARTS - parse_vcf
implements a function to process input vcf file 
for variants matched splicing events
"""

import os
import sys
from cyvcf2 import VCF
from tqdm import tqdm
import argparse as ap
from pkg_resources import resource_filename
import tempfile
import datetime
from . import feature_utils
from . import build_feature_SE
from . import build_feature_ASS
from . import iDARTS_pred

import logging
logger = logging.getLogger('iDARTS.parse_vcf')

# chromosomes index
CHR_INDEX_DICT = {'1':1, '2':2, '3':3, '4':4, '5':5, '6':6, '7':7, '8':8, '9':9, '10':10, '11':11, '12':12, 
  '13':13, '14':14, '15':15, '16':16, '17':17, '18':18, '19':19, '20':20, '21':21, '22':22, 'X':23, 'Y':24}

SE_exon_triplets = resource_filename('iDARTS.resources.AS_events', 'gencode.v26lift37.SE.exon_triplets.txt')
A5SS_events = resource_filename('iDARTS.resources.AS_events', 'gencode.v26lift37.A5SS.txt')
A3SS_events = resource_filename('iDARTS.resources.AS_events', 'gencode.v26lift37.A3SS.txt')
TEMP_FILE_NAME = next(tempfile._get_candidate_names())

def fall_in_interval(pos, interval):
    if pos >= interval[0] and pos < interval[1]:
        return True
    else:
        return False

def get_file_num_lines(filename):
    fp = open(filename)
    count = 0
    for line in fp:
        count += 1
    fp.close()
    return count

def sort_snp(snp1, snp2):
    sp = snp1.split('_')
    sq = snp2.split('_')
    if sp[0] == sq[0]:
        if int(sp[1]) >= int(sq[1]):
            return 1
        return -1
    elif sp[0] in CHR_INDEX_DICT and sq[0] in CHR_INDEX_DICT:
        if CHR_INDEX_DICT[sp[0]] >= CHR_INDEX_DICT[sq[0]]:
            return 1
        else:
            return -1
    elif sp[0] in CHR_INDEX_DICT:
        return 1
    elif sq[0] in CHR_INDEX_DICT:
        return -1
    else:
        if sp[0] >= sq[0]:
            return 1
        else:
            return -1

def check_vcf_file(vcf_path):
    if os.path.exists(vcf_path + '.tbi'):
        return True
    sys.stderr.write("\nCould not retrieve index file for '{}'".format(vcf_path))
    sys.stderr.write("\n- VCF needs to be sorted and index first. \n- Example code: tabix -p vcf example.vcf.gz\n")
    sys.exit(0)

def parse_vcf_SE(args): # parse vcf files for exon skipping events
    if args.input == None:
        input_fn = SE_exon_triplets
    else:
        input_fn = args.input

    out_dir = os.path.dirname(os.path.abspath(args.outFileNamePrefix)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    feature_utils.makedirs(out_dir)
    output_fn = out_dir + '{}.vcf.SE.event'.format(TEMP_FILE_NAME)

    vcf_path = args.vcf_path
    vcf = VCF(vcf_path, strict_gt = True)

    CHROM_INDEX = 3  # start from 0 in input file, rMATS SE format
    STRAND_INDEX = 4
    UPSTREAM_ES_INDEX = 7
    UPSTREAM_EE_INDEX = 8
    EXON_START_INDEX = 5
    EXON_END_INDEX = 6
    DOWNSTREAM_ES_INDEX = 9
    DOWNSTREAM_EE_INDEX = 10

    file_num_lines = get_file_num_lines(input_fn)
    fp = open(input_fn)
    fw = open(output_fn, 'w')
    header = fp.readline().strip()
    fw.write(header + '\t' + 'SNP_ID\tPos0base\tRef\tAlt\n')
    for line in tqdm(fp, total = file_num_lines - 1):
        sp = line.strip().split('\t')
        chrom = sp[CHROM_INDEX]
        if chrom[3:] not in vcf.seqnames:
            continue
        exon_interval = (int(sp[EXON_START_INDEX]) - 300, int(sp[EXON_END_INDEX]) + 300)
        pos_list = []
        ref_list = []
        alt_list = []
        rs_list = []
        for record in vcf('{}:{}-{}'.format(chrom[3:], exon_interval[0] - 1, exon_interval[1] + 1)):
            if record.is_snp:
                _pos = record.POS - 1
                _ref = record.REF
                _alt = str(record.ALT[0])
                _rs = record.ID
                if not fall_in_interval(_pos, exon_interval):
                    continue
                pos_list.append(_pos)
                ref_list.append(_ref)
                alt_list.append(_alt)
                if _rs:
                    rs_list.append(_rs)
                else:
                    rs_list.append('NA')
        output_line = '\t'.join(sp[1:])
        if len(pos_list) == 0:
            continue
        out_line_list = []
        for rs, pos, ref, alt in zip(rs_list, pos_list, ref_list, alt_list):
            rsID = '{}_{}_{}_{}_b37'.format(chrom.replace('chr', ''), pos + 1, ref, alt)
            ID = '|'.join(sp[3:11]) + '@' + rsID
            annotate_snp_info = '{}\t{}\t{}\t{}'.format(rs, pos, ref, alt)
            out_line_list.append(ID + '\t' + output_line + '\t' + annotate_snp_info)
        out_line_list = sorted(list(set(out_line_list)))
        fw.write('\n'.join(out_line_list) + '\n')
    fp.close()
    fw.close()

def parse_vcf_ASS(args): # parse vcf files for alternative 3' or 5' splicing events
    if args.input == None:
        if args.event_type == 'A3SS':
            input_fn = A3SS_events
        else:
            input_fn = A5SS_events
    else:
        input_fn = args.input

    out_dir = os.path.dirname(os.path.abspath(args.outFileNamePrefix)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    feature_utils.makedirs(out_dir)
    output_fn = out_dir + '{}.vcf.{}.event'.format(TEMP_FILE_NAME, args.event_type)

    vcf_path = args.vcf_path
    event_type = args.event_type
    vcf = VCF(vcf_path, strict_gt = True)

    CHROM_INDEX = 3  # start from 0 in input file
    STRAND_INDEX = 4
    LONGEXONSTART = 5
    LONGEXONEND = 6
    SHORTEXONSTART = 7
    SHORTEXONEND = 8
    FLANKINGEXONSTART = 9
    FLANKINGEXONEND = 10

    file_num_lines = get_file_num_lines(input_fn)
    fp = open(input_fn)
    fw = open(output_fn, 'w')
    header = fp.readline().strip()
    fw.write(header + '\t' + 'SNP_ID\tPos0base\tRef\tAlt\n')
    for line in tqdm(fp, total = file_num_lines - 1):
        sp = line.strip().split('\t')
        chrom = sp[CHROM_INDEX]
        if chrom[3:] not in vcf.seqnames:
            continue
        strand = sp[STRAND_INDEX]
        if event_type == 'A3SS' and strand == '+':
            long_exon_interval = (int(sp[LONGEXONSTART]) - 300, int(sp[LONGEXONEND]))
        if event_type == 'A3SS' and strand == '-':
            long_exon_interval = (int(sp[LONGEXONSTART]), int(sp[LONGEXONEND]) + 300)
        if event_type == 'A5SS' and strand == '+':
            long_exon_interval = (int(sp[LONGEXONSTART]), int(sp[LONGEXONEND]) + 300)
        if event_type == 'A5SS' and strand == '-':
            long_exon_interval = (int(sp[LONGEXONSTART]) - 300, int(sp[LONGEXONEND]))
        pos_list = []
        ref_list = []
        alt_list = []
        rs_list = []
        for record in vcf('{}:{}-{}'.format(chrom[3:], long_exon_interval[0] - 1, long_exon_interval[1] + 1)):
            if record.is_snp:
                _pos = record.POS - 1
                _ref = record.REF
                _alt = str(record.ALT[0])
                _rs = record.ID
                if not fall_in_interval(_pos, long_exon_interval):
                    continue
                pos_list.append(_pos)
                ref_list.append(_ref)
                alt_list.append(_alt)
                if _rs:
                    rs_list.append(_rs)
                else:
                    rs_list.append('NA')
        output_line = '\t'.join(sp[1:])
        if len(pos_list) == 0:
            continue
        out_line_list = []
        for rs, pos, ref, alt in zip(rs_list, pos_list, ref_list, alt_list):
            rsID = '{}_{}_{}_{}_b37'.format(chrom.replace('chr', ''), pos + 1, ref, alt)
            ID = '|'.join(sp[3:11]) + '@' + rsID
            annotate_snp_info = '{}\t{}\t{}\t{}'.format(rs, pos, ref, alt)
            out_line_list.append(ID + '\t' + output_line + '\t' + annotate_snp_info)
        out_line_list = sorted(list(set(out_line_list)))
        fw.write('\n'.join(out_line_list) + '\n')
    fp.close()
    fw.close()

def get_featureNpred(args):
    out_dir = os.path.dirname(os.path.abspath(args.outFileNamePrefix)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    input_vcf_event = out_dir + '{}.vcf.{}.event'.format(TEMP_FILE_NAME, args.event_type)

    output_event_ref_feature = out_dir + '{}.vcf.{}.event.ref.feature'.format(TEMP_FILE_NAME, args.event_type)
    output_event_mutate_feature = out_dir + '{}.vcf.{}.event.mutate.feature'.format(TEMP_FILE_NAME, args.event_type)

    output_event_ref_pred = out_dir + '{}.vcf.{}.event.ref.pred'.format(TEMP_FILE_NAME, args.event_type) 
    output_event_mutate_pred = out_dir + '{}.vcf.{}.event.mutate.pred'.format(TEMP_FILE_NAME, args.event_type) 

    build_feature_args = ap.Namespace()
    build_feature_args.event_type = args.event_type
    build_feature_args.input = input_vcf_event

    pred_args = ap.Namespace()
    pred_args.event_type = args.event_type
    pred_args.expr = args.expr

    logger.info('Build features for events with reference alleles')
    build_feature_args.output = output_event_ref_feature
    build_feature_args.mutate = 'False'
    if args.event_type == 'SE':
        build_feature_SE.main(build_feature_args)
    else:
        build_feature_ASS.main(build_feature_args)

    logger.info('Build features for events with mutate alleles')
    build_feature_args.output = output_event_mutate_feature
    build_feature_args.mutate = 'True'
    if args.event_type == 'SE':
        build_feature_SE.main(build_feature_args)
    else:
        build_feature_ASS.main(build_feature_args)

    logger.info('Predict PSI values for events with reference alleles')
    pred_args.input = output_event_ref_feature
    pred_args.output = output_event_ref_pred
    iDARTS_pred.make_prediction(pred_args)

    logger.info('Predict PSI values for events with alternative alleles')
    pred_args.input = output_event_mutate_feature
    pred_args.output = output_event_mutate_pred
    iDARTS_pred.make_prediction(pred_args)

def get_vcf_with_iDARTS_pred(args):
    logger.info('Predict deltaPSI as the effects of variants on splicing')
    out_dir = os.path.dirname(os.path.abspath(args.outFileNamePrefix)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME)
    output_event_mutate_feature = out_dir + '{}.vcf.{}.event.mutate.feature'.format(TEMP_FILE_NAME, args.event_type)
    output_event_ref_pred = out_dir + '{}.vcf.{}.event.ref.pred'.format(TEMP_FILE_NAME, args.event_type) 
    output_event_mutate_pred = out_dir + '{}.vcf.{}.event.mutate.pred'.format(TEMP_FILE_NAME, args.event_type)

    dict_snp2delta_psi = {}
    dict_snp_id2rsID = {}

    fp = open(output_event_mutate_pred)
    mutate_feature_fp = open(output_event_mutate_feature)
    ref_pred_fp = open(output_event_ref_pred)

    tissue_list = fp.readline().strip().split('\t')[1:]
    event_header = mutate_feature_fp.readline().strip().split('\t')
    ref_pred_fp.readline()

    fw = open(os.path.abspath(args.outFileNamePrefix) + '.iDARTS_prediction', 'w')
    header = ['SNP_ID', 'max_deltaPSI','Pos0base', 'Ref', 'Alt'] + event_header[1:event_header.index('SNP_ID')] + \
             [tissue + '_' + 'ref_PSI' for tissue in tissue_list] + \
             [tissue + '_' + 'deltaPSI' for tissue in tissue_list]
    fw.write('\t'.join(header) + '\n')
    for line in fp:
        sp = line.strip().split('\t')
        sq = mutate_feature_fp.readline().strip().split('\t')
        ref_pred_sp = ref_pred_fp.readline().strip().split('\t')
        chrom = sq[event_header.index('chr')]
        event = '|'.join(sq[event_header.index('chr'):event_header.index('SNP_ID')])
        delta_psi_list = [float(mutate) - float(ref) for mutate, ref in zip(sp[1:], ref_pred_sp[1:])]
        max_delta_psi = max(delta_psi_list, key = abs)
        rsID, pos, ref, alt = sq[event_header.index('SNP_ID')], sq[event_header.index('Pos0base')], \
                                sq[event_header.index('Ref')], sq[event_header.index('Alt')]
        snp_id = '{}_{}_{}_{}_b37'.format(chrom.replace('chr', ''), int(pos) + 1, ref, alt)
        fw.write('{}\t{}\t{}\t{}\t{}\t'.format(rsID, round(max_delta_psi, 4), pos, ref, alt))
        fw.write('\t'.join(sq[1:event_header.index('SNP_ID')]) + '\t')
        fw.write('\t'.join([str(round(float(ref), 4)) for ref in ref_pred_sp[1:]]) + '\t')
        fw.write('\t'.join([str(round(delta_psi, 4)) for delta_psi in delta_psi_list]) + '\n')
        dict_snp_id2rsID[snp_id] = rsID
        if snp_id in dict_snp2delta_psi:
            if abs(dict_snp2delta_psi[snp_id]) >= abs(max_delta_psi):
                continue
        dict_snp2delta_psi[snp_id] = max_delta_psi
    fp.close()
    fw.close()

    # sort SNP ID
    sorted_snp_id = sorted(dict_snp2delta_psi.keys(), cmp = sort_snp)
    fw = open(os.path.abspath(args.outFileNamePrefix) + '.vcf', 'w')
    fw.write('##fileformat=VCFv4.2\n')
    fw.write('##fileDate={}\n'.format(datetime.date.today().strftime('20%y%m%d')))
    fw.write('##reference=GRCh37/hg19\n')
    fw.write('##INFO=<ID=iDARTS,Number=.,Type=String,Description="iDARTS variant prediction. Maximum absolute deltaPSI across all tissues">\n')
    fw.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
    for snp_id in sorted_snp_id:
        delta_psi = dict_snp2delta_psi[snp_id]
        sp = snp_id.split('_')
        fw.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sp[0], sp[1], dict_snp_id2rsID[snp_id], sp[2], sp[3], '.', '.',
            'iDARTS=' + args.event_type + '|' + str(round(delta_psi, 4))))
    fw.close()

def parser(args):
    check_vcf_file(args.vcf_path)
    event_type = args.event_type
    if event_type == 'SE':
        parse_vcf_SE(args)
    if event_type == 'A3SS' or event_type == 'A5SS':
        parse_vcf_ASS(args)
    get_featureNpred(args)
    get_vcf_with_iDARTS_pred(args)
    ## remove tmp files
    logger.info('Deleting tmp files')
    feature_utils.remove_file(os.path.dirname(os.path.abspath(args.outFileNamePrefix)) + '/_iDARTS_{}_tmp/'.format(TEMP_FILE_NAME))