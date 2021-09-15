# -*- coding: UTF-8 -*-

"""
iDARTS - parse_vcf
"""

import os
import sys
from . import config
from cyvcf2 import VCF
from tqdm import tqdm

import logging
logger = logging.getLogger('iDARTS.get_resources')

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

def parse_vcf_SE(args): # parse vcf files for exon skipping events
    input_fn = args.input
    output_fn = args.output
    vcf_path = args.vcf_path
    vcf = VCF(vcf_path, strict_gt=True)

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
        exon_interval = (int(sp[EXON_START_INDEX]) - 300, int(sp[EXON_END_INDEX]) + 300)
        upstream_interval = (int(sp[UPSTREAM_ES_INDEX]), int(sp[UPSTREAM_EE_INDEX]) + 300)
        downstream_interval = (int(sp[DOWNSTREAM_ES_INDEX]) - 300, int(sp[DOWNSTREAM_EE_INDEX]))
        pos_list = []
        ref_list = []
        alt_list = []
        rs_list = []
        for interval in [exon_interval, upstream_interval, downstream_interval]:
            for record in vcf('{}:{}-{}'.format(chrom[3:], interval[0] - 1, interval[1] + 1)):
                if record.is_snp:
                    _pos = record.POS - 1
                    _ref = record.REF
                    _alt = str(record.ALT[0])
                    _rs = record.ID
                    if not fall_in_interval(_pos, interval):
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
    input_fn = args.input
    output_fn = args.output
    vcf_path = args.vcf_path
    event_type = args.event_type
    vcf = VCF(vcf_path, strict_gt=True)

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
        strand = sp[STRAND_INDEX]
        if event_type == 'A3SS' and strand == '+':
            long_exon_interval = (int(sp[LONGEXONSTART]) - 300, int(sp[LONGEXONEND]))
            flanking_exon_interval = (int(sp[FLANKINGEXONSTART]), int(sp[FLANKINGEXONEND]) + 300)
        if event_type == 'A3SS' and strand == '-':
            long_exon_interval = (int(sp[LONGEXONSTART]), int(sp[LONGEXONEND]) + 300)
            flanking_exon_interval = (int(sp[FLANKINGEXONSTART]) - 300, int(sp[FLANKINGEXONEND]))
        if event_type == 'A5SS' and strand == '+':
            long_exon_interval = (int(sp[LONGEXONSTART]), int(sp[LONGEXONEND]) + 300)
            flanking_exon_interval = (int(sp[FLANKINGEXONSTART]) - 300, int(sp[FLANKINGEXONEND]))
        if event_type == 'A5SS' and strand == '-':
            long_exon_interval = (int(sp[LONGEXONSTART]) - 300, int(sp[LONGEXONEND]))
            flanking_exon_interval = (int(sp[FLANKINGEXONSTART]), int(sp[FLANKINGEXONEND]) + 300)
        pos_list = []
        ref_list = []
        alt_list = []
        rs_list = []
        for interval in [long_exon_interval, flanking_exon_interval]:
            for record in vcf('{}:{}-{}'.format(chrom[3:], interval[0] - 1, interval[1] + 1)):
                if record.is_snp:
                    _pos = record.POS - 1
                    _ref = record.REF
                    _alt = str(record.ALT[0])
                    _rs = record.ID
                    if not fall_in_interval(_pos, interval):
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

def parse_vcf_RI(args):
    pass

def parser(args):
    event_type = args.event_type
    if event_type == 'SE':
        parse_vcf_SE(args)
    if event_type == 'RI':
        #parse_vcf_RI(args)
        logger.info('RI still under development')
        sys.exit(0)
    if event_type == 'A3SS' or event_type == 'A5SS':
        #parse_vcf_ASS(args)
        logger.info('A3SS and A5SS still under development')
        sys.exit(0)