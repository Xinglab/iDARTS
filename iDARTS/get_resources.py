# -*- coding: UTF-8 -*-

"""
iDARTS - get_resources
Implements an internal downloading module for 
getting data from internet. Data includes training
data, cis feature files, etc.
This module depends on the url and md5sum stored in ``resources/download.yaml``
"""

import os
from subprocess import *
from . import config

import logging
logger = logging.getLogger('iDARTS.get_resources')

def makedirs(_dir):
    try:
        os.stat(_dir)
    except:
        os.makedirs(_dir)

def remove(fn):
    cmd = 'rm -rf ' + fn
    call(cmd, shell = True)

def get_md5(fn):
    cmd = 'md5sum ' + fn
    process = Popen(cmd, shell = True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    md5value = stdout.strip().split()[0]
    return md5value

def download(url, md5value):
    try_count = 0
    download_flag = False
    while try_count < 10:
        download_cmd = 'wget -c ' + url
        call(download_cmd, shell = True)
        fn = url.split('/')[-1]
        try_count += 1
        md5check = get_md5(fn)
        if md5check == md5value:
            download_flag = True
            return download_flag
        else:
            remove(fn)

def wigFix2bigwig(chrom):
    HG19_GENOME = config.HG19_GENOME
    cmd = 'wigToBigWig {0}.phastCons46way.wigFix.gz {1} {0}.phastCons46way.bw'.format(chrom, HG19_GENOME)
    print(cmd)
    call(cmd, shell = True)

def gzip_decompress(fn):
    cmd = 'gzip -d ' + fn
    print(cmd)
    call(cmd, shell = True)

def download_resources(outdir):
    cwd_dir = os.getcwd()
    os.chdir(outdir)
    data_config = config.DOWNLOAD_CONFIG
    for d in data_config:
        url_list = data_config[d]['url']
        md5_list = data_config[d]['md5sum']
        url_list = url_list if isinstance(url_list, list) else [url_list]
        md5_list = md5_list if isinstance(md5_list, list) else [md5_list]
        for url, md5 in zip(url_list, md5_list):
            fn = url.split('/')[-1]
            download_flag = download(url, md5)
            if download_flag:
                if d == 'phastCons46way':
                    wigFix2bigwig(fn.split('.')[0])
                    remove(fn)
                if d == 'hg19_fasta':
                    gzip_decompress(fn)
    os.chdir(cwd_dir)

def check_resources_md5(outdir):
    cwd_dir = os.getcwd()
    os.chdir(outdir)
    logger.info('Checking downloaded data md5value')
    CHECK_MD5_CONFIG = config.CHECK_MD5_CONFIG
    data_config = config.DOWNLOAD_CONFIG
    data2md5_dict = {}
    for d in data_config:
        url_list = data_config[d]['url']
        url_list = url_list if isinstance(url_list, list) else [url_list]
        for url in url_list:
            fn = url.split('/')[-1]
            if d == 'hg19_fasta':
                out_fn = 'hg19.fa'
                if not os.path.exists(out_fn):
                    return False
                data2md5_dict[out_fn] = get_md5(out_fn)
            if d == 'phastCons46way':
                out_fn = '.'.join(fn.split('.')[0:2]) + '.bw'
                if not os.path.exists(out_fn):
                    return False
                data2md5_dict[out_fn] = get_md5(out_fn)
    resources_download_flag = True
    for data in data2md5_dict:
        if data2md5_dict[data] != CHECK_MD5_CONFIG[data]:
            return False
    os.chdir(cwd_dir)
    return resources_download_flag

def parser(args):
    outdir = args.out_dir
    iDARTS_DIRECTORY = config.iDARTS_DIRECTORY
    makedirs(iDARTS_DIRECTORY)
    if outdir == None:
        logger.info('No path provided, resource will be downloaded in ~/.iDARTS/resources/')
        outdir = os.path.join(iDARTS_DIRECTORY, 'resources')
    else:
        outdir = os.path.abspath(outdir)
    makedirs(outdir)
    resources_download_flag = check_resources_md5(outdir)
    if not resources_download_flag:
        download_resources(outdir)
        resources_download_flag = check_resources_md5(outdir)
    if resources_download_flag:
        logger.info('resource download successfully')
        resource_config = iDARTS_DIRECTORY + 'resource.yaml'
        fw = open(resource_config, 'w')
        fw.write('data_download: True\n')
        fw.write('data_dir: ' + outdir + '\n')
        fw.close()
    else:
        logger.info('resource download unsuccessfully')
    return