# -*- coding: UTF-8 -*-

"""
iDARTS - get_resources
Implements an internal downloading module for 
getting data from internet. Data includes trans-RBP expressions, 
cis feature files, etc.
This module depends on the url and md5sum stored in ``resources/download.yaml``
"""

import os
from subprocess import *
import yaml
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

def download_resources(outdir):
    cwd_dir = os.path.abspath(os.getcwd())
    os.chdir(outdir)
    data_config = config.DOWNLOAD_CONFIG
    url, md5 = data_config['resources']['url'], data_config['resources']['md5sum']
    logger.info('Downloading resources')
    download_flag = download(url, md5)
    if download_flag:
        logger.info('Unzipping resources')
        cmd = 'tar -xvzf resources.tar.gz'
        return_code = call(cmd, shell = True)
        if return_code == 0:
            remove('resources.tar.gz')
    os.chdir(cwd_dir)

def check_resources_md5(outdir):
    cwd_dir = os.path.abspath(os.getcwd())
    os.chdir(outdir)
    logger.info('Checking data md5value')
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
                    os.chdir(cwd_dir)
                    return False
                data2md5_dict[out_fn] = get_md5(out_fn)
            if d == 'phastCons46way':
                out_fn = '.'.join(fn.split('.')[0:2]) + '.bw'
                if not os.path.exists(out_fn):
                    os.chdir(cwd_dir)
                    return False
                data2md5_dict[out_fn] = get_md5(out_fn)
    for data in data2md5_dict:
        if data2md5_dict[data] != CHECK_MD5_CONFIG[data]:
            logger.info('check ' + data + ' md5 integrity check failed')
            os.chdir(cwd_dir)
            return False
        else:
            logger.info('check ' + data + ' md5 integrity check passed')
    os.chdir(cwd_dir)
    return True

def get_resource_download_status():
    if os.path.exists(config.iDARTS_DIRECTORY + 'resource.yaml'):
        DATA_CONFIG = yaml.safe_load(open(config.iDARTS_DIRECTORY + 'resource.yaml', 'r'))
        if DATA_CONFIG['data_download'] == True and os.path.exists(DATA_CONFIG['data_dir']):
            logger.info('resources already downloaded: ' + DATA_CONFIG['data_dir'])
            return True
    return False

def parser(args):
    if args.force == 'False' and get_resource_download_status():
        return True
    if args.out_dir == None:
        logger.info('No download directory specified, resources will be downloaded in ~/.iDARTS/')
        outdir = os.path.join(config.iDARTS_DIRECTORY, 'resources')
    else:
        outdir = os.path.join(os.path.abspath(args.out_dir), 'resources')
    makedirs(outdir)
    resources_download_flag = check_resources_md5(outdir)
    if not resources_download_flag:
        download_resources(args.out_dir)
        resources_download_flag = check_resources_md5(outdir)
    if resources_download_flag:
        logger.info('resources download successfully')
        resource_config = config.iDARTS_DIRECTORY + 'resource.yaml'
        makedirs(config.iDARTS_DIRECTORY)
        fw = open(resource_config, 'w')
        fw.write('data_download: True\n')
        fw.write('data_dir: ' + outdir + '\n')
        fw.close()
    else:
        logger.info('resources download unsuccessfully')
        resource_config = config.iDARTS_DIRECTORY + 'resource.yaml'
        makedirs(config.iDARTS_DIRECTORY)
        fw = open(resource_config, 'w')
        fw.write('data_download: False\n')
        fw.write('data_dir: ' + outdir + '\n')
        fw.close()