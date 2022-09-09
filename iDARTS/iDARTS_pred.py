# -*- coding: UTF-8 -*-

"""
iDARTS - predict
Predict splicing levels given sequence features and expression profiles
"""

import os
import h5py
from tqdm import tqdm
import numpy as np
from pkg_resources import resource_filename
from . import config
from keras.models import model_from_json

import logging
logger = logging.getLogger('iDARTS.predict')

def get_file_num_lines(filename):
    fp = open(filename)
    count = 0
    for line in fp:
        count += 1
    fp.close()
    return count

def load_dl_model(event_type):
    model_list = []
    model_config_fopen = open(config.MODEL_DIR + '/model.{}.json'.format(event_type))
    model_config = model_config_fopen.read()
    model_config_fopen.close()
    for model_num in xrange(1, 6):
        model = model_from_json(model_config)
        model.load_weights(config.MODEL_DIR + '/model.{}.{}.hd5'.format(event_type, model_num))
        model_list.append(model)
    return model_list

def get_dnn_pred(X, model_list):
    y = np.mean(np.array([model.predict(X).flatten() for model in model_list]), axis = 0)
    return y

def get_rbp_list():
    rbp_list = []
    fp = open(config.RBP_LIST_FN)
    fp.readline()
    for line in fp:
        sp = line.strip().split('\t')
        rbp_list.append(sp[0])
    fp.close()
    return rbp_list

def get_user_exp_profiles(expr_fn, exp_max_value):
    rbp_list = get_rbp_list()
    rbp_dict = {rbp:True for rbp in rbp_list}
    fp = open(expr_fn)
    header = fp.readline().strip().split('\t')
    tissue_list = header[1].split(',')
    tissue_num = len(tissue_list)
    dict_rbp2exp = {}
    for line in fp:
        sp = line.strip().split('\t')
        gene_id = sp[0].split('.')[0]
        if gene_id in rbp_dict:
            dict_rbp2exp[gene_id] = [float(x) for x in sp[1].split(',')]
    fp.close()
    if len(dict_rbp2exp) != len(rbp_list):
        logger.info("RBP expression profile generated failed, please check RBP Gene lists")
        sys.stderr.write("RBP expression profile generated failed, please check RBP Gene lists")
        sys.exit(0)
    tissue2exp_matrix = []
    for rbp_id in rbp_list:
        tissue2exp_matrix.append(dict_rbp2exp[rbp_id])
    tissue2exp_matrix = np.array(tissue2exp_matrix)
    tissue2exp_matrix = tissue2exp_matrix.T / exp_max_value
    user_exp_dict = {}
    user_exp_dict['tissue2exp_matrix'] = tissue2exp_matrix
    user_exp_dict['tissue_list'] = tissue_list
    user_exp_dict['tissue_num'] = tissue_num
    return user_exp_dict

def make_prediction(args):
    event_type = args.event_type
    input_file = args.input
    out_file = args.output
    model_list = load_dl_model(event_type)
    event_config = h5py.File(config.MODEL_DIR + '/{}.config.hd5'.format(event_type), 'r')
    feature_max_value = event_config.get('feature_max_value')[:]
    if args.expr != None:
        exp_max_value = event_config.get('exp_max_value')[:]
        expr_fn = args.expr
        user_exp_dict = get_user_exp_profiles(expr_fn, exp_max_value)
        tissue_list = user_exp_dict['tissue_list']
        tissue_num = user_exp_dict['tissue_num']
        tissue2exp_matrix = user_exp_dict['tissue2exp_matrix']
    else:
        tissue_list = list(event_config.get('tissue_list')[:])
        tissue_num = len(tissue_list)
        tissue2exp_matrix = event_config.get('tissue2exp_matrix')[:]
    num_lines = get_file_num_lines(input_file) - 1
    fp = open(input_file)
    header = fp.readline().strip().split('\t')
    if event_type == 'A3SS' or event_type == 'A5SS':
        feature_start_idx = header.index('9G8.C_F')
    else:
        feature_start_idx = header.index('LogLen.C1')
    fw = open(out_file, 'w')
    fw.write(header[0] + '\t' + '\t'.join(tissue_list) + '\n')
    for line in tqdm(fp, total = num_lines):
        sp = line.strip().split('\t')
        feature = np.array([float(x) for x in sp[feature_start_idx:]]) / feature_max_value
        X = np.concatenate((np.repeat([feature], tissue_num, axis = 0), tissue2exp_matrix), axis = 1)
        y = get_dnn_pred(X, model_list)
        fw.write(sp[0] + '\t' + '\t'.join([str(x) for x in y]) + '\n')
    fp.close()

def parser(args):
    make_prediction(args)