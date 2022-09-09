# -*- coding: UTF-8 -*-

"""
iDARTS - build_feature_utils
"""

import os
import logging
from subprocess import *
import string
import yaml
import sys
from multiprocessing import Pool
import psutil

from . import config

logger = logging.getLogger('iDARTS.buile_feature')

def parser(args):
    resources_download_flag = False
    if os.path.exists(config.iDARTS_DIRECTORY + 'resource.yaml'):
        data_config = yaml.safe_load(open(config.iDARTS_DIRECTORY + 'resource.yaml', 'r'))
        if 'data_download' in data_config and data_config['data_download']:
            resources_download_flag = True
    if resources_download_flag:
        event_type = args.event_type
        if event_type == 'SE':
            from . import build_feature_SE
            build_feature_SE.main(args)
        if event_type == 'A3SS' or event_type == 'A5SS':
            from . import build_feature_ASS
            build_feature_ASS.main(args)
    else:
        logger.info('resources need to be downloaded first\nUsage:\n...iDARTS get_resources')
        sys.exit(0)