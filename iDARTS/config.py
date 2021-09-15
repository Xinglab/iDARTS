# -*- coding: UTF-8 -*-

"""Here are general configurations for the Darts_DNN package, including 
version control, trained model parameter, etc.
"""

from pkg_resources import resource_filename
import os
import yaml

CURRENT_VERSION = "v1.0"
HG19_GENOME = resource_filename('iDARTS.resources', 'hg19.genome')
MODEL_DIR = resource_filename('iDARTS.resources', 'model')
RBP_LIST_FN = resource_filename('iDARTS.resources', 'RBP.GeneName.txt')
DOWNLOAD_CONFIG = yaml.safe_load(open(resource_filename('iDARTS.resources', 'download.yaml'), 'r'))
CHECK_MD5_CONFIG = yaml.safe_load(open(resource_filename('iDARTS.resources', 'check_md5.yaml'), 'r'))
iDARTS_DIRECTORY = os.path.expanduser('~') + '/.iDARTS/'