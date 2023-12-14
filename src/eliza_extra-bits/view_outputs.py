
import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
pd.options.mode.chained_assignment = None  # default='warn'

parent_dir = "/data/eliza/pytcpl/"
sys.path.append(parent_dir)

from src.pipeline.pipeline_constants import ROOT_DIR, AEIDS_LIST_PATH
from src.pipeline.pipeline_helper import load_config

config, config_path = load_config()
global CONFIG
CONFIG = config

with open(AEIDS_LIST_PATH, 'r') as file:
    aeids_list = [line.strip() for line in file]

path = os.path.join(ROOT_DIR, '../../data/output/54.parquet.gzip')
sample_df = pd.read_parquet(path_sample)