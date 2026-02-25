# -*- coding: utf-8 -*-
"""
Created on Thu May 15 11:26:07 2025

@author: spikesorter
"""

import spikeinterface as si  # import core only
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
import spikeinterface.sorters as ss
# import spikeinterface.postprocessing as spost
# import spikeinterface.qualitymetrics as sqm bbvb   
# import spikeinterface.comparison as sc
# import spikeinterface.exporters as sexp
# import spikeinterface.curation as scur
# import spikeinterface.widgets as sw

import subprocess  
import numpy as np
from kilosort.io import load_ops
from kilosort.data_tools import (
    mean_waveform, cluster_templates, get_good_cluster, get_cluster_spikes,
    get_spike_waveforms, get_best_channels
    )
import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams
rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False
gray = .5 * np.ones(3)
from pathlib import Path
from datetime import datetime
import os
 
global_job_kwargs = dict(n_jobs=8, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)

base_folder = Path('D:\Kelton\probe_data\KW041')

for fold in os.listdir(base_folder):
    sglx_dir = base_folder / fold

    results_dir = base_folder / sglx_dir / 'kilosort4/sorter_output'
    ops = load_ops(results_dir / 'ops.npy')
    # t = (np.arange(ops['nt']) / ops['fs']) * 1000
    fig = plt.figure(figsize=(12,5), dpi=100)
    grid = gridspec.GridSpec(1, 1, figure=fig, hspace=0.5, wspace=0.5)
    
    ax = fig.add_subplot(grid[0,0])
    ax.plot(np.arange(0, ops['Nbatches'])*2, ops['dshift']);
    ax.set_xlabel('time (sec.)')
    ax.set_ylabel('drift (um)')
    fig.savefig(base_folder / sglx_dir / 'kilosort4/drift.png')
    plt.close('all')
