# -*- coding: utf-8 -*-
"""
Created on Thu Aug 15 16:25:59 2024

@author: spikesorter
"""

# import spikeinterface.full as si
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
#from kilosort.data_tools import (
    #mean_waveform, cluster_templates, get_good_cluster, get_cluster_spikes,
    #get_spike_waveforms, get_best_channels
    #)
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

base_folder = Path('D:\Fox\FG038')

for fold in os.listdir(base_folder):
    sglx_dir = base_folder / fold
    os.chdir(sglx_dir)
        
    try:
        tmp = os.listdir(sglx_dir)
        tmp.index('kilosort4')
        print("kilosort4 folder already found for", sglx_dir)
        continue
    except:
        print("Starting SI kilosort4 analysis for", sglx_dir)
    
    # Load raw data, DON'T load sync channel here or subsequent operations won't work!
    raw_rec = se.read_spikeglx(sglx_dir, stream_name='imec0.ap', load_sync_channel = False)

    # Preprocessing chain
    # High pass filter
    rec1 = spre.highpass_filter(raw_rec, freq_min=400.)
    # Align analog-digital converter sweeps (catGT equivalent)
    rec2 = spre.phase_shift(rec1)
    # Remove bad channels
    #bad_channel_ids, channel_labels = spre.detect_bad_channels(rec1)
    #rec3 = rec2.remove_channels(bad_channel_ids)
    rec3 = rec2;
    #print(f"Bad channel IDs {bad_channel_ids}")
    # Global median reference
    rec = spre.common_reference(rec3, operator='median', reference='global')
    
    subF = fold[:-3]        #Accounts for _g0 tag
    # subprocess.run(f"CatGT -dir={base_folder} -run={subF} -g=0 -t=0 -prb_fld -lf -lffilter=butter,12,0,300 -prb=0 -save=2,0,0,384", shell=True)   #Export sync chan only for size reasons
    # subprocess.run(f"CatGT -dir={base_folder} -run={subF} -g=[0,1] -t=0 -prb_fld -lf -lffilter=butter,12,0,300 -prb=0 -save=2,0,0,384 -save=2,0,1,0:383", shell=True)       #Export both sync chan and LFP albeit separately
    subprocess.run(f"CatGT -dir={base_folder} -run={subF} -g=[0,1] -t=0 -prb_fld -ap -prb=0 -apfilter=butter,12,300,9000 -dest=D:\Fox\FG038\combofolderD3\test -save=2,0,0,384 -save=2,0,1,0:383", shell=True)       #Export both sync chan and AP albeit separately
         
