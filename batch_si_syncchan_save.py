# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 16:43:39 2024

@author: spikesorter
"""

import spikeinterface as si  # import core only
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
import spikeinterface.sorters as ss

import numpy as np
# import matplotlib.pyplot as plt
from pathlib import Path
import os

global_job_kwargs = dict(n_jobs=8, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)

base_folder = Path('D:\Kelton\probe_data\KW012')

for fold in os.listdir(base_folder):
    sglx_dir = base_folder / fold
    
    try:
        tmp = os.listdir(sglx_dir)
        tmp.index('.lf.bin')
        print("LFP output folder already found for", sglx_dir)
        continue
    except:
        print("Starting SI preprocessing for", sglx_dir)
        
    # Load raw data, DON'T load sync channel here or subsequent operations won't work!
    sync_rec = se.read_spikeglx(sglx_dir, stream_name='imec0.ap', load_sync_channel = False)
    if np.size(np.unique(sync_rec.get_channel_groups())) == 1:
        # Ignore 1-shank probes
        print("Skipping preprocessing saves for NPX 1.0")
        continue
    else:
        sync_rec = se.read_spikeglx(sglx_dir, stream_name='imec0.ap', load_sync_channel = True)
        sync_rec = sync_rec.remove_channels([sync_rec.channel_ids[:-1]]) #Get sync channel
        job_kwargs = dict(n_jobs=-1, chunk_duration='1s', progress_bar=True)
        sync_rec.save(folder=base_folder / sglx_dir / 'sync_out_NPX2', format='binary', **job_kwargs)
