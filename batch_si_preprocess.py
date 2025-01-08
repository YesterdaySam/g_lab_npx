# -*- coding: utf-8 -*-
"""
Created on Tue Sep 17 10:02:25 2024

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

base_folder = Path('D:\Fox\FG021')

for fold in os.listdir(base_folder):
    sglx_dir = base_folder / fold
    
    try:
        tmp = os.listdir(sglx_dir)
        tmp.index('preprocess')
        print("binary output folder already found for", sglx_dir)
        continue
    except:
        print("Starting SI preprocessing for", sglx_dir)
        
    # Load raw data, DON'T load sync channel here or subsequent operations won't work!
    raw_rec = se.read_spikeglx(sglx_dir, stream_name='imec0.ap', load_sync_channel = False)

    # Preprocessing chain
    # High pass filter
    rec1 = spre.highpass_filter(raw_rec, freq_min=400.)
    # Align analog-digital converter sweeps (catGT equivalent)
    rec2 = spre.phase_shift(rec1)
    # Global median reference
    rec = spre.common_reference(rec2, operator='median', reference='global')
    
    # Run KS4 within SpikeInterface
    # ss.get_default_sorter_params('kilosort4')
    if np.size(np.unique(rec.get_channel_groups())) == 1:
        # Ignore 1-shank probes
        print("Skipping preprocessing saves for NPX 1.0")
        continue
    else:
        # Handle multishank probes or other groupings
        job_kwargs = dict(n_jobs=-1, chunk_duration='1s', progress_bar=True)
        rec.save(folder=base_folder / sglx_dir / 'preprocess', format='binary', **job_kwargs)
