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
# import spikeinterface.qualitymetrics as sqm
# import spikeinterface.comparison as sc
# import spikeinterface.exporters as sexp
# import spikeinterface.curation as scur
# import spikeinterface.widgets as sw

# import numpy as np
# import matplotlib.pyplot as plt
from pathlib import Path
import os
 
global_job_kwargs = dict(n_jobs=8, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)

base_folder = Path('D:\Kelton\probe_data\KW005')

for fold in os.listdir(base_folder):
    sglx_dir = base_folder / fold
    
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
    # Global median reference
    rec = spre.common_reference(rec2, operator='median', reference='global')
    
    # Run KS4 within SpikeInterface
    # ss.get_default_sorter_params('kilosort4')
    sorting = ss.run_sorter('kilosort4', rec, output_folder = base_folder / sglx_dir / 'kilosort4', verbose=True)
    # Load in created params file and change key values
    paramsFile = base_folder / sglx_dir / 'kilosort4/sorter_output/params.py'
    execfile(paramsFile)
    n_channels_dat = rec.get_num_channels() + 1     # Account for Sync channel!
    rec_signals_dict = {e["stream_name"]: e for e in raw_rec.neo_reader.signals_info_list}
    dat_path = rec_signals_dict[raw_rec.stream_id]["bin_file"]  # Write file path to original ap.bin file 

    # Write over params file
    f = open(base_folder / sglx_dir / 'kilosort4/sorter_output/params.py','w')
    writedict = {'n_channels_dat':n_channels_dat, 'offset':offset, 'sample_rate':sample_rate, 'dtype':dtype, 'hp_filtered':hp_filtered, 'dat_path':dat_path}
    for name, val in writedict.items():
        f.write("%s = %s\n" % (name, repr(val)))
    f.close()
