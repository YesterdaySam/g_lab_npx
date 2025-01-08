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

import subprocess
import numpy as np
# import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
import os
 
global_job_kwargs = dict(n_jobs=8, chunk_duration="1s")
si.set_global_job_kwargs(**global_job_kwargs)

base_folder = Path('D:\Kelton\probe_data\KW022')

for fold in os.listdir(base_folder):
    sglx_dir = base_folder / fold
    os.chdir(sglx_dir)
     
    # try:
    #     tmp = os.listdir(sglx_dir)
    #     tmp.index('kilosort4')
    #     print("kilosort4 folder already found for", sglx_dir)
    #     continue
    # except:
    #     print("Starting SI kilosort4 analysis for", sglx_dir)
    
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
    
    # Run KS4 within SpikeInterface
    # ss.get_default_sorter_params('kilosort4')
    if np.size(np.unique(rec.get_channel_groups())) == 1:
        # Handle 1-shank probes
        sorting = ss.run_sorter('kilosort4', rec, output_folder = base_folder / sglx_dir / 'kilosort4', verbose=True)
        paramsFile = base_folder / sglx_dir / 'kilosort4/sorter_output/params.py'
    else:
        # Handle multishank probes by calling CatGT to save extra .lf.bin file with sync chanel for later extraction by matlab without capping RAM
        # sorting = ss.run_sorter_by_property('kilosort4', rec, grouping_property='group', folder=base_folder / sglx_dir / 'kilosort4', verbose=True)
        sorting = ss.run_sorter('kilosort4', rec, output_folder = base_folder / sglx_dir / 'kilosort4', verbose=True)
        paramsFile = base_folder / sglx_dir / 'kilosort4/sorter_output/params.py'
        
        os.chdir(r'C:\Users\spikesorter\Documents\CatGTWinApp\CatGT-win')
        subF = fold[:-3]        #Accounts for _g0 tag
        # subprocess.run(f"CatGT -dir={base_folder} -run={subF} -g=0 -t=0 -prb_fld -lf -lffilter=butter,12,0,300 -prb=0 -save=2,0,0,384", shell=True)   #Export sync chan only for size reasons
        subprocess.run(f"CatGT -dir={base_folder} -run={subF} -g=0 -t=0 -prb_fld -lf -lffilter=butter,12,0,300 -prb=0 -save=2,0,0,384 -save=2,0,10,0:383", shell=True)       #Export both sync chan and LFP albeit separately
         
        # Experimental SI version of saving sync channel, don't use because can't read out later
        # sync_rec = se.read_spikeglx(sglx_dir, stream_name='imec0.ap', load_sync_channel = True)
        # sync_rec = sync_rec.remove_channels([sync_rec.channel_ids[:-1]]) #Get sync channel
        # job_kwargs = dict(n_jobs=-1, chunk_duration='1s', progress_bar=True)
        # sync_rec.save(folder=base_folder / sglx_dir / 'sync_out_NPX2', format='binary', **job_kwargs)
        
        # paramsFile = base_folder / sglx_dir / 'kilosort4/0/sorter_output/params.py'
    
    # Load in created params file and change key values
    # paramsFile = base_folder / sglx_dir / 'kilosort4/sorter_output/params.py'
    execfile(paramsFile)
    n_channels_dat = raw_rec.get_num_channels() + 1     # Account for Sync channel! And use raw_rec not rec after bad chan removal, otherwise Phy reads in wrong channel #
    rec_signals_dict = {e["stream_name"]: e for e in raw_rec.neo_reader.signals_info_list}
    dat_path = rec_signals_dict[raw_rec.stream_id]["bin_file"]  # Write file path to original ap.bin file 

    # Write over params file
    f = open(base_folder / sglx_dir / 'kilosort4/sorter_output/params.py','w')
    writedict = {'n_channels_dat':n_channels_dat, 'offset':offset, 'sample_rate':sample_rate, 'dtype':dtype, 'hp_filtered':hp_filtered, 'dat_path':dat_path}
    for name, val in writedict.items():
        f.write("%s = %s\n" % (name, repr(val)))
    f.close()
    
    print("Finished sort at " + str(datetime.now()))

os.chdir(base_folder)