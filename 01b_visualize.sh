#! /usr/bin/env bash
# Time-stamp: <07-05-2026 m.utrosa@bcbl.eu>

# ---- 1.8 mm isotropics
# acq-BLOCK1
# acq-BLOCK2
# acq-BLOCK3
# acq-BLOCK4
# acq-FUNCLOC

# T1 (MP2RAGE) denoised
python -m scripts.import.visualize_freeview \
  --bids_dir /home/mutrosa/mutrosa/Documents/devLoc/data_MRI/sourcedata/raw/ \
  --sub 05 --ses 02 --modalities bold --acq BLOCK1

