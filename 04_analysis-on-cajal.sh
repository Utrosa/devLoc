#! /usr/bin/env bash
# Time-stamp: <20-10-2026 m.utrosa@bcbl.eu>
# 1st level analysis: task-timDev, GLM
set -eo pipefail

# Find anaconda installation
source /bcbl/home/home_g-m/mutrosa/miniconda3/envs/etc/profile.d/conda.sh

# Activate the environment specific to the project
conda activate nipypee

# Run the python code
python analysis_task-timDev.py
