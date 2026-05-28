#! /usr/bin/env python
# Time-stamp: <2026-04-29 m.utrosa@bcbl.eu>
# conda activate localizer_fMRI
'''
Calculates participant's performance from logfiles (events.tsv) and input parameter
combination .csv, which in necessary for compensation.
'''

# ---------- PREP
import ast
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt

# ---------- FUNCTIONS
def extract_design_data(trials_dir, combo, block, target_column):
    
    # Load the csv with trials
    filename = f'ses-{combo:003d}_exp_parameter_combo.csv'
    df_raw = pd.read_csv(Path(trials_dir) / filename)

    # Correct column data types: lists should not be strings.
    # TODO: correct for it when creating the dataframe. Not here.
    list_cols = ["freq_dev", "freq_dev_type", "freq_loc", "freq_diff", "freq_diff_abs"]
    for col in list_cols:
      df_raw[col] = df_raw[col].apply(
          lambda x: ast.literal_eval(x) if isinstance(x, str) else x
      )

    # Remove silent trials
    df_raw.dropna(inplace=True)

    # Select the part of the dataframe relevant for the current block
    df_block = df_raw[df_raw["block_no"] == block]

    # Reset index to enable comparison with subject's dataframe later
    df_block.reset_index(drop=True, inplace=True)

    return df_block[target_column]

def extract_participant_data(project_dir, subject, session, task, block):

   # Find files for the given task
   project_path = Path(project_dir)
   files = [f for f in project_path.rglob("*") if f.is_file() and task in f.name]

   # Sort the files in ascending order based on the timestamp
   files.sort()

   # Select the file corresponding to the block
   file = files[block - 1]

   # Load the file by remove the first to lines (software info & date stamp)
   df = pd.read_csv(str(file), sep='\t', skiprows=[0,1])
   
   # Return participant's resonses
   return df["trial_type"]

# ---------- USAGE
if __name__ == "__main__":
    no_blocks = int(input("Enter the number of blocks:"))
    comboID   = int(input("Enter the index of exp_parameter_combo.csv:"))
    sesID = int(input("Enter the session ID:"))
    subID = int(input("Enter the subject ID:"))

    per_corr_exp = []
    for i in range(no_blocks):
        block_idx = i + 1

        # Get responses by design and by participant
        freq_dev_no_design = extract_design_data(
            trials_dir = "/home/mutrosa/Documents/projects/auditory_paradigms/detection_accuracy/selected_trials",
            combo = comboID,
            block = block_idx,
            target_column = "freq_dev_no")

        freq_dev_no_participant = extract_participant_data(
            project_dir = f"/home/mutrosa/Documents/projects/devLoc/data_logs/bids_output_{comboID:003d}",
            subject = subID,
            session = sesID,
            task = "freqDev",
            block = block_idx)

        # Calculate percentage correct per block
        # Make sure that the data is of the same type
        resp_design_int  = freq_dev_no_design.astype(int)
        resp_design_str  = resp_design_int.astype(str)
        resp_subject_str = freq_dev_no_participant.astype(str)

        # Count no. of correct responses
        no_same = resp_design_str == resp_subject_str
        no_corr = no_same.sum() # True values are summed
        no_all  = len(freq_dev_no_design)
        per_corr_block = no_corr / no_all
        per_corr_exp.append(per_corr_block)

        # Print block performance in the terminal
        print(f"\nPercentage correct in block no. {block_idx}: {round(per_corr_block*100)} %")

    # Print update about overall performance
    print(f"\n*** Percentage correct in the experiment: {round((sum(per_corr_exp) / len(per_corr_exp))*100)} % ***")

# TODO: ADD A CALCULATION FOR HOW MUCH EXTRA MONEY DOES THIS MEAN