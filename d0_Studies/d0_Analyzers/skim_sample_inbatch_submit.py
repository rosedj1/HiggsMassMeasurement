"""Skim Sample on SLURM in Batches

This code skims a root file in batches,
placing muons which pass selections into a MyMuonCollection
for each event range specified.
Each MyMuonCollection is then pickled ('.pkl').

NOTE:
- All events will be skimmed over, but user can specify events per batch.
- This code uses a template skimmer, like:
d0_Studies/d0_Analyzers/skim_sample_inbatch_template.py
- 
Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-19
Updated:
"""
import os
import shutil
import numpy as np
from ROOT import TFile
from Utils_ROOT.ROOT_fns import get_max_evts
from Utils_Python.SlurmManager import SLURMSubmitter
from Utils_Python.Utils_Files import replace_value, make_dirs, check_overwrite

year = "2016"
prod_mode = "DY2mu"
outfile_prefix = "MC2016DY_skim_fullstats_d0max0p01"
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/DY_JPsi_Upsilon/DY_2016.root"
fullpath_main_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/skim_sample_inbatch_template.py"

outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/output"
outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/copies"
outpkl_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/"

verbose = 1
overwrite = 0
evts_per_batch = 5E6
partition = "hpg2-compute"
mem = (8, 'gb')
nodes = 4
burst = True

def main():
    for d in (outtxt_dir, outcopies_dir, outpkl_dir):
        make_dirs(d, verbose=verbose)
    max_evts = get_max_evts(infile_path, path_to_tree="Ana/passedEvents")
    evt_beg = 0
    while evt_beg < max_evts:
        evt_end = int(evt_beg + evts_per_batch - 1)  # Since starting at event 0.
        if evt_end >= max_evts:
            evt_end = int(max_evts - 1)  # Last event is one fewer than total n_evts.
        # Create main copy.
        base = os.path.basename(fullpath_main_script)
        base = base.replace(".py", f"_copy_{evt_beg}evts{evt_end}")
        fullpath_main_script_copy = os.path.join(outcopies_dir, f"{base}.py")
        check_overwrite(fullpath_main_script_copy, overwrite)
        shutil.copyfile(fullpath_main_script, fullpath_main_script_copy)
        full_file_name = f"{outfile_prefix}_{base}"
        # Modify main copy.
        replace_value("INFILE_PATH",    infile_path,     fullpath_main_script_copy)
        replace_value("PROD_MODE",      prod_mode,       fullpath_main_script_copy)
        replace_value("OUTPKL_DIR",     outpkl_dir,      fullpath_main_script_copy)
        replace_value("FILENAME_FULL",  full_file_name,  fullpath_main_script_copy)
        replace_value("OVERWRITE",      overwrite,       fullpath_main_script_copy)
        replace_value("VERBOSE",        verbose,         fullpath_main_script_copy)
        replace_value("N_EVT_BEG",      evt_beg,         fullpath_main_script_copy)
        replace_value("N_EVT_END",      evt_end,         fullpath_main_script_copy)
        # Prep and submit SLURM script copy.
        sbmtr = SLURMSubmitter(verbose=False)
        sbmtr.prep_directives(
            job_name=full_file_name,
            output_txt=os.path.join(outtxt_dir, f"{full_file_name}.log"),
            email="rosedj1@ufl.edu",
            time="08:00:00",
            acct="avery",
            burst=burst,
            mem=mem,
            partition=partition,
            nodes=nodes,
        )
        cmdtup = (f"time python {fullpath_main_script_copy}")
        fullpath_slurm_copy = os.path.join(outcopies_dir, f"{full_file_name}.sbatch")
        result = sbmtr.make_slurm_script(fullpath_slurm_copy, cmdtup, overwrite=overwrite)
        if result == 0:
            sbmtr.submit_script(fullpath_slurm_copy)
        # Get ready for next batch.
        evt_beg = int(evt_end + 1)

if __name__ == "__main__":
    main()