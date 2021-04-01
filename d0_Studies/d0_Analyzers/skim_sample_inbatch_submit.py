"""Skim Sample on SLURM in Batches

This code skims a root file in batches,
placing muons which pass selections into a MyMuonCollection
for each event range specified.
Each MyMuonCollection is then pickled ('.pkl').

NOTE:
- All events will be skimmed over, but user can specify events per batch.
- This code uses a template skimmer:
d0_Studies/d0_Analyzers/skim_sample_inbatch_template.py

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-19
Updated: 2021-03-31
"""
import os
import shutil
import numpy as np
from ROOT import TFile
from Utils_ROOT.ROOT_fns import get_max_evts
from Utils_Python.SlurmManager import SLURMSubmitter
from Utils_Python.Utils_Files import replace_value, make_dirs, check_overwrite

year = "2018"
prod_mode = "DY2mu"
# outfile_prefix = "MC2016DY_skim_fullstats_0p01_d0_1000p0_nogenmatching"
outfile_prefix = "MC2018DY_skim_fullstats_nogenmatching"
# infile_path = "/cmsuf/data/store/user/t2/users/ferrico/SingleBS_studies/After/DY_2016_0.root"
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/SingleBS_studies/After/DY_2018_2.root"
# infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/DY_JPsi_Upsilon/DY_2016.root"
fullpath_main_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/skim_sample_inbatch_template.py"

outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY/output/"
outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY/copies/"
outpkl_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2018/DY/"

# Event selection.
inv_m_lim = [60.0, 120.0]
eta_lim = [0.0, 2.4]
pT_lim = [5.0, 1000.0]
d0_lim = [0.0, 1000.0]  # No d0 cuts now. Apply d0 cuts when separating muons into KB2Ds.
dR_max = 100 # Make it huge to investigate no gen matching. 0.001

# If `None` then run over entire file.
specified_evt_beg = None # 122000000
specified_evt_end = None

verbose = 1
overwrite = 0
evts_per_batch = 2E6  # The smaller, the faster it will be processed on SLURM.
# SLURM directives.
partition = "hpg2-compute"
mem = (8, 'gb')
nodes = 4
burst = True

def main():
    # Make a subdir for each area.
    new_outtxt_dir = os.path.join(outtxt_dir, outfile_prefix)
    new_outcopies_dir = os.path.join(outcopies_dir, outfile_prefix)
    new_outpkl_dir = os.path.join(outpkl_dir, outfile_prefix)
    assert all(year in f for f in (outfile_prefix, infile_path, new_outtxt_dir, new_outcopies_dir, new_outpkl_dir))
    for d in (new_outtxt_dir, new_outcopies_dir, new_outpkl_dir):
        make_dirs(d, verbose=verbose)
    max_evts = get_max_evts(infile_path, path_to_tree="Ana/passedEvents") if specified_evt_end is None else specified_evt_end
    evt_beg = specified_evt_beg if specified_evt_beg is not None else 0
    while evt_beg < max_evts:
        evt_end = int(evt_beg + evts_per_batch - 1)  # Since starting at event 0.
        if evt_end >= max_evts:
            evt_end = int(max_evts - 1)  # Last event is one fewer than total n_evts.
        # Create main copy.
        base = os.path.basename(fullpath_main_script)
        base = base.replace(".py", f"_copy_{evt_beg}evts{evt_end}")
        fullpath_main_script_copy = os.path.join(new_outcopies_dir, f"{base}.py")
        check_overwrite(fullpath_main_script_copy, overwrite)
        shutil.copyfile(fullpath_main_script, fullpath_main_script_copy)
        full_file_name = f"{outfile_prefix}_{base}"
        # Modify main copy.
        replace_value("INFILE_PATH",   infile_path,    fullpath_main_script_copy)
        replace_value("PROD_MODE",     prod_mode,      fullpath_main_script_copy)
        replace_value("OUTPKL_DIR",    new_outpkl_dir, fullpath_main_script_copy)
        replace_value("FILENAME_FULL", full_file_name, fullpath_main_script_copy)
        replace_value("OVERWRITE",     overwrite,      fullpath_main_script_copy)
        replace_value("VERBOSE",       verbose,        fullpath_main_script_copy)
        replace_value("N_EVT_BEG",     evt_beg,        fullpath_main_script_copy)
        replace_value("N_EVT_END",     evt_end,        fullpath_main_script_copy)
        replace_value("INV_M_LIM",     inv_m_lim,      fullpath_main_script_copy)
        replace_value("ETA_LIM",       eta_lim,        fullpath_main_script_copy)
        replace_value("PT_LIM",        pT_lim,         fullpath_main_script_copy)
        replace_value("D0_LIM",        d0_lim,         fullpath_main_script_copy)
        replace_value("DR_MAX",        dR_max,         fullpath_main_script_copy)
        # Prep and submit SLURM script copy.
        sbmtr = SLURMSubmitter(verbose=False)
        sbmtr.prep_directives(
            job_name=full_file_name,
            output_txt=os.path.join(new_outtxt_dir, f"{full_file_name}.log"),
            email="rosedj1@ufl.edu",
            time="08:00:00",
            acct="avery",
            burst=burst,
            mem=mem,
            partition=partition,
            nodes=nodes,
        )
        cmdtup = (f"time python {fullpath_main_script_copy}")
        fullpath_slurm_copy = os.path.join(new_outcopies_dir, f"{full_file_name}.sbatch")
        result = sbmtr.make_slurm_script(fullpath_slurm_copy, cmdtup, overwrite=overwrite)
        if result == 0:
            sbmtr.submit_script(fullpath_slurm_copy)
        # Get ready for next batch.
        evt_beg = int(evt_end + 1)

if __name__ == "__main__":
    main()