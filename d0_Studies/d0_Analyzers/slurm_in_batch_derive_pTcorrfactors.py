"""SLURM Script Duplicator and Submitter

This code makes a copy of a "main template" script that you wish to run on
SLURM. It actually makes many copies of the main script and of the SLURM
submissions script, each copy only differing from the next by the eta bin
used.

User can put in a whole list of eta values, and each bin will be considered:

Example: full_eta_ls = [0.2, 0.4, 0.6, 0.8]
Then this code will produce a copy of the main script
and of the SLURM submission script using eta_bin = [0.2, 0.4].
The next copy will use eta_bin = [0.4, 0.6], etc.

NOTE:
    The main template script should have replacement strings that start with
    "REPLACE_". Check the `replace_vals_in_files` function to see what values
    get replaced.

Requires an input pickled dict of KinBin2Ds.
- You can make the input dict from roch_vs_noroch_kb2dmaker_inclusivehistplotter.py

Author: Jake Rosenzweig
Updated: 2021-03-10
"""
from Utils_Python.Utils_Files import replace_value, make_dirs
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
# from d0_Studies.d0_Utils.d0_cls import OrganizerKB2D
import subprocess
import shutil
import os
import sys

#-----------------------#
#--- User Parameters ---#
#-----------------------#
# job_name = "RC_vs_NoRC_itergaussfits_fullstats_pT75then200GeV_extendedxaxis"
overwrite = 1
iters = 5
regions = 12
verbose = 0
max_n_evts = -1
print_out_every = 100000

fit_whole_range_first_iter = False  # False gives more consistent fits (with no outlier data).
use_data_in_xlim = True
binned_fit = False  # Bin q*d0 axis in: dpT/pT vs. q*d0 plot.
switch_to_binned_fit = 5000
min_muons_per_qd0_bin = 3000

job_name = "MC2016DY_deriveAdHocpTcorrfactors_fullstats"

template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_template.py"
# template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py"
template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_slurm.sbatch"

inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_fullstats_muoncoll_withkb3dbins.pkl"
# inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/MC2018ggH_KB2D_onlymuoninfo_fullstats_pTbinsroundnumbers_75then200GeV.pkl"

outdir_pkl    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/"
outdir_pdf    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/"
outdir_copies = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/scriptcopies/"
outdir_txt    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/output/"

full_eta_ls = equal_entry_bin_edges_eta_mod1_wholenum[7:]  #[2.1, 2.2]
# full_pT_ls = [30.0, 40.0, 50.0]  #[2.1, 2.2]
full_pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum #[50.0, 75.0, 100.0, 200.0]  #[2.1, 2.2]

delete_copies = False
#------------------------#
#--- Script functions ---#
#------------------------#
def make_name_from_ls(ls, val):
    """Parse a 2-element list of eta values into a title string."""
    name = f"{ls[0]}{val}{ls[-1]}"
    name = name.replace(".", "p")
    return name

def prep_area(dir_tup):
    """Make sure all directories exist. If not, create them."""
    for d in dir_tup:
        make_dirs(d)

def make_filepaths_of_copies(eta_ls, pT_ls, template_main, template_slurm, outdir_copies):
    """Return the full paths of copies of the main script and the slurm script.
    
    Parameters
    ----------
    outdir_copies : str
        Path to directory where copies of scripts will be stored.
    """
    eta_name = make_name_from_ls(eta_ls, "eta")
    pT_name = make_name_from_ls(pT_ls, "pT")
    suffix = f"_{eta_name}_{pT_name}.py"
    filename_main = os.path.split(template_main)[1]
    filename_slurm = os.path.split(template_slurm)[1]
    fullpath_copy_main = os.path.join(outdir_copies, filename_main.replace(".py", f"_copy_{eta_name}_{pT_name}.py"))
    fullpath_copy_slurm = os.path.join(outdir_copies, filename_slurm.replace(".sbatch", f"_copy_{eta_name}_{pT_name}.sbatch"))
    return (fullpath_copy_main, fullpath_copy_slurm)

def replace_vals_in_files(eta_ls, pT_ls, job_name, fullpath_new_main_script,
                          outdir_pkl, outdir_pdf, outdir_txt,
                          template_tup=()):
    """Replace values in scripts corresponding to eta_ls and submit new SLURM script.

    NOTE: Works for a 2-element eta_ls: [eta_min, eta_max]

    Parameters
    ----------
    template_tup : tuple
        Contains file paths to all template files that will have their capital
        words replaced.
    """
    eta_name = make_name_from_ls(eta_ls, "eta")
    pT_name = make_name_from_ls(pT_ls, "eta")
    # Replace phrases in slurm script.
    for template in template_tup:
        replace_value("REPLACE_OVERWRITE", overwrite, template)
        replace_value("REPLACE_ITERS", iters, template)
        replace_value("REPLACE_REGIONS", regions, template)
        replace_value("REPLACE_VERBOSE", verbose, template)
        replace_value("REPLACE_fit_whole_range_first_iter", fit_whole_range_first_iter, template)
        replace_value("REPLACE_use_data_in_xlim", use_data_in_xlim, template)
        replace_value("REPLACE_binned_fit", binned_fit, template)
        replace_value("REPLACE_switch_to_binned_fit", switch_to_binned_fit, template)
        replace_value("REPLACE_min_muons_per_qd0_bin", min_muons_per_qd0_bin, template)
        replace_value("REPLACE_max_n_evts", max_n_evts, template)
        replace_value("REPLACE_print_out_every", print_out_every, template)
        replace_value("REPLACE_ETA_LS", eta_ls, template)
        replace_value("REPLACE_PT_LS", pT_ls, template)
        replace_value("REPLACE_ETA_NAME", eta_name, template)
        replace_value("REPLACE_PT_NAME", pT_name, template)
        replace_value("REPLACE_JOB_NAME", job_name, template)
        replace_value("REPLACE_NEW_SCRIPT", fullpath_new_main_script, template)
        replace_value("REPLACE_INPATH_PKL", inpath_pkl, template)
        replace_value("REPLACE_OUTDIR_PKL", outdir_pkl.rstrip('/'), template)
        replace_value("REPLACE_OUTDIR_PDF", outdir_pdf.rstrip('/'), template)
        replace_value("REPLACE_OUTDIR_TXT", outdir_txt.rstrip('/'), template)

def print_info(fullpath_copy_main_script, fullpath_copy_slurm_script, 
                outdir_copies, outdir_pkl, outdir_txt, outdir_pdf):
    """Print the filepath locations."""
    print(f"[INFO] Main script copy created:\n{fullpath_copy_main_script}")
    print(f"[INFO] SLURM script copy created:\n{fullpath_copy_slurm_script}")
    print(f"[INFO] Dir to copies: {outdir_copies}")
    print(f"[INFO] Dir to pickle: {outdir_pkl}")
    print(f"[INFO] Dir to txt:    {outdir_txt}")
    print(f"[INFO] Dir to pdf:    {outdir_pdf}")

if __name__ == "__main__":
    for eta_min, eta_max in zip(full_eta_ls[:-1], full_eta_ls[1:]):
    # for eta_min, eta_max in zip(full_eta_ls[:-1], full_eta_ls[1:]):
        # print(f"...Preparing work area.")
        prep_area([outdir_copies, outdir_pkl, outdir_txt, outdir_pdf])
        eta_ls = [eta_min, eta_max]
        fullpath_copy_main_script, fullpath_copy_slurm_script = make_filepaths_of_copies(eta_ls, full_pT_ls, template_script_main, template_script_slurm, outdir_copies)
        print_info(fullpath_copy_main_script, fullpath_copy_slurm_script, 
                   outdir_copies, outdir_pkl, outdir_txt, outdir_pdf)
        
        shutil.copyfile(template_script_main, fullpath_copy_main_script)
        shutil.copyfile(template_script_slurm, fullpath_copy_slurm_script)
        replace_vals_in_files(eta_ls, full_pT_ls, job_name, fullpath_copy_main_script,
                          outdir_pkl, outdir_pdf, outdir_txt,
                          template_tup=(fullpath_copy_main_script, fullpath_copy_slurm_script))
        print(f"...Submitting SLURM script for eta range: {eta_ls}")
        output = subprocess.run(["sbatch", fullpath_copy_slurm_script])
        if delete_copies:
            print(f"[INFO] Removing copies of scripts.")
            os.remove(new_main)
            os.remove(new_slurm)
