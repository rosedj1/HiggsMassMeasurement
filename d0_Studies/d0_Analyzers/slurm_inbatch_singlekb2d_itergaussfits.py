"""SLURM Script Duplicator and Submitter

This code submits a N SLURM jobs,
where N = number of pT bins in 1 eta bin.

#--- OLD DOCSTRING BELOW ---#

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
Updated: 2021-03-15
"""
from Utils_Python.Utils_Files import replace_value, make_dirs
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
# from d0_Studies.d0_Utils.d0_cls import OrganizerKB2D
import subprocess
import shutil
import os
import sys
#-----------------------#
#--- User Parameters ---#
#-----------------------#
overwrite = 0
iters = 5
verbose = 0
fit_whole_range_first_iter = False  # False gives more consistent fits (with no outlier data).
use_data_in_xlim = 1
binned_fit = False
switch_to_binned_fit = 999999

job_name_base = "MC2016DY_individKB2D_withitergaussfitsonKB3Ds"  # Also prefix for outfile.
eta_range = equal_entry_bin_edges_eta_mod1_wholenum[12:]
full_pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum[:-1]  # I applied a cut of pT < 200 GeV. Oops.

template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_singlekb2d_itergaussfits_template.py"
# template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py"
template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_to_slurm_template.sbatch"
# template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/RochCorrAnalyzers/roch_vs_noroch_slurm.sbatch"

inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/kb2d_dictsnofitinfo/MC2016DY_fullstats_muoncoll_withkb3dbins__ETAPART_PTPART.pkl"

# /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/output/DeriveCorrections/2016DY/individualKB2Ditergaussfits/MC2016DY_kb2d_ETARANGE_PTRANGE_error.log
outdir_copies = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/scriptcopies/{job_name_base}"
outdir_pkl    = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/{job_name_base}"
outdir_txt    = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/output/{job_name_base}"
outdir_pdf    = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/{job_name_base}"

#------------------------#
#--- Script functions ---#
#------------------------#
def prep_area(dir_tup):
    """Make sure all directories exist. If not, create them."""
    for d in dir_tup:
        make_dirs(d)

def make_filepaths_of_copies(eta_name, pT_name, template_main, template_slurm, outdir_copies):
    """Return the full paths of copies of the main script and the slurm script.
    
    Parameters
    ----------
    outdir_copies : str
        Path to directory where copies of scripts will be stored.
    """
    suffix = f"_{eta_name}_{pT_name}"
    filename_main = os.path.split(template_main)[1]
    filename_slurm = os.path.split(template_slurm)[1]
    fullpath_copy_main = os.path.join(outdir_copies, filename_main.replace(".py", f"_copy_{suffix}.py"))
    fullpath_copy_slurm = os.path.join(outdir_copies, filename_slurm.replace(".sbatch", f"_copy_{suffix}.sbatch"))
    return (fullpath_copy_main, fullpath_copy_slurm)

def replace_vals_in_files(eta_range, pT_range, inpkl_path, outpkl_path,
                          job_name_base, fullpath_new_main_script,
                          eta_name, pT_name,
                          template_tup=()):
    """Replace values in scripts corresponding to eta_range and submit new SLURM script.

    NOTE: Works for a 2-element eta_range: [eta_min, eta_max]

    Parameters
    ----------
    template_tup : tuple
        Contains file paths to all template files that will have their capital
        words replaced.
    """
    # Replace phrases in slurm script.
    for template in template_tup:
        # In argument order.
        replace_value("REPLACE_ETA_LS", eta_range, template)
        replace_value("REPLACE_PT_RANGE", pT_range, template)
        replace_value("REPLACE_INPKL_PATH", inpkl_path, template)
        replace_value("REPLACE_OUTPATH_PKL", outpkl_path, template)
        replace_value("REPLACE_JOB_NAME", job_name_base, template)
        replace_value("REPLACE_OUTDIR_TXT", outdir_txt.rstrip('/'), template)
        replace_value("REPLACE_NEW_SCRIPT", fullpath_new_main_script, template)
        replace_value("REPLACE_ETA_NAME", eta_name, template)
        replace_value("REPLACE_PT_NAME", pT_name, template)
    # eta_name = make_name_from_ls(eta_range, "eta")
    # pT_name = make_name_from_ls(pT_range, "pT")
        # Globals.
        # replace_value("REPLACE_OVERWRITE", overwrite, template)
        # replace_value("REPLACE_ITERS", iters, template)
        # replace_value("REPLACE_VERBOSE", verbose, template)
        # replace_value("REPLACE_fit_whole_range_first_iter", fit_whole_range_first_iter, template)
        # replace_value("REPLACE_use_data_in_xlim", use_data_in_xlim, template)
        # replace_value("REPLACE_binned_fit", binned_fit, template)
        # replace_value("REPLACE_switch_to_binned_fit", switch_to_binned_fit, template)

def print_info(fullpath_copy_main_script, fullpath_copy_slurm_script, 
                outdir_copies, outdir_pkl, outdir_txt, outdir_pdf):
    """Print the filepath locations."""
    print(f"[INFO] Main script copy created:\n{fullpath_copy_main_script}")
    print(f"[INFO] SLURM script copy created:\n{fullpath_copy_slurm_script}")
    print(f"[INFO] Dir to copies: {outdir_copies}")
    print(f"[INFO] Dir to pickle: {outdir_pkl}")
    print(f"[INFO] Dir to txt:    {outdir_txt}")
    print(f"[INFO] Dir to pdf:    {outdir_pdf}")

def main():
    eta_name = make_name_from_ls(eta_range, "eta")
    for pT_min, pT_max in zip(full_pT_ls[:-1], full_pT_ls[1:]):
        pT_range = [pT_min, pT_max]
        pT_name = make_name_from_ls(pT_range, "pT")
        print(f"...Preparing work area for: eta={eta_name}, pT={pT_name}")
        prep_area([outdir_copies, outdir_pkl, outdir_txt])
        fullpath_copy_main_script, fullpath_copy_slurm_script = make_filepaths_of_copies(eta_name, pT_name, template_script_main, template_script_slurm, outdir_copies)
        print_info(fullpath_copy_main_script, fullpath_copy_slurm_script, 
                   outdir_copies, outdir_pkl, outdir_txt, outdir_pdf)
        shutil.copyfile(template_script_main, fullpath_copy_main_script)
        shutil.copyfile(template_script_slurm, fullpath_copy_slurm_script)
        print(fullpath_copy_main_script)
        inpkl_path = inpkl_path_template.replace("ETAPART", eta_name).replace("PTPART", pT_name)
        outpkl_path = os.path.join(outdir_pkl, f"{job_name_base}_{eta_name}_{pT_name}.pkl")
        replace_vals_in_files(eta_range=eta_range, pT_range=pT_range,
                              inpkl_path=inpkl_path, outpkl_path=outpkl_path,
                              job_name_base=job_name_base, fullpath_new_main_script=fullpath_copy_main_script,
                              eta_name=eta_name, pT_name=pT_name,
                              template_tup=(fullpath_copy_main_script, fullpath_copy_slurm_script)
                              )
        print(f"...Submitting SLURM script for: eta_range={eta_range} pT_range={pT_range}")  #GOOD
        output = subprocess.run(["sbatch", fullpath_copy_slurm_script])  #GOOD
            
if __name__ == "__main__":
    main()
