"""SLURM Script Duplicator and Submitter

This code submits SLURM jobs, where each SLURM job performs iterated Gaussian
fits on a single KB2D. Fits are done on dpT/pT and q*d0 distributions. Fit
statistics are saved. The updated KB2D is saved as a new '.pkl'.

Iterated fits can hang for a very long time (memory issues?). A workaround is
what this script does: only work on 1 KB2D at a time.

The kb2d.muon_ls and kb3d.muon_ls are erased to save memory.

TODO:
- Implement SlurmManager class.

This code makes a copy of a "main template" script that you wish to run on
SLURM. It actually makes many copies of the main script and of the SLURM
submissions script, each copy only differing from the next by the eta bin
used.

User can put in a whole list of eta values, and each bin will be processed:

Example: eta_ls = [0.2, 0.4, 0.6, 0.8]
This code will produce a copy of the main script
and of the SLURM submission script using eta_bin = [0.2, 0.4].
This job will be submitted to SLURM.
Then the next copy will use eta_bin = [0.4, 0.6], etc.

NOTE:
    The main template script should have replacement strings that are in all
    caps: "REPLACE_". Check the `replace_vals_in_files` function to see what values
    get replaced.

Requires an input pickled dict of KinBin2Ds. You can make the '.pkl' with:
- save_kb2ds_separate_dicts.py
- roch_vs_noroch_kb2dmaker_inclusivehistplotter.py

Author: Jake Rosenzweig
Created: <2021-03-15
Updated: 2021-04-14
"""
from Utils_Python.Utils_Files import replace_value, make_dirs
from d0_Studies.kinematic_bins import eta_bins_geofitvsVXBS, pT_bins_geofitvsVXBS, equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
from Utils_Python.SlurmManager import SLURMSubmitter
# from d0_Studies.d0_Utils.d0_cls import OrganizerKB2D
import subprocess
import shutil
import os
import sys
#-----------------------#
#--- User Parameters ---#
#-----------------------#
year = "2016"
overwrite = 1
verbose = 1

iters = 5
fit_with_zero_interc = True
regions = 12
min_muons_per_qd0_bin = 100
fit_whole_range_first_iter = False  # False gives more consistent fits (with no outlier data).
use_smart_window = True
use_data_in_xlim = 1
binned_fit = False
switch_to_binned_fit = 999999999

job_name_base = "individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc_redo"  # Also prefix for outfile.
# job_name_base = "MC2016DY_individKB2D_withitergaussfitsonKB3Ds_nogenmatching_0p0_d0_0p01"  # Also prefix for outfile.
# eta_ls = eta_bins_geofitvsVXBS #equal_entry_bin_edges_eta_mod1_wholenum
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum
# full_pT_ls = pT_bins_geofitvsVXBS  #[20.0, 30.0, 40.0, 50.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
full_pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum

template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_singlekb2d_itergaussfitsonkb3ds_template.py"
# template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py"
template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_to_slurm_template.sbatch"
# template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/RochCorrAnalyzers/roch_vs_noroch_slurm.sbatch"

# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dictsnofitinfo/MC2016DY_fullstats_muoncoll_withkb3dbins_nogenmatching_0p0_d0_0p01__ETAPART_PTPART.pkl"
inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/skim_fullstats_verify/pickles/kb2d_dicts/ETAPART_PTPART.pkl"

# /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/output/DeriveCorrections/2016DY/individualKB2Ditergaussfits/MC2016DY_kb2d_ETARANGE_PTRANGE_error.log
outdir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu"

# SLURM directives.
partition = "hpg2-compute"
mem = (8, 'gb')
nodes = 1
burst = True
time = "01:00:00"  # hh:mm:ss

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
                          job_name_base, outdir_txt, fullpath_new_main_script,
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
        replace_value("VERBOSE", verbose, template)
        replace_value("OVERWRITE", overwrite, template)
        replace_value("REPLACE_ETA_LS", eta_range, template)
        replace_value("REPLACE_PT_RANGE", pT_range, template)
        replace_value("REPLACE_INPKL_PATH", inpkl_path, template)
        replace_value("REPLACE_OUTPATH_PKL", outpkl_path, template)
        replace_value("REPLACE_JOB_NAME", job_name_base, template)
        replace_value("REPLACE_OUTDIR_TXT", outdir_txt.rstrip('/'), template)
        replace_value("REPLACE_NEW_SCRIPT", fullpath_new_main_script, template)
        replace_value("REPLACE_ETA_NAME", eta_name, template)
        replace_value("REPLACE_PT_NAME", pT_name, template)
        replace_value("BINNED_FIT", binned_fit, template)
        replace_value("ITERS", iters, template)
        replace_value("FIT_WITH_ZERO_INTERC", fit_with_zero_interc, template)
        replace_value("REGIONS", regions, template)
        replace_value("MIN_MUONS_PER_QD0_BIN", min_muons_per_qd0_bin, template)
        replace_value("USE_SMART_WINDOW", use_smart_window, template)

def print_info(fullpath_copy_main_script, fullpath_copy_slurm_script, 
                outdir_copies, outdir_pkl, outdir_txt):
    """Print the filepath locations."""
    print(f"[INFO] Main script copy created:\n{fullpath_copy_main_script}")
    print(f"[INFO] SLURM script copy created:\n{fullpath_copy_slurm_script}")
    print(f"[INFO] Dir to copies: {outdir_copies}")
    print(f"[INFO] Dir to pickle: {outdir_pkl}")
    print(f"[INFO] Dir to outtxt: {outdir_txt}")

def main():
    outdir_copies = os.path.join(outdir, f"{job_name_base}/copies")
    outdir_pkl    = os.path.join(outdir, f"{job_name_base}/pickles")
    outdir_txt    = os.path.join(outdir, f"{job_name_base}/output")
    assert all(year in f for f in (inpkl_path_template, outdir_copies, outdir_pkl, outdir_txt))
    prep_area([outdir_copies, outdir_pkl, outdir_txt])
    for eta_min, eta_max in zip(eta_ls[:-1], eta_ls[1:]):
        eta_range = [eta_min, eta_max]
        eta_name = make_name_from_ls(eta_range, "eta")
        for pT_min, pT_max in zip(full_pT_ls[:-1], full_pT_ls[1:]):
            pT_range = [pT_min, pT_max]
            pT_name = make_name_from_ls(pT_range, "pT")
            print(f"...Preparing work area for: eta={eta_name}, pT={pT_name}")
            fullpath_copy_main_script, fullpath_copy_slurm_script = make_filepaths_of_copies(eta_name, pT_name, template_script_main, template_script_slurm, outdir_copies)
            print_info(fullpath_copy_main_script, fullpath_copy_slurm_script, 
                    outdir_copies, outdir_pkl, outdir_txt)
            shutil.copyfile(template_script_main, fullpath_copy_main_script)
            shutil.copyfile(template_script_slurm, fullpath_copy_slurm_script)
            print(fullpath_copy_main_script)
            inpkl_path = inpkl_path_template.replace("ETAPART", eta_name).replace("PTPART", pT_name)
            outpkl_path = os.path.join(outdir_pkl, f"{eta_name}_{pT_name}.pkl")
            replace_vals_in_files(eta_range=eta_range, pT_range=pT_range,
                                inpkl_path=inpkl_path, outpkl_path=outpkl_path,
                                job_name_base=job_name_base, outdir_txt=outdir_txt,
                                fullpath_new_main_script=fullpath_copy_main_script,
                                eta_name=eta_name, pT_name=pT_name,
                                template_tup=(fullpath_copy_main_script, fullpath_copy_slurm_script)
                                )
            print(f"...Submitting SLURM script for: eta_range={eta_range} pT_range={pT_range}")

            # Prep and submit SLURM script copy.
            # sbmtr = SLURMSubmitter(verbose=False)
            # sbmtr.prep_directives(
            #     job_name=full_file_name,
            #     REPLACE_OUTDIR_TXT/REPLACE_JOB_NAME_REPLACE_ETA_NAME_REPLACE_PT_NAME.log = os.path.join(new_outtxt_dir, f"{full_file_name}.log"
            #     output_txt=os.path.join(new_outtxt_dir, f"{full_file_name}.log"),
            #     email="rosedj1@ufl.edu",
            #     time=time,
            #     acct="avery",
            #     burst=burst,
            #     mem=mem,
            #     partition=partition,
            #     nodes=nodes,
            # )
            # cmdtup = (f"time python {fullpath_main_script_copy}")
            # fullpath_slurm_copy = os.path.join(new_outcopies_dir, f"{full_file_name}.sbatch")
            # result = sbmtr.make_slurm_script(fullpath_slurm_copy, cmdtup, overwrite=overwrite)
            # if result == 0:
            #     sbmtr.submit_script(fullpath_slurm_copy)

            output = subprocess.run(["sbatch", fullpath_copy_slurm_script])
            
if __name__ == "__main__":
    main()
