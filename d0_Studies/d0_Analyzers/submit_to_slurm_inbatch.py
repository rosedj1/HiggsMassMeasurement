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
    The main template script should have the following replacement strings:
    - `REPLACE_ETA_LS` : wherever you want a 2-elem list of eta vals to go.
    - `REPLACE_JOB_NAME` : will be used as file name of pdf and output.txt files
    Example SLURM script should have: 
    - `REPLACE_JOB_NAME`
    - `REPLACE_ETA_NAME`
    - `REPLACE_NEW_FILE` : file path to main template script

Author: Jake Rosenzweig
Updated: 2021-02-05
"""
from Utils_Python.Utils_Files import replace_value, make_dirs
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
# from d0_Studies.d0_Utils.d0_cls import OrganizerKB2D
import subprocess
import shutil
import os
import sys

class SlurmManager:
    """A class to handle the details of working with SLURM."""

    def __init__(self): #, path_main, path_sbatch):
        self.file_copies_ls = []
        pass

    def copy_file(self, src, dst):
        """Copy `src` file to `dst`. Store the paths."""
        try:
            shutil.copyfile(src, dst)
            self.file_copies_ls.append(dst)
        except:
            print(f"[WARNING] Could not copy\n{src}to\n{dst}")

    def write_copies_to_file(self, outf):
        """Save the file paths of file_copies_ls to file `outf`."""
        # with open(outf) as f:
        #     f.
        pass

#-----------------------#
#--- User Parameters ---#
#-----------------------#
# job_name = "RC_vs_NoRC_itergaussfits_fullstats_pT75then200GeV_extendedxaxis"
job_name = "RC_vs_NoRC_itergaussfits_fullstats_pT75then200GeV_extendedxaxis"

template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/RochCorrAnalyzers/roch_vs_noroch_itergausfit_template.py"
# template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py"
template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/RochCorrAnalyzers/roch_vs_noroch_slurm.sbatch"

outdir_copies = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/output/copies/"
outdir_pkl    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/"
outdir_txt    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/output/"
outdir_pdf    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/plots/"

full_eta_ls = equal_entry_bin_edges_eta_mod1_wholenum

delete_copies = False
#------------------------#
#--- Script functions ---#
#------------------------#
def make_eta_name(eta_ls):
    """Parse a 2-element list of eta values into a title string."""
    eta_str = f"{eta_ls[0]}eta{eta_ls[1]}"
    eta_str = eta_str.replace(".", "p")
    return eta_str

def prep_area(dir_tup):
    """Make sure all directories exist. If not, create them."""
    for d in dir_tup:
        make_dirs(d)

def make_filepaths_of_copies(eta_ls, template_main, template_slurm, outdir_copies):
    """Return the full paths of copies of the main script and the slurm script.
    
    Parameters
    ----------
    outdir_copies : str
        Path to directory where copies of scripts will be stored.
    """
    eta_name = make_eta_name(eta_ls)
    suffix = f"_{eta_name}.py"
    filename_main = os.path.split(template_main)[1]
    filename_slurm = os.path.split(template_slurm)[1]
    fullpath_copy_main = os.path.join(outdir_copies, filename_main.replace(".py", f"_copy_{eta_name}.py"))
    fullpath_copy_slurm = os.path.join(outdir_copies, filename_slurm.replace(".sbatch", f"_copy_{eta_name}.sbatch"))
    return (fullpath_copy_main, fullpath_copy_slurm)

def replace_vals_in_files(eta_ls, job_name, fullpath_new_main_script,
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
    eta_name = make_eta_name(eta_ls)
    # Replace phrases in slurm script.
    for template in template_tup:
        replace_value("REPLACE_ETA_LS", eta_ls, template)
        replace_value("REPLACE_ETA_NAME", eta_name, template)
        replace_value("REPLACE_JOB_NAME", job_name, template)
        replace_value("REPLACE_NEW_SCRIPT", fullpath_new_main_script, template)
        replace_value("REPLACE_OUTDIR_PKL", outdir_pkl.rstrip('/'), template)
        replace_value("REPLACE_OUTDIR_PDF", outdir_pdf.rstrip('/'), template)
        replace_value("REPLACE_OUTDIR_TXT", outdir_txt.rstrip('/'), template)

if __name__ == "__main__":
    for eta_min, eta_max in zip(full_eta_ls[:-1], full_eta_ls[1:]):
        print(f"...Preparing work area.")
        prep_area([outdir_copies, outdir_pkl, outdir_txt, outdir_pdf])
        eta_ls = [eta_min, eta_max]
        print(f"[INFO] Making copies of scripts...")
        fullpath_copy_main_script, fullpath_copy_slurm_script = make_filepaths_of_copies(eta_ls, template_script_main, template_script_slurm, outdir_copies)
        shutil.copyfile(template_script_main, fullpath_copy_main_script)
        shutil.copyfile(template_script_slurm, fullpath_copy_slurm_script)
        replace_vals_in_files(eta_ls, job_name, fullpath_copy_main_script,
                          outdir_pkl, outdir_pdf, outdir_txt,
                          template_tup=(fullpath_copy_main_script, fullpath_copy_slurm_script))
        print(f"...Submitting SLURM script for eta range: {eta_ls}")
        output = subprocess.run(["sbatch", fullpath_copy_slurm_script])
        if delete_copies:
            print(f"[INFO] Removing copies of scripts.")
            os.remove(new_main)
            os.remove(new_slurm)
