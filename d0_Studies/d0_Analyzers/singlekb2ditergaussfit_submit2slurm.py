"""MyMuon Reco pT Corrector

The only way to evaluate m2mu_corr or m4mu_corr is by extracting muons from a
file. Thus far you can't open up a pickled MyMuonCollection and apply pT corr
to uncorrected m4mu vals.

Idea: class Event
"""
import os
import shutil
from glob import glob
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, replace_value, make_dirs
from Utils_Python.SlurmManager import SLURMSubmitter

verbose = 1
overwrite = 0
delete_kb2d_muon_ls = 1

inpkl_kb2d_path_glob = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p01_d0_1000p0/*eta*_*pT*0.pkl"
fullpath_main_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/singlekb2ditergaussfit_slurm_template.py"

new_filename_prefix = "unbinnedfit_widerwindow_fitwholerangefirstiter"
outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/beforeafterGeoFitCorr/0p01_d0_1000p0/copies/"
outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/beforeafterGeoFitCorr/0p01_d0_1000p0/output/"

switch_to_binned_fit = 99E9

# SLURM directives.
partition = "bigmem" #"hpg2-compute"
mem = (128, 'gb')
nodes = 4
burst = True
time = "02:00:00"  # hh:mm:ss

# # Perform iterated Gaussian fits on each KinBin2D.
outpkl_dir = os.path.dirname(inpkl_kb2d_path_glob)
kb2d_pkl_ls = glob(inpkl_kb2d_path_glob)
# template_basename = template.rstrip(".py")
for inpkl_path in kb2d_pkl_ls:
    # Create main copy.
    base = os.path.basename(fullpath_main_script).rstrip(".py")
    kb2d_bin_key = os.path.basename(inpkl_path).rstrip(".pkl")
    if len(new_filename_prefix) > 0:
        base = f"{new_filename_prefix}_{base}"
    file_name_copy = f"{base}_copy_{kb2d_bin_key}"
    template_copy = os.path.join(outcopies_dir, f"{file_name_copy}.py")
    make_dirs(outcopies_dir)
    make_dirs(outtxt_dir)
    shutil.copyfile(fullpath_main_script, template_copy)
    outpkl_suffix = f"{kb2d_bin_key}_withkb2dfits.pkl"
    outpkl_filename = f"{new_filename_prefix}_{outpkl_suffix}" if len(new_filename_prefix) > 0 else outpkl_suffix
    outpkl_path = os.path.join(outpkl_dir, outpkl_filename)
    # Modify main copy.
    replace_value("INPKL_PATH", inpkl_path, template_copy)
    replace_value("OUTPKL_PATH", outpkl_path, template_copy)
    replace_value("OVERWRITE", overwrite, template_copy)
    replace_value("DELETE_KB2D_MUON_LS", delete_kb2d_muon_ls, template_copy)
    replace_value("SWITCH2BINNED", switch_to_binned_fit, template_copy)
    
    sbmtr = SLURMSubmitter(verbose=verbose)
    sbmtr.prep_directives(
        job_name=file_name_copy,
        output_txt=os.path.join(outtxt_dir, f"{file_name_copy}.log"),
        email="rosedj1@ufl.edu",
        time=time,
        acct="avery",
        burst=burst,
        mem=mem,
        partition=partition,
        nodes=nodes,
    )
    cmdtup = (f"python {template_copy}")
    fullpath_slurm_copy = os.path.join(outcopies_dir, f"{file_name_copy}.sbatch")
    result = sbmtr.make_slurm_script(fullpath_slurm_copy, cmdtup, overwrite=overwrite)
    if result == 0:
        sbmtr.submit_script(fullpath_slurm_copy)

print(f"[INFO] Copies stored at:\n{outcopies_dir}")
print(f"[INFO] Output stored at:\n{outtxt_dir}")
print(f"[INFO] Pickle stored at:\n{outpkl_dir}")