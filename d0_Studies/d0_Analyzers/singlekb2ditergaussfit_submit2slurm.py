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

inpkl_kb2d_path_glob = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/*eta*_*pT*0.pkl"
# inpkl_kb2d_path_glob = [
#     "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/0p0eta2p4_5p0pT10p0.pkl",
#     "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/0p0eta2p4_10p0pT15p0.pkl",
#     "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/0p0eta2p4_15p0pT20p0.pkl",
#     "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/0p0eta2p4_50p0pT60p0.pkl",
#     "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/0p0eta2p4_60p0pT100p0.pkl",
#     "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/0p0eta2p4_100p0pT200p0.pkl"
# ]
fullpath_main_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/singlekb2ditergaussfit_slurm_template.py"

new_filename_prefix = "unbinnedfit"
outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/beforeafterGeoFitCorr/copies/"
outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/beforeafterGeoFitCorr/output/"

switch_to_binned_fit = 99E9

# SLURM directives.
partition = "bigmem" #"hpg2-compute"
mem = (128, 'gb')
nodes = 4
burst = True
time = "02:00:00"  # hh:mm:ss

# inpkl_path_kb2d          = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p0_d0_0p01.pkl"
# inpkl_path_pTcorrfactors = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/GeoFit/GeoFitTcorrfact_derivedfromMC2016_3etabins_0p0eta2p4_newerformat.pkl"
# outdir_kb2d_dicts = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/"

# outdir_kb2d_dicts
# GOOD muon_coll           = open_pkl(inpkl_path_kb2d)
# GOOD pT_corr_factor_dict = open_pkl(inpkl_path_pTcorrfactors)

# GOOD muon_coll.apply_pTcorr_to_all_muons(pT_corr_factor_dict, use_GeoFit_algo=1, force_zero_intercept=True, verbose=True)
# # Analyze each KB2D one at a time on SLURM.
# GOOD muon_coll.save_KB2Ds_separate_dcts(outdir=outdir_kb2d_dicts, file_prefix="", overwrite=overwrite, verbose=verbose)

# # Perform iterated Gaussian fits on each KinBin2D.
new_outpkl_dir = os.path.dirname(inpkl_kb2d_path_glob)
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
    outpkl_path = os.path.join(new_outpkl_dir, outpkl_filename)
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
# .do_itergausfit(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0,
#                     x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0,
#                     fit_whole_range_first_iter=fit_whole_range_first_iter,
#                     iters=iters, num_sigmas=num_sigmas,
#                     switch_to_binned_fit=switch_to_binned_fit,
#                     verbose=verbose, alarm_level=alarm_level,
#                     use_mu_pT_corr=use_mu_pT_corr, only_draw_last=only_draw_last)