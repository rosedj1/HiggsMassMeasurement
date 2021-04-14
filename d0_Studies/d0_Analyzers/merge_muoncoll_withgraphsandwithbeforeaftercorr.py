"""Post-Fit MyMuonCollection Merger

Merges two MyMuonCollections with the same KB2Ds together.
mucoll_1 : owns KB2Ds with before/after pT correction fits.
mucoll_2 : owns KB2Ds used to derive pT corrections.

Merge mucoll_2 into mucoll_1.

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-04-13
Updated: 
"""
# import os
# from pprint import pprint
# from glob import glob
# from Utils_ROOT.Printer import CanvasPrinter
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, make_dirs
# from d0_Studies.d0_Utils.d0_fns import get_pT_part, parse_etapT_key
# from Utils_Python.printing import print_header_message
# from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
# from ParticleCollections import MyMuonCollection
# from natsort import natsorted

overwrite = 0
# This MyMuonCollection owns KB2Ds used to derive pT corrections.
inpkl_mucoll_derivecorr = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/muoncoll_itergaussfitsonKB3Ds_redo.pkl"
# This MyMuonCollection owns KB2Ds with before/after pT correction fits.
inpkl_mucoll_applycorr = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/muoncoll_withKB2Dfits_beforeaftercorr.pkl"

mucoll_derivecorr = open_pkl(inpkl_mucoll_derivecorr)
mucoll_applycorr = open_pkl(inpkl_mucoll_applycorr)

def merge_attr(kb2d_keep, kb2d_toss):
    """Set attributes of kb2d_keep to those of kb2d_toss."""
    kb2d_keep.iters = kb2d_toss.iters
    kb2d_keep.use_mu_pT_corr = kb2d_toss.use_mu_pT_corr
    kb2d_keep.fit_stats_dict_dpTOverpT = kb2d_toss.fit_stats_dict_dpTOverpT
    kb2d_keep.frame_dpTOverpT = kb2d_toss.frame_dpTOverpT
    kb2d_keep.n_entries = kb2d_toss.n_entries
    kb2d_keep.fit_stats_dict_qd0 = kb2d_toss.fit_stats_dict_qd0
    kb2d_keep.frame_qd0 = kb2d_toss.frame_qd0
    kb2d_keep.bestfit_mean_beforecorr = kb2d_toss.bestfit_mean_beforecorr
    kb2d_keep.bestfit_meanerr_beforecorr = kb2d_toss.bestfit_meanerr_beforecorr
    kb2d_keep.bestfit_std_beforecorr = kb2d_toss.bestfit_std_beforecorr
    kb2d_keep.bestfit_stderr_beforecorr = kb2d_toss.bestfit_stderr_beforecorr
    kb2d_keep.fit_stats_dict_dpTOverpT_corr = kb2d_toss.fit_stats_dict_dpTOverpT_corr
    kb2d_keep.frame_dpTOverpT_corr = kb2d_toss.frame_dpTOverpT_corr
    kb2d_keep.bestfit_mean_aftercorr = kb2d_toss.bestfit_mean_aftercorr
    kb2d_keep.bestfit_meanerr_aftercorr = kb2d_toss.bestfit_meanerr_aftercorr
    kb2d_keep.bestfit_std_aftercorr = kb2d_toss.bestfit_std_aftercorr
    kb2d_keep.bestfit_stderr_aftercorr = kb2d_toss.bestfit_stderr_aftercorr
    kb2d_keep.sigma_perc_improve = kb2d_toss.sigma_perc_improve
    kb2d_keep.sigma_perc_improve_err = kb2d_toss.sigma_perc_improve_err
    return (kb2d_keep, kb2d_toss)

def main():
    outfile = inpkl_mucoll_applycorr.replace(".pkl", "_merged.pkl")
    check_overwrite(outfile, overwrite=overwrite)
    print("...Merging before/after MuonColl into Derive MuonColl...")
    for key in mucoll_derivecorr.KinBin2D_dict.keys():
        msg = f"key=`{key}` not found in\n{inpkl_mucoll_applycorr}"
        assert key in mucoll_applycorr.KinBin2D_dict.keys(), msg
        kb2d_keep = mucoll_derivecorr.KinBin2D_dict[key]
        kb2d_toss = mucoll_applycorr.KinBin2D_dict[key]
        kb2d_keep, kb2d_toss = merge_attr(kb2d_keep, kb2d_toss)
        mucoll_applycorr.KinBin2D_dict[key] = kb2d_keep
    save_to_pkl(mucoll_derivecorr, outfile)

if __name__ == "__main__":
    main()