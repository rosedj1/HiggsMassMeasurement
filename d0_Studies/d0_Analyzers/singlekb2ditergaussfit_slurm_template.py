"""MyMuon Reco pT Corrector

Apply pT corrections stored in a '.pkl' to a list of individual pickled KB2Ds.
The muons will the store their pT_reco and pT_corr values.

NOTE: This script doesn't yet produce m4mu_corr.
The only way to evaluate m2mu_corr or m4mu_corr is by extracting muons from a
file. Thus far you can't open up a pickled MyMuonCollection and apply pT corr
to uncorrected m4mu vals.
Idea: class Event, which can store all 4 muons per event (hence m4mu).

This script is controlled by:
d0_Studies/d0_Analyzers/singlekb2ditergaussfit_submit2slurm.py

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: <2021-04-12
Updated: 2021-04-14
"""
from Utils_Python.Utils_Files import open_pkl, save_to_pkl
from d0_Studies.d0_Utils.d0_fns import correct_muon_pT

overwrite = OVERWRITE
verbose = VERBOSE
inpkl_path_kb2d = "INPKL_PATH"
inpkl_path_pTcorrfactordict = "INPKL_PTCORRFACTORDICT"
outpkl_path_kb2d = "OUTPKL_PATH"
delete_kb2d_muon_ls = DELETE_KB2D_MUON_LS
switch_to_binned_fit = SWITCH2BINNED
use_smart_window = USE_SMART_WINDOW

# eta_binedge_ls = ETA_BINEDGE_LS
# pT_binedge_ls = PT_BINEDGE_LS

bins_dpTOverpT = BINS_DPTOVERPT
bins_qd0 = BINS_QD0
x_lim_dpTOverpT = X_LIM_DPTOVERPT
x_lim_qd0 = X_LIM_QD0
iters = ITERS

kb2d = open_pkl(inpkl_path_kb2d)
pT_corr_factor_dct = open_pkl(inpkl_path_pTcorrfactordict)
# Correct muon pTs.
for mu in kb2d.muon_ls:
    mu.pT_corr = correct_muon_pT(mu.eta, mu.pT, mu.charge, mu.d0, 
                    pT_corr_factor_dct, 
                    # eta_binedge_ls=eta_binedge_ls,
                    # pT_binedge_ls=pT_binedge_ls,
                    detection="auto",
                    verbose=verbose)

# Without pT corr. Also fit the q*d0 dist.
kb2d.do_itergausfit(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0, 
                    x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0, 
                    fit_whole_range_first_iter=True, 
                    iters=iters, num_sigmas=2.5, marker_color=None, line_color=None, 
                    switch_to_binned_fit=switch_to_binned_fit, verbose=verbose, alarm_level="warning", 
                    use_mu_pT_corr=False, only_draw_last=False, fit_qd0_dist=True, use_smart_window=use_smart_window) 
# WITH pT corr.
kb2d.do_itergausfit(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0, 
                    x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0, 
                    fit_whole_range_first_iter=True, 
                    iters=iters, num_sigmas=2.5, marker_color=None, line_color=None, 
                    switch_to_binned_fit=switch_to_binned_fit, verbose=verbose, alarm_level="warning", 
                    use_mu_pT_corr=True, only_draw_last=False, use_smart_window=use_smart_window)
kb2d.overwrite_muon_info(delete_all=True)
save_to_pkl(kb2d, outpkl_path_kb2d, overwrite=overwrite)