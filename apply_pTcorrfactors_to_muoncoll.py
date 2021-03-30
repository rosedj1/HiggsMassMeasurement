"""MyMuon Reco pT Corrector

The only way to evaluate m2mu_corr or m4mu_corr is by extracting muons from a
file. Thus far you can't open up a pickled MyMuonCollection and apply pT corr
to uncorrected m4mu vals.

Idea: class Event
"""
from Utils_Python.Utils_Files import open_pkl, save_to_pkl
from Utils_Python.SlurmManager import SLURMSubmitter

verbose = 1
overwrite = 0

inpkl_path_muoncoll      = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p0_d0_0p01.pkl"
inpkl_path_pTcorrfactors = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/GeoFit/GeoFitTcorrfact_derivedfromMC2016_3etabins_0p0eta2p4_newerformat.pkl"

burst = 1
mem = (8, 'gb')
partition = "hpg2-compute"
nodes = 4

full_file_name
new_outtxt_dir

muon_coll           = open_pkl(inpkl_path_muoncoll)
pT_corr_factor_dict = open_pkl(inpkl_path_pTcorrfactors)
muon_coll.apply_pTcorr_to_all_muons(pT_corr_factor_dict, use_GeoFit_algo=1, force_zero_intercept=True, verbose=True)
# KinBin2Ds now have muons with before pT corr and after.
# Perform iterated Gaussian fits on each KinBin2D.
muon_coll.split_kb2ds_and_save()

sbmtr = SLURMSubmitter(verbose=verbose)
sbmtr.prep_directives(
    job_name=full_file_name,
    output_txt=os.path.join(new_outtxt_dir, f"{full_file_name}.log"),
    email="rosedj1@ufl.edu",
    time="02:00:00",
    acct="avery",
    burst=burst,
    mem=mem,
    partition=partition,
    nodes=nodes,
)
cmdtup = (f"python {fullpath_main_script_copy}")
fullpath_slurm_copy = os.path.join(new_outcopies_dir, f"{full_file_name}.sbatch")
result = sbmtr.make_slurm_script(fullpath_slurm_copy, cmdtup, overwrite=overwrite)
if result == 0:
    sbmtr.submit_script(fullpath_slurm_copy)
# Get ready for next batch.
evt_beg = int(evt_end + 1)

kb2d_on_slurm.do_itergausfit(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0,
                    x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0,
                    fit_whole_range_first_iter=fit_whole_range_first_iter,
                    iters=iters, num_sigmas=num_sigmas,
                    switch_to_binned_fit=switch_to_binned_fit,
                    verbose=verbose, alarm_level=alarm_level,
                    use_mu_pT_corr=use_mu_pT_corr, only_draw_last=only_draw_last)

# muon_coll.do_2D_iter_gaus_fits(
#                             bins_dpTOverpT=100, bins_qd0=100,
#                             x_lim_dpTOverpT=[-0.4,0.4], x_lim_qd0=[-0.01,0.01],
#                             fit_whole_range_first_iter=True,
#                             iters=iters, num_sigmas=2.5,
#                             switch_to_binned_fit=switch_to_binned_fit,
#                             verbose=verbose, alarm_level="warning",
#                             use_mu_pT_corr=do_mu_pT_corr, only_draw_last=only_draw_last)

muon_coll.make_plots_beforeafterpTcorr()
# plot_ls = muon_coll.hist_inclusive_ls
c = ROOT.TCanvas()
c.Print(outpath_pdf + "[")
for kb2d in muon_coll.KinBin2D_dict.values():
    kb2d.overwrite_muon_info(delete_all=True)
    kb2d.make_beforeafterpTcorr_frames()
    kb2d.draw_beforeafterpTcorr()
    c.Print(outpath_pdf)
    # plot_ls.append(kb2d.frame_dpTOverpT)
    # plot_ls.append(kb2d.frame_dpTOverpT_corr)
c.Print(outpath_pdf + "]")
# printer.make_pdf_of_plots(plot_ls, outpath_pdf)
save_to_pkl(muon_coll, outpkl_path)