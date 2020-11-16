"""Muon pT Corrector

Select m4mu events from a ggH sample which pass selections.
Apply pT correction factors per MyMuon, from a pT corr dict.
"""
import pickle
import ROOT as r
from Utils_ROOT.Printer import CanvasPrinter
from ParticleCollections import MyMuonCollection
from Utils_Python.Utils_Files import check_overwrite
#---------------------------#
#----- User Parameters -----#
#---------------------------#
# suffix = "slurm_test02"
suffix = "test25"
infileroot_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2018.root"
# inpkl_pT_corr_factor_dict = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/BestSoFar/MC2018ggH_combined_pT_corr_factor_dct.pkl"
inpkl_pT_corr_factor_dict = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/small_test_pT_corr_factor_dict.pkl"

outpath_rootfile = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorr_vals_{suffix}.root"
outpath_pdf      = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Plots/applypTcorrplots/MC2018_dpTOverpTdist_beforeafterpTcorr_{suffix}.pdf"
outpath_pkl      = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/AfterApplyCorr/afterapplycorr_kb2d_dict_{suffix}.pkl"

verbose = False
overwrite = True
do_mu_pT_corr = True
only_draw_last = True
max_n_evts = 10000
print_out_every = 5000
iters = 3
switch_to_binned_fit = 10000
make_m4mu_rootfile = False
#----------------------------#
#----- Script Functions -----#
#----------------------------#
if __name__ == "__main__":
    for f in [outpath_rootfile, outpath_pdf, outpath_pkl]:
        check_overwrite(f, overwrite)
    printer = CanvasPrinter(show_plots=False, gridOn=True, canv=None)
    # Get pT corr factors.
    with open(inpkl_pT_corr_factor_dict, "rb") as p:
        pT_corr_factor_dict = pickle.load(p)
    # Make muon container.
    muon_collection = MyMuonCollection(["ggH"])
    # Extract muons and apply pT corrections.
    muon_collection.extract_muons_from_Higgs_file(infileroot_path,
                                                  n_evts=max_n_evts,
                                                  print_out_every=print_out_every,
                                                  do_mu_pT_corr=do_mu_pT_corr, 
                                                  pT_corr_factor_dict=pT_corr_factor_dict,
                                                  verbose=verbose)
    # Save m4mu and m4mu_corr info in separate root file for DSCB fits.
    if make_m4mu_rootfile:
        muon_collection.write_m4muinfo_to_rootfile(outpath_rootfile, overwrite=overwrite)
    # Make inclusive kinematic plots. Stores MuonCollection.hist_inclusive_ls
    muon_collection.make_inclusive_kinematic_plots()
    # Put all muons in respective KinBin2Ds.
    # Then make and fill all KinBin2D hists.
    muon_collection.create_KinBins_from_pT_corr_dict(pT_corr_factor_dict)
    muon_collection.sort_muons(eta_ls=None, pT_ls=None,
                               pT_corr_factor_dict=pT_corr_factor_dict,
                               n_bins_dpTOverpT=100, x_lim_dpTOverpT=[-0.4, 0.4],
                               n_bins_qd0=100, x_lim_qd0=[-0.01, 0.01], 
                               verbose=verbose)
    # KinBin2Ds now have muons with before pT corr and after.
    # Perform iterated Gaussian fits on each KinBin2D.
    muon_collection.do_2D_iter_gaus_fits(
                                bins_dpTOverpT=100, bins_qd0=100,
                                x_lim_dpTOverpT=[-0.4,0.4], x_lim_qd0=[-0.01,0.01],
                                fit_whole_range_first_iter=True,
                                iters=iters, num_sigmas=2.5,
                                switch_to_binned_fit=switch_to_binned_fit,
                                verbose=verbose, alarm_level="warning",
                                use_mu_pT_corr=do_mu_pT_corr, only_draw_last=only_draw_last)
    # muon_collection.make_plots_beforeafterpTcorr()
    # plot_ls = muon_collection.hist_inclusive_ls
    # canv_ls = []
    c = r.TCanvas()
    c.Print(outpath_pdf + "[")
    for kb2d in muon_collection.KinBin2D_dict.values():
        kb2d.overwrite_muon_info(delete_all=False)
        kb2d.make_beforeafterpTcorr_frames()
        kb2d.draw_beforeafterpTcorr()
        c.Print(outpath_pdf)
        # canv_ls.append(kb2d.canv_beforeafterpTcorr)
        # plot_ls.append(kb2d.frame_dpTOverpT)
        # plot_ls.append(kb2d.frame_dpTOverpT_corr)
        # c.Clear()
    c.Print(outpath_pdf + "]")
    # printer.make_pdf_from_canv(canv_ls, outpath_pdf)
    # printer.make_pdf_of_plots(plot_ls, outpath_pdf)
    muon_collection.save_to_pkl(muon_collection.KinBin2D_dict, outpath_pkl)