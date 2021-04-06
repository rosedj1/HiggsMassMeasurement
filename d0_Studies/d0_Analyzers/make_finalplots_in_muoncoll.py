# from ROOT import TCanvas
import os
from Utils_Python.Utils_Files import open_pkl, make_dirs, check_overwrite
from Utils_ROOT.Printer import CanvasPrinter

# outtxt_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/output/"
#--- |d0| < 0.01 cm:
# mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p0_d0_0p01.pkl")
# mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p0_d0_0p01_1000bins_neg0p2xlim0p2.pkl")
mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p0_d0_0p01_unbinned.pkl")
#--- |d0| > 0.01 cm:
# mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p01_d0_1000p0_unbinned.pkl")
# mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p01_d0_1000p0_unbinned_widerwindow_fitwholerangefirstiter.pkl")

# outpdf_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/kb2ditergaussfits_beforeafterGeoFitcorr/unbinnedIGFs_widerwindow_fitwholerangefirstiter__0p01_d0_1000p0.pdf"
outpdf_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/kb2ditergaussfits_beforeafterGeoFitcorr/unbinnedIGFs_widerwindow_fitwholerangefirstiter__0p0_d0_0p01.pdf"
overwrite = 0
make_pdf = 0  # If no, then just print values.

# kb2d = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_nogenmatching_0p0_d0_0p01/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_nogenmatching_0p0_d0_0p01_0p0eta2p4_10p0pT15p0.pkl")
make_dirs(os.path.dirname(outpdf_path))

# printer.make_pdf_of_plots()
# kb3d.frame_dpTOverpT.Draw()

# Let's make some plots.

# mu_coll.make_plots_beforeafterpTcorr()
# plot_ls = mu_coll.hist_inclusive_ls
# c = TCanvas()
if make_pdf:
    check_overwrite(outpdf_path, overwrite=overwrite)
    printer = CanvasPrinter(show_plots=0)
    printer.make_plots_pretty()
    c = printer.canv
    c.Print(outpdf_path + "[")
# with open(outtxt_path, "w") as f:
print(
    f"{'eta_min':<10}, {'eta_max':<10}, "
    f"{'pT_min':<10}, {'pT_max':<10}, "
    f"{'bf mu before_corr (GeV)':<25}, {'bf mu_err before_corr (GeV)':<25}, "
    f"{'bf mu after_corr (GeV)':<25}, {'bf mu_err after_corr (GeV)':<25}, "
    f"{'bf sig before_corr (GeV)':<25}, {'bf sig_err before_corr (GeV)':<25}, "
    f"{'bf sig after_corr (GeV)':<25}, {'bf sig_err after_corr (GeV)':<25}, "
    f"{'sigma_perc_improve (%)':<25}, {'sigma_perc_improve_err (%)':<25}, "
    f"{'chi2_beforecorr_ls[-1]':<25}, {'chi2_aftercorr_ls[-1]':<25}, "
    f"{'n_muons':<10}"
    )
for ct, kb2d in enumerate(mu_coll.KinBin2D_dict.values()):
    # kb2d.overwrite_muon_info(delete_all=True)
    if make_pdf:
        kb2d.make_beforeafterpTcorr_frames()
        kb2d.draw_beforeafterpTcorr()
        c.Print(outpdf_path)
    # plot_ls.append(kb2d.frame_dpTOverpT)
    # plot_ls.append(kb2d.frame_dpTOverpT_corr)
    print(
        f"{kb2d.eta_min:<10}, {kb2d.eta_max:<10}, "
        f"{kb2d.pT_min:<10}, {kb2d.pT_max:<10}, "
        f"{kb2d.fit_stats_dict_dpTOverpT['mean_ls'][-1]:<25}, {kb2d.fit_stats_dict_dpTOverpT['mean_err_ls'][-1]:<25}, "
        f"{kb2d.fit_stats_dict_dpTOverpT_corr['mean_ls'][-1]:<25}, {kb2d.fit_stats_dict_dpTOverpT_corr['mean_err_ls'][-1]:<25}, "
        f"{kb2d.fit_stats_dict_dpTOverpT['std_ls'][-1]:<25}, {kb2d.fit_stats_dict_dpTOverpT['std_err_ls'][-1]:<25}, "
        f"{kb2d.fit_stats_dict_dpTOverpT_corr['std_ls'][-1]:<25}, {kb2d.fit_stats_dict_dpTOverpT_corr['std_err_ls'][-1]:<25}, "
        f"{kb2d.sigma_perc_improve:<25}, {kb2d.sigma_perc_improve_err:<25}, "
        f"{kb2d.fit_stats_dict_dpTOverpT['chi2_ls'][-1]:<25}, {kb2d.fit_stats_dict_dpTOverpT_corr['chi2_ls'][-1]:<25}, "
        f"{kb2d.n_entries:<10}"
        )
if make_pdf:
    c.Print(outpdf_path + "]")
# printer.make_pdf_of_plots(plot_ls, outpdf_path)
# save_to_pkl(mu_coll, outpkl_path, overwrite=overwrite)