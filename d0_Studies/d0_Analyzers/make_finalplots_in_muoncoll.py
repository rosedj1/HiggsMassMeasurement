"""Make a PDF of all iterated Gaussian fits from KB3Ds.

Author: Jake Rosenzweig
Created: <2021-04-12
Updated: 2021-04-14
"""
import os
import ROOT
from Utils_Python.Utils_Files import open_pkl, check_overwrite, make_dirs
from Utils_ROOT.Printer import CanvasPrinter
from d0_Studies.d0_Utils.d0_dicts import color_dict_RooFit
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum

overwrite = 1
# scale_by_1divpT = 0
# scale_by_avgOf1divpT = 0
# scale_by_muOf1divpT = 0
print_vals = 0  # Eventually write vals to file.

year = "2016"
# filename = "MC2018DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds.pdf"
filename = "MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds_redo.pdf"
# infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY/pickles/final_muon_coll/MC2018DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds_new.pkl"
# infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc/pickles/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc_2p1eta2p2_5p0pT7p0.pkl"
infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/muoncoll_withKB2Dfits_beforeaftercorr_merged.pkl"
outpdf_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/ApplyCorr/MC2016DY2mu/"

eta_ls = equal_entry_bin_edges_eta_mod1_wholenum

# if __name__ == "__main__":
#     pass
outpdf_path = os.path.join(outpdf_dir, filename if '.pdf' in filename else f"{filename}.pdf")
assert all(year in f for f in (infile, outpdf_dir, filename))
make_dirs(outpdf_dir)
check_overwrite(outpdf_path, overwrite=overwrite)

print("...Opening MyMuonCollection.")
mu_coll = open_pkl(infile)
printer = CanvasPrinter(show_plots=0)
printer.make_plots_pretty()

# Draw all KB2D dpT/pT vs. q*d0 multigraphs.
printer.canv.Print(outpdf_path + "[")
def make_and_draw_mgs(eta_ls, x_lim=None, y_lim=None,
                              scale_by_1divpT=False,
                              scale_by_avgOf1divpT=False,
                              scale_by_muOf1divpT=False,
                              draw_leg=True):
    """For each eta bin, make and draw a mg of dpT/pT vs. q*d0.
    
    NOTE:
    """
    assert sum((scale_by_1divpT, scale_by_avgOf1divpT, scale_by_muOf1divpT)) <= 1
    mg_ls = []
    for eta_min, eta_max in zip(eta_ls[:-1], eta_ls[1:]):
        mg_ls.append(
            mu_coll.make_multigraph(eta_min, eta_max, y_lim, scale_by_1divpT=scale_by_1divpT,
                                                scale_by_avgOf1divpT=scale_by_avgOf1divpT,
                                                scale_by_muOf1divpT=scale_by_muOf1divpT)
        )

    pave_ls = []
    for mg in mg_ls:
        pave_ls.append(
            mu_coll.draw_mg_and_fits(mg, x_lim=x_lim,
                                scale_by_1divpT=scale_by_1divpT,
                                scale_by_avgOf1divpT=scale_by_avgOf1divpT,
                                scale_by_muOf1divpT=scale_by_muOf1divpT,
                                draw_leg=draw_leg)
        )

        printer.canv.Print(outpdf_path)
    return (mg_ls, pave_ls)
print("...Drawing dpT/pT vs. q*d0 plots.")
mg_ls, pave_ls = make_and_draw_mgs(eta_ls, x_lim=[-0.012, 0.006], y_lim=[-0.08, 0.12],
                              scale_by_1divpT=False,
                              scale_by_avgOf1divpT=False,
                              scale_by_muOf1divpT=False,
                              draw_leg=True)
print("...Drawing dpT/pT * 1/<pT> vs. q*d0 plots.")
mg_ls, pave_ls = make_and_draw_mgs(eta_ls, x_lim=[-0.016, 0.006], y_lim=[-0.003, 0.005],
                              scale_by_1divpT=True,
                              scale_by_avgOf1divpT=False,
                              scale_by_muOf1divpT=False,
                              draw_leg=True)
print("...Drawing dpT/pT plots before and after pT corr.")
if print_vals:
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
for kb2d in mu_coll.KinBin2D_dict.values():
    kb2d.make_beforeafterpTcorr_frames()
    kb2d.draw_beforeafterpTcorr()
    printer.canv.Print(outpdf_path)         
    if print_vals:
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
# Eventually make 2D tables here.
printer.canv.Print(outpdf_path + "]")
print("Done.")