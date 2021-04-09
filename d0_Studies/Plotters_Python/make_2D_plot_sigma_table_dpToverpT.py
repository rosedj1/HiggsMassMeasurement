"""
# PURPOSE: 
#  This code reads in values from a pickled dictionary, of the format:
#      Example (Key : Val) pair:
#         ("h_2.3eta2.4_100.0pT150.0_combined_dpToverpTcorr" : 
#             {
#                 'mean_ls'     : [iterfit1, iterfit2, ...],
#                 'mean_err_ls' : [iterfit1, iterfit2, ...],
#                 'std_ls'      : [iterfit1, iterfit2, ...],
#                 'std_err_ls'  : [iterfit1, iterfit2, ...],
#             }
#         )
#  calculates the % improvement per [eta, pT] bin, and puts the improvement
# into a 2D histogram. Finally it draws the histogram to a PDF.
# NOTES: 
#  The values contained in the pickled dict are the best-fit sigma values
#  from iterative Gaussian fits of kinematic distributions, like qd0 or dpT/pT.
#  The keys in an "improvement dictionary" come in pairs: 
#    there is a corrected key and an uncorrected key. Example:
#    h_2.3eta2.4_150.0pT200.0_combined_dpToverpTcorr <--> h_2.3eta2.4_150.0pT200.0_combined_dpToverpT
#  We must be clever to grab these associated values and calculate the % diff. 
#  The main purpose of this code is to make the 2D plot without having to 
#  reproduce all the kinematic distributions and perform the iter. Gaus. fits.
#  FIXME: Some of the variables are deprecated, but verify before eliminating.
#  - Be sure to check all the variables in --- User Parameters ---.
#  - If you get blank cells, try modifying your color_lim range. 
# SYNTAX: python this_script.py
# AUTHOR: Jake Rosenzweig
# CREATED: <=2020-07-24
# UPDATED: 2021-04-08
"""
import os
import pickle
import numpy as np
import ROOT as rt

from Utils_Python.Utils_StatsAndFits import prop_err_on_dsigoversig
from Utils_Python.Utils_Files import check_overwrite, open_pkl, make_dirs
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
from d0_Studies.Plotters_Python.plotter_2D import make_hist_and_errhist
from Utils_ROOT.ROOT_classes import set_TH2F_errs#make_TH2F

#--- User Parameters ---#
year = "2017"
process_derive = r"DY+J/#psi"
method = "AdHoc"
# inpath_sigmadict = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/qd0/MC{year}JpsiDY_2D_plot_qd0dist_gausiterfitsigmas_4unbinnedfits_0p0eta2p4_5p0pT1000p0GeV.pkl"
# inpath_sigmadict = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/20200709_MC{year}_combinesamples_applycorr_fullstats_2p00sigmas__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# inpath_sigmadict = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/pkl_and_csv/MC{year}JpsiDY_gausiterfitsigmas_final_delta_pToverRecpT_5unbinnedfits_0p0eta2p4_5p0pT1000p0GeV.pkl"
inpath_sigmadict = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/FitStats/MC{year}_fitstats_combinesamples_applycorr_fullstats_2p00sigmas__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
# FOR QD0: outfile_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/qd0/MC{year}JpsiDY_2Dplot_qd0dist_gausiterfitsigmas_4unbinnedfits_0p0eta2p4_5p0pT1000p0GeV_final.pdf"
outfile_path =     f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/plots/2Dtables_pTvseta/sigmadpToverpTimprov_MC{year}JpsiDYmuons_6unbinnedfits_0p0eta2p4_5p0pT1000p0GeV_narrowcolor.pdf"
# outfile_path = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/2D_tables_etavspT/MC{year}JpsiDY_2Dplot_dpToverpTimprovement_gausiterfitsigmas_6unbinnedfits_0p0eta2p4_5p0pT1000p0GeV_final.pdf"
overwrite = 0
use_percentage = True

#--- Below is possibly deprecated ---#
# scale_factor = 1.
# kinem = "dpToverpT" #"qd0", used to determine how to access different dictionaries.
# `10000.` should be used for qd0BS (cm -> um). 
# `100.`` should be used for delta_pT (convert to percent)
# WARNING: the dictionary which stores the values may already have scale factor applied. 
#--- Above is possibly deprecated ---#
title_2D_plot_improvement = r"Improvement (#left|#Delta#sigma#right|/#sigma)"
title_2D_plot_before_corr = r"#sigma_{Gauss}(#Deltap_{T}/p_{T}), Unbinned Iter. Fits #it{before corr.}"
if use_percentage:
    title_2D_plot_improvement += r' [%]'
    title_2D_plot_before_corr = title_2D_plot_before_corr.replace(r"),", r") [%],")
# "Unbinned Iter. Fit #sigma_{Gaus}(qd_{0} dist.) (#mum)"
#"Unbinned Iter. Fit #sigma_{Gaus}(#Deltap_{T}/p_{T} dist.) (%)"
title_extra = r" using %s method (MC %s %s muons)" % (method, year, process_derive)

# title_2D_plot_isthisevenused = "Unbinned Iter. Fit #sigma_{Gaus}(qd_{0} dist.) (#mum)"

# Fixing up the titles a little more.
title_2D_plot_improvement += title_extra
# title_2D_plot_before_corr += title_extra
title_2D_plot_after_corr = title_2D_plot_before_corr.replace("before", "after") + title_extra

text_size = 0.85
z_label_size = 0.015
n_contour = 200
color_map = rt.kBird #rt.kTemperatureMap
color_lim_beforeafter_corr = [0, 5.0]
color_lim_improv = [0, 12] #[8, 30] #[0, 12]  # Can also be: None

# Number format: e.g. "6.2f" could be 12.34
# <padding>.<decimals><type>
formatting_beforeafter = ".3f"
formatting_improv = ".3f"
# formatting_beforeafter = "6.3f"
# formatting_improv = "6.2f"

eta_binedge_ls = equal_entry_bin_edges_eta_mod1_wholenum
pT_binedge_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum
#--- Local Functions ---#
def validate(year, inpath_sigmadict, outfile_path):
    """Make sure specified paths and variables are consistent."""
    assert all(year in f for f in [inpath_sigmadict, outfile_path])

def get_associated_keys(dict_sigmas, eta_min, eta_max, pT_min, pT_max):
    """
    Given an eta range and pT range, retrieve the 2 associated keys
    (corr and uncorr).
    """
    # print(f"Testing: eta_min={eta_min}, eta_max={eta_max}, pT_min={pT_min}, pT_max={pT_max}")
    keys = [k for k in dict_sigmas.keys() if (f"_{pT_min}" in k) and (f"{pT_max}_" in k) and (f"_{eta_min}" in k) and (f"{eta_max}_" in k)]
    assert len(keys) == 2
    key_ls = sorted(keys)  # Goes: [h_2.3eta2.4_1..._dpToverpT, h_2.3eta2.4_1..._dpToverpTcorr]
    return key_ls

def get_sigma_and_err(key, dict_sigmas):
    """Return the best-fit iterated Gaussian fit sigma from the dict."""
    sig = dict_sigmas[key]['std_ls'][-1]  # Last element is the last/best fit value.
    sig_err = dict_sigmas[key]['std_err_ls'][-1]
    return (sig, sig_err)

def fill_hists_with_dict_dpToverpTvals(eta_binedge_ls, pT_binedge_ls,
                                      h2_before_corr, h2_before_corr_err,
                                      h2_after_corr, h2_after_corr_err,
                                      h2_improv, h2_improv_err,
                                      dict_sigmas, use_percentage=True):
    """
    Fill three sets of hists with the dpT/pT fit stats:
    - Set 1: sigma(dpT/pT) before correction and its error
    - Set 2: sigma(dpT/pT) after correction and its error
    - Set 3: Improvement of sigma(dpT/pT) and its error
    
    NOTE:
    Finds the 2 associated keys (corr and uncorr), 
    calculates the % improvement, and puts the value 
    into a specific cell of the 2D hist.
    The y bins are eta vals and the x bins are pT vals.
    
    Example (Key : Val) pair:
        ("h_2.3eta2.4_100.0pT150.0_combined_dpToverpTcorr" : 
            {
                'mean_ls'     : [iterfit1, iterfit2, ...],
                'mean_err_ls' : [iterfit1, iterfit2, ...],
                'std_ls'      : [iterfit1, iterfit2, ...],
                'std_err_ls'  : [iterfit1, iterfit2, ...],
            }
        )

    Parameters
    ----------
    scale : bool
        If True, dpT/pT 2D tables will show percentages on cells.
    """
    scale = 100.0 if use_percentage else 1.0
    for eta_min,eta_max in zip(eta_binedge_ls[:-1], eta_binedge_ls[1:]):
        for pT_min,pT_max in zip(pT_binedge_ls[:-1], pT_binedge_ls[1:]):
            key, key_corr = get_associated_keys(dict_sigmas, eta_min, eta_max, pT_min, pT_max)
            sigma, sigma_err = get_sigma_and_err(key, dict_sigmas)
            sigma_corr, sigma_corr_err = get_sigma_and_err(key_corr, dict_sigmas)
            improv = (sigma_corr - sigma) / sigma
            improv_err = prop_err_on_dsigoversig(sigma, sigma_corr, sigma_err, sigma_corr_err)
            midpt_pT = np.mean((pT_min, pT_max))
            midpt_eta = np.mean((eta_min, eta_max))
            # print(f"...Filling 2D hist: midpt_pT={midpt_pT}, midpt_eta={midpt_eta}, sigma={sigma}")
            h2_before_corr.Fill(midpt_pT, midpt_eta, sigma * scale)
            h2_before_corr_err.Fill(midpt_pT, midpt_eta, sigma_err * scale)
            # print(f"...Filling 2D hist: midpt_pT={midpt_pT}, midpt_eta={midpt_eta}, sigma_corr={sigma_corr}")
            h2_after_corr.Fill(midpt_pT, midpt_eta, sigma_corr * scale)
            h2_after_corr_err.Fill(midpt_pT, midpt_eta, sigma_corr_err * scale)
            # print(f"...Filling 2D hist: midpt_pT={midpt_pT}, midpt_eta={midpt_eta}, improvement(%)={improv}")
            h2_improv.Fill(midpt_pT, midpt_eta, -1*(improv * scale))
            h2_improv_err.Fill(midpt_pT, midpt_eta, -1*(improv_err * scale))

def get_ranges(key):
    """
    Parse a str key and convert it to the corresponding eta_range and pT_range.
    Example of key:
        2p2eta2p3_5p0pT7p0_sigmawitherr_delta_pToverRecpT
    After conversion:
        eta_range = [2.2, 2.3]
        pT_range = [5.0, 7.0]        
    """
    eta_str, pT_str = key.split("_")[0:2]
    eta_range_str = eta_str.split("eta")
    pT_range_str = pT_str.split("pT")
    eta_range = [float(num.replace("p",".")) for num in eta_range_str]
    pT_range = [float(num.replace("p",".")) for num in pT_range_str]
    return eta_range, pT_range

def make_pdf(h_before_corr, h_after_corr, h_improve, outfile_path):
    rt.gStyle.SetPalette(color_map) 
    rt.gROOT.SetBatch(rt.kTRUE)
    rt.gStyle.SetOptStat(0)
    c1 = rt.TCanvas()
    c1.SetLogx(True)
    c1.Print(outfile_path + "[")
    rt.gStyle.SetPaintTextFormat(formatting_beforeafter)

    h_before_corr.SetTitleOffset(1.3)
    h_before_corr.SetMarkerSize(text_size)
    h_before_corr.Draw("colz text e1")
    c1.Print(outfile_path)

    h_after_corr.SetTitleOffset(1.3)
    h_after_corr.SetMarkerSize(text_size)
    h_after_corr.Draw("colz text e1")
    c1.Print(outfile_path)
    rt.gStyle.SetPaintTextFormat(formatting_improv)

    h_improve.SetTitleOffset(1.3)
    h_improve.SetMarkerSize(text_size)
    h_improve.Draw("colz text e1")
    c1.Print(outfile_path)
    c1.Print(outfile_path + "]")
    # c1.SaveAs(outfile_path)

if __name__ == "__main__":
    make_dirs(os.path.dirname(outfile_path))
    validate(year, inpath_sigmadict, outfile_path)
    check_overwrite(outfile_path, overwrite=overwrite)
    dict_sigmas = open_pkl(inpath_sigmadict)
    h2_before_corr, h2_before_corr_err = make_hist_and_errhist("before_corr", title=title_2D_plot_before_corr, 
              n_binsx=pT_binedge_ls, x_label="p_{T}", x_units="GeV",
              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
              z_min=color_lim_beforeafter_corr[0], z_max=color_lim_beforeafter_corr[1],
              z_label_size=z_label_size, n_contour=n_contour, extra_name="err")
    h2_after_corr, h2_after_corr_err = make_hist_and_errhist("after_corr", title=title_2D_plot_after_corr, 
              n_binsx=pT_binedge_ls, x_label="p_{T}", x_units="GeV",
              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
              z_min=color_lim_beforeafter_corr[0], z_max=color_lim_beforeafter_corr[1],
              z_label_size=z_label_size, n_contour=n_contour, extra_name="err")
    h2_improv, h2_improv_err = make_hist_and_errhist("improvement", title=title_2D_plot_improvement, 
              n_binsx=pT_binedge_ls, x_label="p_{T}", x_units="GeV",
              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
              z_min=color_lim_improv[0], z_max=color_lim_improv[1],
              z_label_size=z_label_size, n_contour=n_contour, extra_name="err")

# if "dpToverpT" in kinem:
    fill_hists_with_dict_dpToverpTvals(eta_binedge_ls, pT_binedge_ls,
                                       h2_before_corr, h2_before_corr_err,
                                       h2_after_corr, h2_after_corr_err,
                                       h2_improv, h2_improv_err,
                                       dict_sigmas, use_percentage=use_percentage)
    set_TH2F_errs(h2_before_corr, h2_before_corr_err)
    set_TH2F_errs(h2_after_corr, h2_after_corr_err)
    set_TH2F_errs(h2_improv, h2_improv_err)
    print(f"Making PDF of 2D hists...")
    make_pdf(h2_before_corr, h2_after_corr, h2_improv, outfile_path)