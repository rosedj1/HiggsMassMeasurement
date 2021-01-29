"""
# PURPOSE: 
#  This code reads in values from a pickled dictionary, calculates the 
#  % improvement per [eta, pT] bin, and puts the improvement into 
#  a 2D histogram. Finally it draws the histogram to a PDF.
# NOTES: 
#  The values contained in the pickled dict are the best-fit sigma values
#  from iterative Gaussian fits of kinematic distributions, like qd0 or dpT/pT.
#  The keys in an "improvement dictionary" come in pairs: 
#    there is a corrected key and an uncorrected key. Example:
#    h_2.3eta2.4_150.0pT200.0_combined_dpToverpTcorr <--> h_2.3eta2.4_150.0pT200.0_combined_dpToverpT
#  We must be clever to grab these associated values and calculate the % diff. 
#  The main purpose of this code is to make the 2D plot without having to 
#  reproduce all the kinematic distributions and perform the iter. Gaus. fits.
#  FIXME: Some of the variables are deprecated, but make sure: scale_factor, etc.
#  - Be sure to check all the variables in --- User Parameters ---.
#  - If you get blank cells, try modifying your color_lim range. 
# SYNTAX: python this_script.py
# AUTHOR: Jake Rosenzweig
# EDITED: 2020-07-24
"""
import pickle
import numpy as np
import ROOT as r

from Utils_Python.Utils_Files import check_overwrite
#--- User Parameters ---#
year = "2017"
# inpath_sigmadict = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/qd0/MC{year}JpsiDY_2D_plot_qd0dist_gausiterfitsigmas_4unbinnedfits_0p0eta2p4_5p0pT1000p0GeV.pkl"
inpath_sigmadict = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/20200709_MC{year}_combinesamples_applycorr_fullstats_2p00sigmas__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
outfile_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/2D_tables_etavspT/MC{year}JpsiDY_2Dplot_dpToverpTimprovement_gausiterfitsigmas_6unbinnedfits_0p0eta2p4_5p0pT1000p0GeV_final.pdf"
overwrite = True

kinem = "qd0"  # dpToverpT, used to determine how to access different dictionaries.
# `10000.` should be used for qd0BS (cm -> um). 
# `100.`` should be used for delta_pT (convert to percent)
# WARNING: the dictionary which stores the values may already have scale factor applied. 
scale_factor = 1.
title_2D_plot_improvement = r"Improvement in #sigma_{Gaus}(#Deltap_{T}/p_{T} dist.) [%] from Unbinned Iter. Fit"
title_2D_plot_before_corr = r"#sigma_{Gaus}(#Deltap_{T}/p_{T} dist.) [GeV] from Unbinned Iter. Fit #it{before corr.}"
# "Unbinned Iter. Fit #sigma_{Gaus}(qd_{0} dist.) (#mum)"
#"Unbinned Iter. Fit #sigma_{Gaus}(#Deltap_{T}/p_{T} dist.) (%)"
title_extra = ", using ad hoc p_{T} corr. (MC %s)" % year

# Fixing up the titles a little more.
title_2D_plot_improvement += title_extra
title_2D_plot_before_corr += title_extra
title_2D_plot_after_corr = title_2D_plot_before_corr.replace("before", "after")

color_lim_beforeafter_corr = [0, 0.05]
color_lim_improv = [0, 12] #[8, 30] #[0, 12]  # Can also be: None

# Number format: e.g. "6.2f" could be 12.34
# <padding>.<decimals><type>
formatting_beforeafter = "6.3f"
formatting_improv = "6.2f"

eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
pT_binedge_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]
#--- Local Functions ---#
def get_dict_of_sigmas(inpath_dict):
    """Return the dictionary of sigma values from pickle file."""
    with open(inpath_dict, "rb") as pkl:
        dct = pickle.load(pkl)
    return dct

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

def get_sigma(key, dict_sigmas):
    """Return the best-fit iterated Gaussian fit sigma from the dict."""
    return dict_sigmas[key]['std_ls'][-1]  # Last element is the last/best fit value.

def calc_improvement(sigma, sigma_corr):
    """Calculate the absolute relative difference."""
    return abs(sigma_corr - sigma) / sigma

def fill_hist_with_dict_dpToverpTvals(eta_binedge_ls, pT_binedge_ls, h_2d_improv, h_2d_before_corr, h_2d_after_corr, dict_sigmas):
    """
    Find the 2 associated keys (corr and uncorr), 
    calculate the % improvement, and put the value 
    into a specific cell of the 2D hist.

    Also puts the sigma before-correction and sigma after-correction
    into their corresponding hists.
    
    Example (Key : Val) pair:
        ("h_2.3eta2.4_100.0pT150.0_combined_dpToverpTcorr" : 4-tuple)
    """
    for eta_min,eta_max in zip(eta_binedge_ls[:-1], eta_binedge_ls[1:]):
        for pT_min,pT_max in zip(pT_binedge_ls[:-1], pT_binedge_ls[1:]):

            key, key_corr = get_associated_keys(dict_sigmas, eta_min, eta_max, pT_min, pT_max)
            sigma = get_sigma(key, dict_sigmas)
            sigma_corr = get_sigma(key_corr, dict_sigmas)
            improv = calc_improvement(sigma, sigma_corr) * 100.0

            midpt_pT = (pT_min + pT_max) / 2.0
            midpt_eta = (eta_min + eta_max) / 2.0
            print(f"...Filling 2D hist: midpt_pT, midpt_eta, sigma=({midpt_pT}, {midpt_eta}, {sigma})")
            h_2d_before_corr.Fill(midpt_pT, midpt_eta, sigma)
            print(f"...Filling 2D hist: midpt_pT, midpt_eta, sigma_corr=({midpt_pT}, {midpt_eta}, {sigma_corr})")
            h_2d_after_corr.Fill(midpt_pT, midpt_eta, sigma_corr)
            print(f"...Filling 2D hist: midpt_pT, midpt_eta, improvement(%)=({midpt_pT}, {midpt_eta}, {improv})")
            h_2d_improv.Fill(midpt_pT, midpt_eta, improv)

def fill_hist_with_dict_qd0vals(h_2d_improv, h_2d_before_corr, h_2d_after_corr, dict_sigmas, scale_factor=scale_factor)):
    """
    Put the value of each key into a specific cell of the 2D hist.
    
    Example (Key : Val) pair:
        ("2p2eta2p3_5p0pT7p0_sigmawitherr_delta_pToverRecpT" : [bestfit_sigma, bestfit_sigma_err])
    """
    for key,sigma_ls in dict_sigmas.items():
        if "binedges" in key:
            continue
        bestfit_sigma = sigma_ls[0]
        bestfit_sigma *= scale_factor

        eta_range, pT_range = get_ranges(key)
        midpt_pT = sum(pT_range) / 2.
        midpt_eta = sum(eta_range) / 2.

        print(f"...Filling 2D hist: midpt_pT, midpt_eta, bestfit_sigma=({midpt_pT}, {midpt_eta}, {bestfit_sigma})")
        h_2d.Fill(midpt_pT, midpt_eta, bestfit_sigma)

def make_2D_hist(name, title, eta_ls, pT_ls, color_lim, formatting="6.2f"):
    """
    Make 2D hist which will hold the best-fit sigma values
    from iterated gaus fits of the kinematic distributions. 
    
    The x-bin edges are the pT bin edges.
    The y-bin edges are the eta bin edges.
    """
    h2d = r.TH2F(name, 
                 title, 
                 len(pT_ls)-1, np.array(pT_ls, dtype=float), 
                 len(eta_ls)-1, np.array(eta_ls, dtype=float) 
                 )

    r.gStyle.SetOptStat(0)

    h2d.SetXTitle("p_{T} [GeV]")
    h2d.SetYTitle("#left|#eta#right|")

    h2d.GetXaxis().SetTitleOffset(1.3)
    h2d.SetContour(200)  # Use N different colors to smooth out color bar gradient.
    h2d.GetZaxis().SetRangeUser(*color_lim)  # 
    h2d.SetLabelSize(0.015, "Z")
    return h2d

def make_pdf(h_before_corr, h_after_corr, h_improve, outfile_path):
    r.gROOT.SetBatch(r.kTRUE)
    c1 = r.TCanvas()
    c1.SetLogx(True)
    c1.Print(outfile_path + "[")
    r.gStyle.SetPaintTextFormat(formatting_beforeafter)
    h_before_corr.Draw("colz text")
    c1.Print(outfile_path)
    h_after_corr.Draw("colz text")
    c1.Print(outfile_path)
    r.gStyle.SetPaintTextFormat(formatting_improv)
    h_improve.Draw("colz text")
    c1.Print(outfile_path)
    c1.Print(outfile_path + "]")
    # c1.SaveAs(outfile_path)

if __name__ == "__main__":
    check_overwrite(outfile_path, overwrite=overwrite)
    dict_sigmas = get_dict_of_sigmas(inpath_sigmadict)
    h_2d_before_corr = make_2D_hist("before_corr", title_2D_plot_before_corr,
                                    eta_binedge_ls, pT_binedge_ls, color_lim_beforeafter_corr, formatting_beforeafter)
    h_2d_after_corr = make_2D_hist("after_corr", title_2D_plot_after_corr, 
                                    eta_binedge_ls, pT_binedge_ls, color_lim_beforeafter_corr, formatting_beforeafter)
    h_2d_improv = make_2D_hist("improvement", title_2D_plot_improvement, 
                                    eta_binedge_ls, pT_binedge_ls, color_lim_improv, formatting_improv)

    if "qd0" in kinem:
        fill_hist_with_dict_qd0vals(h_2d_improv, h_2d_before_corr, h_2d_after_corr, dict_sigmas, scale_factor=scale_factor)
    elif "dpToverpT" in kinem:
        fill_hist_with_dict_vals(eta_binedge_ls, pT_binedge_ls, h_2d_improv, h_2d_before_corr, h_2d_after_corr, dict_sigmas)
    else:
        print(f"Not filling plots. You should probably look into this.")

    print(f"Making PDF of 2D hists...")
    make_pdf(h_2d_before_corr, h_2d_after_corr, h_2d_improv, outfile_path)