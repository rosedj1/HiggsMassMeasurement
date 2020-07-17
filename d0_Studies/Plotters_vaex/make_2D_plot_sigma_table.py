"""
# PURPOSE: 
#  This code reads in values from a pickled dictionary and puts those 
#  values into a 2D histogram. Finally it draws the histogram to a PDF.
# NOTES: 
#  The values contained in the pickled dict are the best-fit sigma values
#  from iterative Gaussian fits of kinematic distributions, like dpT/pT or qd0. 
#  The main purpose of this code is to make the 2D plot without having to 
#  reproduce all the kinematic distributions and perform the iter. Gaus. fits.
#  - Be sure to check all the variables in --- User Parameters ---.
#  - If you get blank cells, try modifying your color_lim range. 
# SYNTAX: python this_script.py
# AUTHOR: Jake Rosenzweig
# EDITED: 2020-07-16
"""
import pickle
import numpy as np
import ROOT as r

from Utils_Python.Utils_Files import check_overwrite

#--- User Parameters ---#
year = "2017"
inpath_sigmadict = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/qd0/MC{year}JpsiDY_2D_plot_qd0dist_gausiterfitsigmas_4unbinnedfits_0p0eta2p4_5p0pT1000p0GeV.pkl"
outfile_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/qd0/MC{year}JpsiDY_2Dplot_qd0dist_gausiterfitsigmas_4unbinnedfits_0p0eta2p4_5p0pT1000p0GeV_final.pdf"
overwrite = False

# Only used to decide how text is displayed on cells and in file-naming.
kinem = "qd0BS" #"delta_pToverRecpT"

# `10000.` should be used for qd0BS (cm -> um). 
# `100.`` should be used for delta_pT (convert to percent)
# WARNING: the dictionary which stores the values may already have scale factor applied. 
scale_factor = 1.
title_2D_plot = "Unbinned Iter. Fit #sigma_{Gaus}(qd_{0} dist.) (#mum)"
#"Unbinned Iter. Fit #sigma_{Gaus}(#Deltap_{T}/p_{T} dist.) (%)"
color_lim = [8, 30]#[0, 12]  # Can also be: None
#--- Local Functions ---#
def get_dict_of_sigmas(inpath_dict):
    """Return the dictionary of sigma values from pickle file."""
    with open(inpath_dict, "rb") as pkl:
        dct = pickle.load(pkl)
    return dct

def extract_binedge_lists(dict_sigmas):
    """Get eta and pT bin edge lists from dictionary."""
    return dict_sigmas["binedges_eta"], dict_sigmas["binedges_pT"]

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
    
def fill_hist_with_dict_vals(h_2d, dict_sigmas, scale_factor=1):
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

def make_2D_hist(kinem, title_2D_plot, eta_ls, pT_ls):
    """
    Make 2D hist which will hold the best-fit sigma values
    from iterated gaus fits of the kinematic distributions. 
    
    The x-bin edges are the pT bin edges.
    The y-bin edges are the eta bin edges.
    """
    hist_2D = r.TH2F("h_2d", title_2D_plot, 
              len(pT_ls)-1, np.array(pT_ls, dtype=float), 
              len(eta_ls)-1, np.array(eta_ls, dtype=float) )

    r.gStyle.SetOptStat(0)

    hist_2D.SetXTitle("p_{T} [GeV]")
    hist_2D.SetYTitle("#left|#eta#right|")

    hist_2D.GetXaxis().SetTitleOffset(1.3)
    hist_2D.SetContour(200)
    if "qd0" in kinem:
        r.gStyle.SetPaintTextFormat("6.0f")  # Number format. Don't measure sub-micron.
    elif "delta_pT" in kinem:
        r.gStyle.SetPaintTextFormat("6.2f")  # Number format.
    return hist_2D

def make_pdf(kinem, hist, outfile_path, preview_plots=False):
    r.gROOT.SetBatch(not preview_plots)
    c1 = r.TCanvas()
    c1.SetLogx(True)
    hist.Draw("colz text")
    c1.SaveAs(outfile_path)

if __name__ == "__main__":
    check_overwrite(outfile_path, overwrite=overwrite)
    dict_sigmas = get_dict_of_sigmas(inpath_sigmadict)
    eta_ls, pT_ls = extract_binedge_lists(dict_sigmas)
    h_2d = make_2D_hist(kinem, title_2D_plot, eta_ls, pT_ls)
    fill_hist_with_dict_vals(h_2d, dict_sigmas, scale_factor=scale_factor)

    print(f"Making PDF of 2D hist...")
    if color_lim is not None:
        h_2d.GetZaxis().SetRangeUser(*color_lim)
    make_pdf(kinem, h_2d, outfile_path)