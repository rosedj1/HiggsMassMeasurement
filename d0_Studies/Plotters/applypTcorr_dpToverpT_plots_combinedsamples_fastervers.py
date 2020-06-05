"""
# Purpose: 
#   This script scans over User-specified samples (DY, Jpsi, Upsilon) and
#   selects muons which pass selection criteria and makes dpT/pT histograms. 
#   A dpT/pT distribution is made for all eta, pT bins specified. 
#   Iterated Gaussian fits are then performed on these distributions.
#   Finally the pT_rec is "corrected" using results from the 
#   best-fit parameters from q*d0 studies. 
# Notes: 
#   All muons scanned over, whether DY or Jpsi, etc., are treated as 
#   "individual/independent muons", i.e. they are all put into the same
#   dpT/pT distributions, depending on the (eta, pT) bin of the muon.
#   Example key of hist_dict = "h_0.0eta0.2_27.0pT38.0_Jpsi_dpToverpT"
# Syntax  : python <script>.py
# Author  : Jake Rosenzweig
# Updated : 2020-06-04
"""
import os
import pickle
import ROOT
import numpy as np

from ROOT import TFile, TH1F, TCanvas, gROOT, kTRUE, TLatex

from Utils_Python.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Selections import Selector
from Utils_Python.Utils_Physics import calc_dphi, calc_dR
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import fixOverlay, setTDRStyle, tdrGrid

from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit
from Utils_ROOT.ROOT_fns import (fill_dpToverpT_hist, make_hist_lookuptable, 
                                 fill_dict_of_dpToverpT_hists)

from d0_Studies.d0_Utils.d0_fns import calc_num_bins, print_header_message

#-----------------------#
#----- User Params -----#
#-----------------------#
infile_path_MC_2017_DY = "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_DY.root"
infile_path_MC_2017_Jpsi = "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Jpsi.root"
infile_path_MC_2017_Upsilon = "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Upsilon.root"

infile_path_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/d0_pT_corrfactors_0p0eta2p4__5p0pT1000p0.pkl"

outdir_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/Testing/"
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/hists_dpToverpT/different_etapT_bins/"

filename = "20200604_combinesamples_applycorr_fastvers_fullstats"
overwrite = False

n_evts = -1
verbose = False

do_iter_gaus_fit = True
do_binned_fit = True
only_draw_last = True
apply_d0_pT_corrections = True

apply_dR_cut = True
apply_m2l_cut = True

iters = 6
num_sigmas = 2.0

color_before_corr = ROOT.kBlue
color_after_corr = ROOT.kRed

# For some dumb reason, I've programmed it 
# to explicitly go "Jpsi", "Upsilon"...
# Due to the zip()!!!
sample_ls = ["Jpsi","Upsilon","DY","combined"]

# Note: To apply pT correction factors, 
# you must have a hist_type that contains "corr".
hist_type_ls = ["dpToverpT", "dpToverpTcorr"]#"m2muplus", "m2muminus", "m2mu"

# Since samples get combined, make sure 
# all samples have same binning info. 
# Everyone uses this binning info
bin_info_dpToverpT = [-0.08, 0.08, 0.001]  # [bin_min, bin_max, bin_width]

eta_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
pT_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]
# eta_ls = [0.0, 0.2]
# pT_ls  = [20.0, 27.0]

draw_stats = False
show_plots_to_screen = False

# Only need to configure this if you introduce new samples,
# or change file paths. 
f_Jpsi = TFile.Open(infile_path_MC_2017_Jpsi)
t_Jpsi = f_Jpsi.Get("passedEvents")
f_Upsilon = TFile.Open(infile_path_MC_2017_Upsilon)
t_Upsilon = f_Upsilon.Get("passedEvents")
f_DY = TFile.Open(infile_path_MC_2017_DY)
t_DY = f_DY.Get("passedEvents")
sample_dict = {
    # "bin_info" is deprecated.
    "Jpsi"    : { "bin_info"   : {"dpToverpT"      : bin_info_dpToverpT,
                                  "dpToverpTcorr"  : bin_info_dpToverpT },
                  "fancy_font" :  {"ROOT" : "J/#psi", "LaTeX" : r'$J/\psi$'},
                  "tree"       : t_Jpsi,
                },

    "Upsilon" : {"bin_info"   :  {"dpToverpT"      : bin_info_dpToverpT,
                                  "dpToverpTcorr"  : bin_info_dpToverpT},
                  "fancy_font" : {"ROOT" : "#Upsilon", "LaTeX" : r'$\Upsilon$'},
                  "tree"       : t_Upsilon,
                },
    "DY"      : {"bin_info"   :  {"dpToverpT"      : bin_info_dpToverpT,
                                  "dpToverpTcorr"  : bin_info_dpToverpT},
            "fancy_font" : {"ROOT" : "Z#rightarrow#mu^+#mu^-", "LaTeX" : r'$Z\rightarrow\mu^+\mu^-$'},
            "tree"       : t_DY,
                },
    "combined" : {"bin_info"   :  {"dpToverpT"      : bin_info_dpToverpT,
                                  "dpToverpTcorr"  : bin_info_dpToverpT},
            "fancy_font" : {"ROOT" : "all muons from DY, J/#psi, #Upsilon", 
                            "LaTeX" : ""},
            "tree"       : None,
                }
}
#----------------------#
#----- Automatons -----#
#----------------------#
# Quick checks.
if "combined" in sample_ls:
    # Must also specify at least one sample (like, "Jpsi"). 
    assert len(sample_ls) > 1
if (apply_d0_pT_corrections):
    msg = "You have to have a hist_type with 'corr' in it."
    assert any(["corr" in x for x in hist_type_ls]), msg
    # Load up the correction factors. 
    with open(infile_path_pkl, "rb") as f:
        pT_corr_factor_dict = pickle.load(f)

# Make directories, if need be.
makeDirs(outdir_plots)
makeDirs(outdir_pkl)

# Handle naming of files. 
full_extra = (
    f"_{num_sigmas:.2f}sigmas"
    f"__{min(eta_ls):.1f}_eta_{max(eta_ls):.1f}"
    f"__{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
)
full_extra = make_str_title_friendly(full_extra)

fullfilename = filename + full_extra

fullpath_pkl = os.path.join(outdir_pkl, fullfilename + ".pkl")
for samp in sample_ls:
    # Check if PDF for this sample exists.
    filename_samp = f"{fullfilename}_{samp}.pdf"
    fullpath_file_samp = os.path.join(outdir_plots, filename_samp)
    check_overwrite(fullpath_file_samp, overwrite)

# See if pickle already exists.
# Could check for PDFs, but trickier since there are many
# produced each time this script is run.
check_overwrite(fullpath_pkl, overwrite)

# Everything seems good. Start the analysis.

# Plotting info.
if not (show_plots_to_screen): gROOT.SetBatch(kTRUE)
tdrStyle = setTDRStyle()
tdrGrid(tdrStyle, gridOn=True)

dpToverpT_str = "(p_{T}^{REC} - p_{T}^{GEN}) / p_{T}^{GEN}"

#---------------------------#
#----- Local Functions -----#
#---------------------------#
class HistDict:
    def __init__(self, hist_dict, sample, hist_type):
        self.hist_dict = hist_dict
        self.sample = sample
        self.hist_type = hist_type
        
    def get_TH1F_from_dict(self, eta_bin_ls, pT_min, pT_max):
        eta_min = eta_bin_ls[0]
        eta_max = eta_bin_ls[1]
        pT_min = pT_bin_ls[0]
        pT_max = pT_bin_ls[1]

        key = f"h_{eta_min}eta{eta_max}_{pT_min}pT{pT_max}_{self.sample}_{self.hist_type}"
        h = self.hist_dict[key]
        assert isinstance(h, ROOT.TH1F)
        return h

class OrgHistDict:
    def __init__(self):
        self.dict_of_HistDicts = {}
    
    def store_HistDict(self, HistDict_):
        sample = HistDict_.sample
        hist_type = HistDict_.hist_type
        key = f"{sample}_{hist_type}"
        self.dict_of_HistDicts[key] = HistDict_ 

    def get_HistDict(self, sample, hist_type):
        key = f"{sample}_{hist_type}"
        return self.dict_of_HistDicts[key]

    def get_kinbin_TH1F(self, eta_bin_ls, pT_bin_ls, sample, hist_type):
        etamin = eta_bin_ls[0]
        etamax = eta_bin_ls[1]
        pTmin = pT_bin_ls[0]
        pTmax = pT_bin_ls[1]

        hdict = self.get_HistDict(sample, hist_type).hist_dict
        key = f"h_{etamin}eta{etamax}_{pTmin}pT{pTmax}_{sample}_{hist_type}"
        return hdict[key]

    def make_ls_TH1F_samekinbin_diffsample(self, eta_bin_ls, pT_bin_ls, sample_ls, hist_type):
        assert "combined" not in sample_ls
        # Get all TH1F that are in this kinbin for all individual samples.
        th1f_ls = []
        for sample in sample_ls:
            th1f_ls.append( self.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, sample, hist_type) )
            # key = f"{sample}_{hist_type}"
            # hdict = self.get_HistDict(sample, hist_type)[key]
        return th1f_ls

    def add_diffsamples_kinbin_TH1F(self, h_comb, h_ls):
        """
        Add all TH1F within one kinbin into single, combined TH1F.
        """
        for hist in h_ls:
            h_comb.Add(hist)

    def fill_combined_kinbin_hist(self, eta_bin_ls, pT_bin_ls, sample_ls, hist_type):
        """
        Based on given parameters, find the corresponding combined TH1F 
        (from the "combined" hist_dict) in this organizer, and fill the TH1F 
        using all entries from individual samples.
        """
        th1f_ls = self.make_ls_TH1F_samekinbin_diffsample(eta_bin_ls, pT_bin_ls, sample_ls, hist_type)
        # Get the empty TH1F that corresponds to all the given parameters.
        h_comb = self.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, "combined", hist_type)
        assert h_comb.GetEntries() == 0

        # Fill the empty, combined hist.
        self.add_diffsamples_kinbin_TH1F(h_comb, th1f_ls)

def make_iter_binedge_ls(ls): 
    """
    Turns a list from this: [1,2,3,4]
    into this: 
      [ [1,2],
        [2,3],
        [3,4] ]

    Good for iterative over neighboring elements.

    ...and I just realized I can simply use a `zip()` while looping:
    for eta_min, eta_max in zip(myls[:-1], myls[1:]): 
        print(eta_min, eta_max) 
    """
    ls_beg = ls[:-1] 
    ls_end = ls[1:] 
    combine = zip(ls_beg, ls_end) 
    return list(combine) 

def get_hist_name_parts(h_name):
    parts = h_name.split("_")
    eta_piece = parts[1]
    pT_piece = parts[2]
    sample = parts[3]
    h_type = parts[4:]

    eta_min = float(eta_piece.split("eta")[0])
    eta_max = float(eta_piece.split("eta")[1])
    pT_min = float(pT_piece.split("pT")[0])
    pT_max = float(pT_piece.split("pT")[1])
    return eta_min, eta_max, pT_min, pT_max, sample, h_type
#----------------#
#----- Main -----#
#----------------#
# Initialize all hist dicts.
org_hist_dict = OrgHistDict()
eta_bin_ls_2D = make_iter_binedge_ls(eta_ls)
pT_bin_ls_2D = make_iter_binedge_ls(pT_ls)

for h_type in hist_type_ls:
    for sample in sample_ls:
        print(f"#----- Making {sample} hist_dict: {h_type} -----#")
        h_dict = make_hist_lookuptable(eta_ls=eta_ls,pT_ls=pT_ls,sample=sample,
                                       hist_type=h_type, bin_info_ls=bin_info_dpToverpT)
        print("hist_dict made.")
        
        # Store them in an organizer.
        org_hist_dict.store_HistDict(HistDict(hist_dict=h_dict, 
                                                hist_type=h_type, 
                                                sample=sample) 
                                    )

# Fill each sample hist dict. 
individ_sample_ls = [x for x in sample_ls if "combined" not in x]
for sample in individ_sample_ls:
    hist_dict_ls = []
    for hist_type in hist_type_ls:
        hist_dict_ls.append( org_hist_dict.get_HistDict(sample, hist_type).hist_dict )

    # Now enter each root file ONLY ONCE (maximize efficiency).
    # Fill the sample hist and corresponding combined hist in parallel.
    print(f"#----- Filling {sample} hist_dict: {hist_type_ls} -----#")
    fill_dict_of_dpToverpT_hists(
        tree=sample_dict[sample]["tree"], 
        hist_dict=hist_dict_ls, 
        hist_type=hist_type_ls,
        sample=sample, 
        n_evts=n_evts,
        eta_binedge_ls=eta_ls, pT_binedge_ls=pT_ls,
        apply_dR_cut=apply_dR_cut, apply_m2l_cut=apply_m2l_cut, verbose=verbose,
        apply_d0_pT_corrections=apply_d0_pT_corrections, 
        pT_corr_factor_dict=pT_corr_factor_dict
    )
    print("Hist dict filled.")

# Each (hist_type, sample) has its own dict of TH1F. 
# Only samples are filled so far, not the combined TH1Fs.
for h_type in hist_type_ls:
    for eta_bin_ls in eta_bin_ls_2D:
        for pT_bin_ls in pT_bin_ls_2D:
            org_hist_dict.fill_combined_kinbin_hist(eta_bin_ls, pT_bin_ls, individ_sample_ls, h_type)
            
            # While we're here, check that number of events is consistent
            # in the combined hists. 
            num_entries = 0
            for sample in individ_sample_ls:
                num_entries += org_hist_dict.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, sample, h_type).GetEntries()
            assert num_entries == org_hist_dict.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, "combined", h_type).GetEntries()

# Combined hists are all filled. 

dict_of_fit_stats_dict = {}
total_entries = 0
total_entries_corr = 0

# for sample in sample_ls:
    # Create PDF name.
filename_sample = f"{fullfilename}.pdf"  # fullfilename is entire filename, except extension.
fullpath_pdf_ = os.path.join(outdir_plots, filename_sample)
# Make a PDF for each sample (including "combined"!)
HD      = org_hist_dict.get_HistDict("combined", "dpToverpT")
HD_corr = org_hist_dict.get_HistDict("combined", "dpToverpTcorr")

canv = TCanvas(filename, filename, 600, 600)
print(f"Opening up big PDF...")
canv.Print(fullpath_pdf_ + "[")

# Plot similar histograms.
# for eta_bin_ls in eta_bin_ls_2D:
#     for pT_bin_ls in pT_bin_ls_2D:
for h,h_corr in zip(HD.hist_dict.values(), HD_corr.hist_dict.values()):
    # h = hdictobj.get_TH1F_from_dict(eta_bin_ls, pT_bin_ls)
    # h_corr = hdictobj.get_TH1F_from_dict(eta_bin_ls, pT_bin_ls)

    # Quick checks.
    h_name_parts = get_hist_name_parts(h.GetName())
    h_corr_name_parts = get_hist_name_parts(h_corr.GetName())
    for part,part_corr in zip(h_name_parts[:-1],h_corr_name_parts[:-1]):
        # Make sure everything agrees except the hist_type. 
        assert part == part_corr

    # Pretty up the plots.
    latex_str = sample_dict[HD.sample]["fancy_font"]["ROOT"]
    title_str = (
        f"MC 2017 {latex_str}    binning=("
        f"{h_name_parts[0]:.1f} < #left|#eta#right| < {h_name_parts[1]:.1f}    "
        )
    title_str += r"%.1f < p_{T} < %.1f GeV), " % (h_name_parts[2], h_name_parts[3])
    title_str += "fit_range=#mu#pm{:.2f}#sigma".format(num_sigmas)
    
    print(f"Information about these hists:\n"
            f"  name={h.GetName()},       name_corr={h_corr.GetName()}"
            f"  entries={h.GetEntries()}, entries_corr={h_corr.GetEntries()}\n")
    
    print(f"...Performing iterative gaus fit on hists...")
    canv.cd()

    x_roofit = ROOT.RooRealVar("x_roofit","x_roofit", -1.0, 1.0)
    xframe = x_roofit.frame(ROOT.RooFit.Title("Some title here?"))
    fit_stats_dict      = RooFit_iterative_gaus_fit(h,  x_roofit,    xframe, iters=iters, num_sigmas=num_sigmas, 
                                                binned_data=do_binned_fit, draw_stats=draw_stats, 
                                                only_draw_last=only_draw_last, 
                                                line_color=color_before_corr, marker_color=color_before_corr+2)
    fit_stats_dict_corr = RooFit_iterative_gaus_fit(h_corr, x_roofit, xframe, iters=iters, num_sigmas=num_sigmas, 
                                                binned_data=do_binned_fit, draw_stats=draw_stats, 
                                                only_draw_last=only_draw_last, 
                                                line_color=color_after_corr, marker_color=color_after_corr+2)
    
    h_corr.SetTitle(title_str)
    h_corr.SetXTitle(dpToverpT_str)
    h_corr.SetYTitle("Events / [%.4f]" % h.GetBinWidth(0) )
    
    bestfit_mean = fit_stats_dict["mean_ls"][-1]
    bestfit_mean_err = fit_stats_dict["mean_err_ls"][-1]
    bestfit_std = fit_stats_dict["std_ls"][-1]
    bestfit_std_err = fit_stats_dict["std_err_ls"][-1]

    bestfit_mean_corr = fit_stats_dict_corr["mean_ls"][-1]
    bestfit_mean_err_corr = fit_stats_dict_corr["mean_err_ls"][-1]
    bestfit_std_corr = fit_stats_dict_corr["std_ls"][-1]
    bestfit_std_err_corr = fit_stats_dict_corr["std_err_ls"][-1]

    latex = ROOT.TLatex()
    latex_corr = ROOT.TLatex()
    latex.SetNDC()
    latex_corr.SetNDC()
    latex.SetTextSize(0.013) 
    latex_corr.SetTextSize(0.013) 
    def shift_y(y, go_down):
        return y - go_down
    x_coord = 0.185
    y_coord = 0.90
    go_down = 0.03
    latex.SetTextColor(color_before_corr)
    latex.DrawText(x_coord, y_coord, "Gaus fit iterations = %d" % iters)
    new_y = shift_y(y_coord, go_down+0.02)
    latex.DrawText(x_coord, new_y, "Before pT corrections:")
    new_y = shift_y(new_y, go_down)
    latex.DrawText(x_coord, new_y, "mean  = %.5E +- %.5E" % (bestfit_mean, bestfit_mean_err))
    new_y = shift_y(new_y, go_down)
    latex.DrawText(x_coord, new_y, "sigma = %.5E +- %.5E" % (bestfit_std, bestfit_std_err))
    new_y = shift_y(new_y, go_down+0.02)
    
    latex_corr.SetTextColor(color_after_corr)
    latex_corr.DrawText(x_coord, new_y, "After pT corrections:")
    new_y = shift_y(new_y, go_down)
    latex_corr.DrawText(x_coord, new_y, "mean  = %.5E +- %.5E" % (bestfit_mean_corr, bestfit_mean_err_corr))
    new_y = shift_y(new_y, go_down)
    latex_corr.DrawText(x_coord, new_y, "sigma = %.5E +- %.5E" % (bestfit_std_corr, bestfit_std_err_corr))
    

    dict_of_fit_stats_dict[h.GetName()] = fit_stats_dict
    dict_of_fit_stats_dict[h_corr.GetName()] = fit_stats_dict_corr

    # Save some info:
    total_entries += h.GetEntries()
    total_entries_corr += h_corr.GetEntries()
    # Make one page in PDF.
    canv.Print(fullpath_pdf_)
# End loop over HistDict obj, analyzed in parallel.
            # End loop over pT_bin_ls.
            # End loop over eta_bin_ls.
print(f"...closing big PDF")
canv.Print(fullpath_pdf_ + "]")

tot_mu = 0
for samp in individ_sample_ls:
    tot_mu += sample_dict[samp]["tree"].GetEntries() * 2.0
print(f"  Total muons in files:   {tot_mu}")
print(f"  Total muons passed sel: {total_entries}")
print(f"  Total muons corrected:  {total_entries_corr}")
with open(fullpath_pkl, "wb") as f:
    pickle.dump(dict_of_fit_stats_dict, f, protocol=2)
print(f"dict of fit_stats_dict saved at:\n{fullpath_pkl}")