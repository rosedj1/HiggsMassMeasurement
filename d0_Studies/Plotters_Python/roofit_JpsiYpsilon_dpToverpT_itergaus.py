import os
import pickle
import ROOT
from ROOT import TFile, TH1F, TCanvas, gROOT, kTRUE

from Utils_Python.Utils_Files import make_dirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Selections import Selector
from Utils_Python.Utils_Physics import calc_dphi, calc_dR
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import fixOverlay, setTDRStyle, tdrGrid

from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit 
from Utils_ROOT.ROOT_fns import fill_dpToverpT_hist

# from d0_Studies.d0_Utils.d0_dicts import color_dict_RooFit
from d0_Studies.d0_Utils.d0_fns import calc_num_bins
#-----------------------#
#----- User Params -----#
#-----------------------#
infile_path_Jpsi_dpToverpT = "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Jpsi.root"
infile_path_Upsilon_dpToverpT = "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Upsilon.root"

outdir_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/Testing/"
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/hists_dpToverpT/"

filename = "20200602_fullstats_dpToverpT"
overwrite = False
verbose = False

n_evts = -1

do_iter_gaus_fit = True
iters = 6
num_sigmas = 1.75

do_BWconvCB_fit = False

sample_ls = ["Jpsi", "Upsilon"]

x_limits_dpToverpT = [-0.08, 0.08, 0.0004]  # [bin_min, bin_max, bin_width]

eta_ls = [0.0, 2.4]
pT_ls  = [5.0, 200.0]

show_plots_to_screen = False
apply_dR_cut = True
apply_m2l_cut = True

dpToverpT_str = "(p_{T}^{REC} - p_{T}^{GEN}) / p_{T}^{GEN}"
#----------------------#
#----- Automatons -----#
#----------------------#
# Sanity checks. 
assert do_iter_gaus_fit != do_BWconvCB_fit

# Make directories, if need be.
make_dirs(outdir_plots)
make_dirs(outdir_pkl)

# Handle naming of files. 
extra = (
    f"_{num_sigmas:.2f}sigmas"
    f"__{min(eta_ls)}_eta_{max(eta_ls)}"
    f"__{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
)
extra = make_str_title_friendly(extra)

fullfilename = filename + extra

fullpath_pdf = os.path.join(outdir_plots, fullfilename + ".pdf")
fullpath_pkl = os.path.join(outdir_pkl, fullfilename + ".pkl")

# See if files already exist.
check_overwrite(fullpath_pdf, overwrite)
check_overwrite(fullpath_pdf, overwrite)

if not (show_plots_to_screen): gROOT.SetBatch(kTRUE)
tdrStyle = setTDRStyle()
tdrGrid(tdrStyle, gridOn=True)

# Access data. 
f_Jpsi    = TFile.Open(infile_path_Jpsi_dpToverpT)
f_Upsilon = TFile.Open(infile_path_Upsilon_dpToverpT)

t_Jpsi    = f_Jpsi.Get("passedEvents")
t_Upsilon = f_Upsilon.Get("passedEvents")

n_tot_Jpsi = t_Jpsi.GetEntries()
n_tot_Upsilon = t_Upsilon.GetEntries()

x_min_dpToverpT = x_limits_dpToverpT[0]
x_max_dpToverpT = x_limits_dpToverpT[1]
n_bins_dpToverpT = calc_num_bins(*x_limits_dpToverpT)
#----------------#
#----- Main -----#
#----------------#
h_Jpsi_dpToverpT = TH1F("h_Jpsi_dpToverpT", "J/#psi mass resonance", n_bins_dpToverpT, x_min_dpToverpT, x_max_dpToverpT)
h_Upsilon_dpToverpT = TH1F("h_Upsilon_dpToverpT", "#Upsilon mass resonance", n_bins_dpToverpT, x_min_dpToverpT, x_max_dpToverpT)

# Pretty up the plots.
h_Jpsi_dpToverpT.SetTitle("MC 2017 J/#psi #Deltap_{T}/p_{T}, fit range = #mu#pm%.1f#sigma"%num_sigmas)
h_Jpsi_dpToverpT.SetXTitle(dpToverpT_str)
h_Jpsi_dpToverpT.SetYTitle("Events / [%.4f GeV]" % h_Jpsi_dpToverpT.GetBinWidth(0) )
h_Upsilon_dpToverpT.SetTitle("MC 2017 #Upsilon #Deltap_{T}/p_{T}, fit range = #mu#pm%.1f#sigma"%num_sigmas)
h_Upsilon_dpToverpT.SetXTitle(dpToverpT_str)
h_Upsilon_dpToverpT.SetYTitle("Events / [%.4f GeV]" % h_Upsilon_dpToverpT.GetBinWidth(0) )

canv = TCanvas(filename, filename, 600, 600)

if eta_ls is not None:
    check_eta = True
    eta_min = eta_ls[0]
    eta_max = eta_ls[1]
if pT_ls is not None:
    check_pT = True
    pT_min = pT_ls[0]
    pT_max = pT_ls[1]

# Select events and fill histograms.
print(f"Running over ")
fill_dpToverpT_hist(t_Jpsi, h_Jpsi_dpToverpT,
                "Jpsi", n_evts,
                bounds_ls_eta=eta_ls, bounds_ls_pT=pT_ls,
                apply_dR_cut=apply_dR_cut, apply_m2l_cut=apply_m2l_cut)
fill_dpToverpT_hist(t_Upsilon, h_Upsilon_dpToverpT, 
                "Upsilon", n_evts,
                bounds_ls_eta=eta_ls, bounds_ls_pT=pT_ls,
                apply_dR_cut=apply_dR_cut, apply_m2l_cut=apply_m2l_cut)
print(f"Opening up big PDF...")
canv.Print(fullpath_pdf + "[")
dict_of_fit_stats_dict = {}

for h in (h_Jpsi_dpToverpT, h_Upsilon_dpToverpT):
    canv.cd()
    h.Sumw2()

    print(
        f"Information about this hist: {h.GetName()}\n"
        f"  h.GetEntries: {h.GetEntries()}\n\n"
        )
    
    print(f"Performing iterative gaus fit on hist: {h.GetName()}")
    if (do_iter_gaus_fit):
        # FIXME: Probably need to draw onto the canvas here!
        fit_stats_dict = RooFit_iterative_gaus_fit(h, 
                            binned_fit=True, switch_to_binned_fit=2000, 
                            iters=iters, num_sigmas=num_sigmas, 
                            n_bins=100, x_lim=None, 
                            xframe=None, x_label="Independent var", title="", units="",
                            marker_color=1, only_draw_last=False, draw=False, verbose=verbose)
    elif (do_BWconvCB_fit):
        # fit_stats_dict = RooFit_BWconvCB_fit(h, canv, iters=iters, num_sigmas=num_sigmas, binned_fit=True)
        pass
    
    dict_of_fit_stats_dict[h.GetName()] = fit_stats_dict

    # Make one page in PDF.
    canv.Print(fullpath_pdf)
    
print(f"...closing big PDF")
canv.Print(fullpath_pdf + "]")

with open(fullpath_pkl, "wb") as f:
    pickle.dump(dict_of_fit_stats_dict, f, protocol=2)
print(f"dict of fit_stats_dict saved at:\n{fullpath_pkl}")