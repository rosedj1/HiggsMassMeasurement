from __future__ import print_function
import ROOT
import os
import sys
from array import array
from pprint import pprint
from collections import OrderedDict
import numpy as np
# Local imports.
sys.path.append("/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/")
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid, fixOverlay
from Utils_Python.Utils_Files import check_overwrite
from Utils_ROOT.ROOT_classes import make_TH1F, make_TGraphErrors, make_TLegend
from Utils_ROOT.ROOT_fns import read_cpp_as_txt, load_cpp_code
from Utils_ROOT.ROOT_StatsAndFits import DSCBFitter, DSCBFitScanner, DSCBFitPlotter

# infile_path = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/MC2018ggH_passFull_fullstats.root"
infile_path = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/MC2018ggH_passFull_fullstats.root"
outfile_dir = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/DSCB_fit_m4mu/"
# outfile_dir = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/plots/test/DSCBscanoutput/"
outfile_prefix = "MC2018ggH_applypTcorrAdHoc_test01"
cpp_DSCB_code_path = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/fit_and_draw_DSCB.C"
# cpp_DSCB_code_path = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/fit_and_draw_DSCB.C"

overwrite = 0
m4mu_min, m4mu_max = 105.0, 140.0
relmasserr_min, relmasserr_max = 0.0, 0.75
m4muErr_binedge_ls = [0.0, 0.75, 1.25, 5.0]  # This is a percent!
# m4muErr_binedge_ls = [0.0, 0.75]
max_lower_edge = 116.5
min_upper_edge = 128.0
step_GeV = 0.5
n_bins = 100
# m4mu_min = 121.2
# m4mu_max = 128.4

def make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, relmasserr_min, relmasserr_max, step_GeV):
    f = "{}_{}m4mu{}_step{:.1f}GeV_{}relmasserr{}".format(outfile_prefix, m4mu_min, m4mu_max, step_GeV, relmasserr_min, relmasserr_max)
    outfile_name = f.replace(".", "p") + ".txt"
    return os.path.join(outfile_dir, outfile_name)

#***** (Probably) useless functions below *****#
def get_params_and_errs(tree, m4mu_min, m4mu_max):
    """Return mean and mean_err from DSCB on m4mu range.
    
    FIXME: Makes a call to deprecated fn: ROOT.DSCB_unbinnedfit
    """
    print("Performing DSCB unbinned fit.")
    result = ROOT.DSCB_unbinnedfit(t, m4mu_min, m4mu_max)
    result_tup = parse_DSCBfit_result(result)
    return result_tup

def make_fit_dict(fit_range_edges, max_lower_edge, min_upper_edge):
    """Return a dict of fit stats from unbinned DSCB...
    
    FIXME: Makes a call to deprecated fn: ROOT.DSCB_unbinnedfit
    """
    fit_dct = OrderedDict()
    # Use the lowest vales in fit_range_edges and
    # the highest vals in its flipped array.
    for minval, maxval in zip(fit_range_edges, fit_range_edges[::-1]):
        if (minval >= max_lower_edge) or (maxval < min_upper_edge):
            break
        print("\n...Filling dict over range: [{}, {}]".format(minval, maxval))
        # mean, mean_err, sigma, sigma_err = get_params_and_errs(t, minval, maxval)
        key = (minval, maxval)
        fit_dct[key] = {"mean":mean, "sigma":sigma, "mean_err":mean_err, "sigma_err":sigma_err}
    return fit_dct

def dscb_fits_GeVsteps(infile_path, outfile_dir, outfile_prefix, m4mu_min, m4mu_max, step_GeV, max_lower_edge, bins=100, overwrite=False):
    """
    Make a txt file with the fit stats from multiple DSCB fits on m4mu dist.
    Subsequent fits will narrow the fit range by `step_GeV`.
    """
    # outfile_txt = make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, step_GeV)
    # check_overwrite(outfile_txt, overwrite)
    # print("Creating file:\n", outfile_txt)
    # print("Loading C++ DSCB unbinned fitting code...")
    # load_cpp_code()
    # print("Opening file:\n", infile_path)
    # f = ROOT.TFile(infile_path)
    # t = f.Get("passedEvents")

    bins = (m4mu_max - m4mu_min) / float(step_GeV)
    fit_range_edges = np.linspace(m4mu_min, m4mu_max, bins+1)
    print("Using step_size = {} GeV".format(step_GeV))
    print("Using bins:\n", fit_range_edges)
    # lower_bounds = np.linspace(105, 120, 4)
    # upper_bounds = np.linspace(130, 140, 4)[::-1]  # Flip it around to zip properly.
    fit_dct = make_fit_dict(t, fit_range_edges, max_lower_edge)
    pprint(fit_dct)
    write_output_file(outfile_txt, fit_dct)
    # End loop over bounds.
    #print("\nFinal means and errors:")
    #pprint(fit_list)
#***** (Probably) useless functions above *****#
def write_output_file(outfile_txt, fit_dct):
    with open(outfile_txt, "w") as f:
        header = "fit_range_min\tfit_range_max\tmean\tmean_err\tsigma\tsigma_err\n"
        f.write(header)
        for key,val in fit_dct.items():
            tup_of_keys = (key[0], key[1], val["mean"], val["mean_err"], val["sigma"], val["sigma_err"])
            print("...Writing these vals to file:", tup_of_keys)
            line = "{:13.8f}\t{:13.8f}\t{:13.8f}\t{:13.8f}\t{:13.8f}\t{:13.8f}\n".format(*tup_of_keys)
            f.write(line)

def make_and_fill_hists(t):
    """
    Return a list of 3 filled histograms:
    - mass4mu
    - mass4lErr
    - mass4mu/mass4lErr * 100%
    """
    # Make histograms.
    h_m4mu = make_TH1F(r"h_m4mu", title=r"Mass", xlabel=r"m_{4#mu}", n_bins=140, x_min=105, x_max=140, units=r"GeV")
    h_dm4mu = make_TH1F(r"h_dm4mu", title=r"Mass Error", xlabel=r"#deltam_{4#mu}", n_bins=100, x_min=0, x_max=5, units=r"GeV")
    h_reldm4mu = make_TH1F(r"h_reldm4mu", title=r"Relative Mass Error", xlabel=r"#deltam_{4#mu}/m_{4#mu} (%)", n_bins=100, x_min=0, x_max=5.0)
    cuts = "passedFullSelection==1 && finalState==1"
    # Fill histograms.
    t.Draw("mass4l >> h_m4mu", cuts, "goff")
    t.Draw("mass4lErr >> h_dm4mu", cuts, "goff")
    t.Draw("mass4lErr/mass4l*100 >> h_reldm4mu", cuts, "goff")
    return [h_m4mu, h_dm4mu, h_reldm4mu]

def make_pdf_mass4mu_mass4lErr_relmass4lErr(tree):
    """Create a PDF with 3 dists: mass4mu, mass4lErr, their ratio."""
    print("...Making and filling hists.")
    hist_ls = make_and_fill_hists(tree)
    print("...Hists made. Making plots.")

    canv = ROOT.TCanvas()
    canv.Print(outfile + "[")
    for h in hist_ls:
        h.Draw("hist")
        canv.Print(outfile)
    canv.Print(outfile + "]")
    print("Done.")

if __name__ == "__main__":
    print("...Preparing your area.")
    outfile_pdf = make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, relmasserr_min, relmasserr_max, step_GeV)
    outfile_pdf = outfile_pdf.replace(".txt", ".pdf")
    check_overwrite(outfile_pdf, overwrite=overwrite)
    tdrStyle = setTDRStyle(show_statsbox=True)
    tdrGrid(tdrStyle, gridOn=False)
    ROOT.gROOT.SetBatch(True)

    print("...Loading DSCB fitter code...")
    load_cpp_code("./RooMyPDF_DSCB_C.so", cpp_DSCB_code_path)
    print("...Opening root file...")
    f = ROOT.TFile(infile_path)
    tree = f.Get("passedEvents")

    print("Using step_size = {} GeV".format(step_GeV))
    # print("Using bins:\n", fit_range_edges)
    print("...Preparing canvas.")
    canv = ROOT.TCanvas()
    canv.Print(outfile_pdf + "[")

    print("Performing DSCB fits.")
    scanner = DSCBFitScanner(relmasserr_min, relmasserr_max)  # Loop over these.
    scanner.do_DSCB_fits_over_mass4mu_range(m4mu_min, m4mu_max, step_GeV,
                                                max_lower_edge, min_upper_edge,
                                                tree, canv, outfile_pdf, n_bins)
    plotter = DSCBFitPlotter()
    plotter.plot_X_vs_GeVfitrange("mean", scanner, canv, outfile_pdf, make_new_page_after=True)
    plotter.plot_X_vs_GeVfitrange("sigma", scanner, canv, outfile_pdf, make_new_page_after=True)
    plotter.plot_X_vs_GeVfitrange("integral", scanner, canv, outfile_pdf, make_new_page_after=True)
    # plotter.plot_X_vs_GeVfitrange("sigma", scanner, canv, outfile_pdf, make_new_page_after=True)
    canv.Print(outfile_pdf + "]")
    print("Done.")