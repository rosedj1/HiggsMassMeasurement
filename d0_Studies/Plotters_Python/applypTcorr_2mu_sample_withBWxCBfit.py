"""
# Purpose: 
#   This script makes one PDF and one root file.
#   The PDF contains:
#   - distribution of pT of a m2mu sample before/after pT correction
#   - distribution of pT shifts due to corrections
#   - distribution of m2mu and BWxCB+Exp fits before/after pT correction
# Notes: 
#   Give the file path to the correct pickled dict of pT correction factors.
#   Example key of pT_corr_dict = "2.3eta2.4_100.0pT150.0".
#   The root file contains 2 TTrees and 5 hists: 
#     (tree1) tree_pT : contains all the pT values which passed selection
#     (tree2) tree_m2mu : contains all m2mu values which passed selection
#     -------
#     (hist1) h_pT_uncorr: uncorrected pT vals.
#     (hist1) h_pT_corr: correct pT vals.
#     (hist1) h_pT_diff: shifts in pT vals.
#     (hist1) h_m2mu: uncorrected m2mu vals.
#     (hist1) h_m2mu_corr: corrected m2mu vals.
#   Only works on 1 sample_name at a time.
# Syntax : python <script>.py
# Author : Jake Rosenzweig
# Edited : 2020-07-22
"""
import ROOT as r
import os
import pickle
from array import array
from Utils_Python.Utils_Files import make_dirs, check_overwrite
from Utils_ROOT.ROOT_Plotting import Root_Hist_GetLastBinRightEdge
from Utils_ROOT.ROOT_StatsAndFits import RooFit_CBxBWplusExp_fit_binned
from d0_Studies.d0_Utils.d0_fns import calc_num_bins, correct_muon_pT
from Samples.sample_info import sample_info_dict
#-------------------#
#--- User Params ---#
#-------------------#
year = "2018"
sample_name = "DY"
# "Jpsi"  
# "DY"    
# "Upsilon
infile_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_{year}/MC_{year}_DY.root"
inpkl_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/MC{year}_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0.pkl"

outplot_dir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/hists_m2mu/"
outfile_dir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/m2mu/"

filename_base = f"MC{year}DY_applycorr_m2mu_dists_fullstats"

n_evts = -1  # Use `-1` to run over all events.
overwrite = False
verbose = False

color_before_corr = r.kBlack
color_after_corr = r.kGreen + 2

bin_info_pT_dist = [0.0, 200.0, 1.0]  # [bin_min, bin_max, bin_width]
bin_info_deltapT_dist = [-5.0, 5.0, 0.025]
bin_info_m2mu_dist = [60.0, 120.0, 0.5]

eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
pT_binedge_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]
#------------------------#
#--- Script Functions ---#
#------------------------#
def prep_area(outpath_dir, filename, ext="", overwrite=False):
    """Create PDF name, and see if the PDF and output dir exist."""
    outpath_file = os.path.join(outpath_dir, f"{filename}{ext}")
    check_overwrite(outpath_file, overwrite)  # See if files already exist.
    make_dirs(outpath_dir)  # Make directories, if need be.
    return outpath_file

def open_pT_corr_factor_dict(inpkl_path):
    """Return the dictionary that holds the pT correction factors."""
    with open(inpkl_path, "rb") as inpkl:
        d = pickle.load(inpkl)
    return d

def make_empty_hists(bin_info_pT_dist, bin_info_deltapT_dist, bin_info_m2mu_dist):
    """Create 2 empty histograms: before pT correction, after pT correction."""
    binedge_min = bin_info_pT_dist[0]
    binedge_min_dpT = bin_info_deltapT_dist[0]
    binedge_min_m2mu = bin_info_m2mu_dist[0]
    binedge_max = bin_info_pT_dist[1]
    binedge_max_dpT = bin_info_deltapT_dist[1]
    binedge_max_m2mu = bin_info_m2mu_dist[1]
    bin_width = bin_info_pT_dist[2]
    bin_width_dpT = bin_info_deltapT_dist[2]
    bin_width_m2mu = bin_info_m2mu_dist[2]
    n_bins = calc_num_bins(*bin_info_pT_dist)
    n_bins_dpT = calc_num_bins(*bin_info_deltapT_dist)
    n_bins_m2mu = calc_num_bins(*bin_info_m2mu_dist)

    title = "p_{T} dist : d_{0} corrections applied to p_{T} of #mu^{#pm} from DY (%s MC)" % year
    title_diff = "Shift in p_{T} (#equiv p_{T}^{reco} - p_{T}^{corr.})"
    title_m2mu = "m_{2#mu} dist : d_{0} corrections applied to p_{T} of #mu^{#pm} from DY (%s MC)" % year
    x_title = "p_{T} [GeV]" 
    x_title_m2mu = "m_{2#mu} [GeV]"
    y_title = f"Events / ({bin_width:.1f} GeV)"
    y_title_diff = f"Events / ({bin_width_dpT:.4f} GeV)"
    y_title_m2mu = f"Events / ({bin_width_m2mu:.1f} GeV)"

    h_pT_uncorr = r.TH1F("h_pT_uncorr", f"{title};{x_title};{y_title}", 
                     n_bins, binedge_min, binedge_max)
    h_pT_corr = r.TH1F("h_pT_corr", f"{title};{x_title};{y_title}", 
                     n_bins, binedge_min, binedge_max)
    h_pT_diff = r.TH1F("h_pT_diff", f"{title_diff};{x_title};{y_title_diff}", 
                     n_bins_dpT, binedge_min_dpT, binedge_max_dpT)
    h_m2mu = r.TH1F("h_m2mu", f"{title_m2mu};{x_title_m2mu};{y_title_m2mu}", 
                     n_bins_m2mu, binedge_min_m2mu, binedge_max_m2mu)
    h_m2mu_corr = r.TH1F("h_m2mu_corr", f"{title_m2mu};{x_title_m2mu};{y_title_m2mu}", 
                     n_bins_m2mu, binedge_min_m2mu, binedge_max_m2mu)

    h_pT_uncorr.Sumw2()
    h_pT_corr.Sumw2()
    h_pT_diff.Sumw2()
    h_m2mu.Sumw2()
    h_m2mu_corr.Sumw2()
    return (h_pT_uncorr, h_pT_corr, h_pT_diff, h_m2mu, h_m2mu_corr)

def initialize_muon(evt, num):
    """
    For a given event (evt), build muon (mu) 1 or 2 
    using the reco values stored in the NTuple. 

    Additionally: 
    - checks ID to make sure it is a muon.
    - Stores: pT_FSR, eta_FSR, phi_FSR, mass_FSR, charge, d0BS 
    """
    mu = r.Math.PtEtaPhiMVector(getattr(evt, f"pT_FSR{num}"), 
                                getattr(evt, f"eta_FSR{num}"),
                                getattr(evt, f"phi_FSR{num}"),
                                getattr(evt, f"m_FSR{num}")
         )
    mu.charge = getattr(evt, f"Id{num}") / -13.
    if abs(mu.charge) != 1:
        print(f"mu.charge ({mu.charge}]) != +-1")
        raise ValueError
    mu.d0 = getattr(evt, f"d0BS{num}")
    return mu

def fill_hist_pT_vals(sample_info_dict, sample_name, tree,
                      h_pT_uncorr, h_pT_corr, h_pT_diff, h_m2mu, h_m2mu_corr,
                      pT_FSR, pT_FSR_corr, ptr_m2mu, ptr_m2mu_corr, tree_pT, tree_m2mu, 
                      eta_binedge_ls, pT_binedge_ls, pT_corr_factor_dict):
    inv_mass_cut_lim = sample_info_dict[sample_name]["inv_mass_cut_lim"]
    for count,evt in enumerate(tree):
        if (count >= n_evts) and (n_evts != -1): break
        if (count % 500000) == 0: print(f"Running over event: {count}")
        if passed_2mu_selections(evt, inv_mass_cut_lim):
            mu1 = initialize_muon(evt, "1")
            mu2 = initialize_muon(evt, "2")
            for mu in [mu1, mu2]:
                pT_uncorr = mu.Pt()
                pT_corr = correct_muon_pT(mu.Eta(), mu.Pt(), mu.charge, mu.d0,
                                          pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls,
                                          verbose=verbose)
                h_pT_uncorr.Fill(pT_uncorr)
                h_pT_corr.Fill(pT_corr)
                h_pT_diff.Fill(pT_uncorr - pT_corr)
                mu.SetPt(pT_corr)

                pT_FSR[0] = pT_uncorr
                pT_FSR_corr[0] = pT_corr
                tree_pT.Fill()
            # End loop over both muons.
            # Now calculate new m2mu.
            m2mu = evt.massZ_FSR
            h_m2mu.Fill(m2mu)
            Z_corr = mu1 + mu2
            m2mu_corr = Z_corr.M()
            h_m2mu_corr.Fill(m2mu_corr)
            ptr_m2mu[0] = m2mu
            ptr_m2mu_corr[0] = m2mu_corr
            tree_m2mu.Fill()

def passed_2mu_selections(evt, inv_mass_cut_lim):
    """
    Return True, if mu1 and mu2 satisfy m2mu_min < massZ_FSR < m2mu_max GeV.
    m2mu_min/_max come from sample_info_dict[sample_name]["inv_mass_cut_lim"].
    """
    m2mu_min = inv_mass_cut_lim[0]
    m2mu_max = inv_mass_cut_lim[1]
    return True if (m2mu_min < evt.massZ_FSR) and (evt.massZ_FSR < m2mu_max) else False
        
def make_pads():
    """Create an upper pad (for hists), and a lower pad (for a ratio plot)."""
    ptop = r.TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0)
    pbot = r.TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25)
    ptop.SetBottomMargin(0)
    pbot.SetTopMargin(0)
    pbot.SetBottomMargin(0.25)
    ptop.Draw()
    pbot.Draw()
    ptop.cd()
    ptop.SetTicks(1,1)
    pbot.cd()
    pbot.SetTicks(1,1)
    return (ptop, pbot)

def make_ratio_leg_line(h_numer, h_denom, x_title, text_size=0.10):
    """Return a filled ratio plot by dividing two histograms."""
    # Make the ratio plot.
    assert h_numer.GetNbinsX() == h_denom.GetNbinsX()
    assert h_numer.GetBinLowEdge(1) == h_denom.GetBinLowEdge(1)
    nbins = h_numer.GetNbinsX()
    x_min = h_numer.GetBinLowEdge(1)
    x_max = Root_Hist_GetLastBinRightEdge(h_numer)
    ratio = r.TH1F("h_ratio", "", nbins, x_min, x_max)
    # Fill the ratio plot.
    ratio.Divide(h_numer, h_denom)
    ratio.GetXaxis().SetLabelSize(text_size)
    ratio.GetYaxis().SetLabelSize(text_size)
    ratio.GetXaxis().SetTitleSize(text_size)
    ratio.GetYaxis().SetTitleSize(text_size)
    ratio.GetYaxis().SetTitleOffset(0.3)  # Default is 0.005.
    ratio.SetXTitle(x_title)
    ratio.SetYTitle("corr. / uncorr.")
    ratio.SetMinimum(0.8)
    ratio.SetMaximum(1.2)
    ratio.SetNdivisions(207, "Y")
    ratio.GetXaxis().SetTickLength(0.12)
    ratio.SetLineColor(color_after_corr)
    # Add reference line at y = 1.
    line = r.TLine(x_min, 1., x_max, 1.)
    # Make legend.
    leg = r.TLegend(0.75, 0.8, 0.9, 0.9)
    leg.AddEntry(h_denom, "Uncorr.", "le")
    leg.AddEntry(h_numer, "Corr.", "le")
    leg.Draw("same")
    return ratio, leg, line

def draw_2hists_and_ratio(ptop, pbot, h_no_corr, h_corr, h_ratio, 
               leg, line, color_before_corr, color_after_corr, tex=None):
    """Draw 2 histograms to the top pad and ratio plot to the bottom pad."""
    r.gROOT.SetBatch(r.kTRUE)
    ptop.cd()
    # c1.SetTicks(1,1)
    r.gStyle.SetOptStat("iouRMe")
    h_no_corr.Draw("hist")
    h_no_corr.SetLineColor(color_before_corr)
    h_corr.Draw("same hist")
    h_corr.SetLineColor(color_after_corr)
    # Stretch the y-axis to show all bin heights.
    max1 = h_no_corr.GetMaximum()
    max2 = h_corr.GetMaximum()
    y_max = max(max1, max2) * 1.12
    h_corr.GetYaxis().SetRangeUser(0.0, y_max)
    h_no_corr.GetYaxis().SetRangeUser(0.0, y_max)
    
    if tex is not None:
        tex.SetNDC()

        def get_stats_text(h_no_corr):
            """Return mean and stdev in text form for make-shift legend."""
            mean_no_corr = h_no_corr.GetMean()
            mean_no_corr_err = h_no_corr.GetMeanError()
            sigma_no_corr = h_no_corr.GetStdDev()
            sigma_no_corr_err = h_no_corr.GetStdDevError()
            text_no_corr = r"#splitline{Mean = %.2f #pm %.1E}" % (mean_no_corr, mean_no_corr_err)
            text_no_corr += r"{StdDev = %.2f #pm %.1E}" % (sigma_no_corr, sigma_no_corr_err)
            return text_no_corr
        text_no_corr = get_stats_text(h_no_corr)
        text_corr = get_stats_text(h_corr)

        tex.SetTextColor(color_before_corr)
        tex.DrawLatex(0.15, 0.8, text_no_corr)
        tex.SetTextColor(color_after_corr)
        tex.DrawLatex(0.15, 0.6, text_corr)

    leg.Draw("same")
    r.gPad.Update()

    # Draw ratio plot.
    pbot.cd()
    h_ratio.Draw()
    line.Draw("same")
    # c1.SetTicks(1,1)
    r.gStyle.SetOptStat(0)
    h_ratio.Draw("same")
    h_ratio.SetLineColor(color_after_corr)
    # statsbox = h_ratio.FindObject("stats")
    # statsbox.SetOptStat(0)
    r.gPad.Update()

# def fix_stats_box_and_draw(hist, canv, dim=1):
#     if dim == 1:
#         hist.Draw("hist")
#     elif dim == 2:
#         hist.Draw("colz")
#     r.gPad.Update()
#     statsbox = hist.FindObject("stats")
#     statsbox.SetX1NDC(0.75)
#     statsbox.SetX2NDC(0.90)
#     statsbox.SetY1NDC(0.75)
#     statsbox.SetY2NDC(0.90)
#     if dim == 1:
#         hist.Draw("hist")
#     elif dim == 2:
#         hist.Draw("colz")
#     canv.Update()

if __name__ == "__main__":
    print(f"Opening file:\n{infile_path}")
    f = r.TFile.Open(infile_path)
    t = f.Get("passedEvents")

    # New file, TTree, and TH1Fs.
    fullpath_pdf = prep_area(outplot_dir, filename_base, ext=".pdf", overwrite=overwrite)
    fullpath_rootfile = prep_area(outfile_dir, filename_base, ext=".root", overwrite=overwrite)
    outf = r.TFile(fullpath_rootfile, "recreate")
    tree_pT = r.TTree("tree_pT", "tree_pT")
    tree_m2mu = r.TTree("tree_m2mu", "tree_m2mu")

    ptr_m2mu = array('f', [0.])
    ptr_m2mu_corr = array('f', [0.])
    pT_FSR = array('f', [0.])
    pT_FSR_corr = array('f', [0.])
    tree_pT.Branch("pT_FSR", pT_FSR, "pT_FSR/F")
    tree_pT.Branch("pT_FSR_corr", pT_FSR_corr, "pT_FSR_corr/F")
    tree_m2mu.Branch("m2mu", ptr_m2mu, "m2mu/F")
    tree_m2mu.Branch("m2mu_corr", ptr_m2mu_corr, "m2mu_corr/F")

    print(f"...Opening pickle:\n{inpkl_path}")
    pT_corr_factor_dict = open_pT_corr_factor_dict(inpkl_path)
    print("...Making empty histograms.")
    h_pT_uncorr, h_pT_corr, h_pT_diff, h_m2mu, h_m2mu_corr = make_empty_hists(bin_info_pT_dist,
                                                                               bin_info_deltapT_dist,
                                                                               bin_info_m2mu_dist)
    print("...Filling histograms and correcting pTs.")
    fill_hist_pT_vals(sample_info_dict, sample_name, t,
                h_pT_uncorr, h_pT_corr, h_pT_diff, h_m2mu, h_m2mu_corr,
                pT_FSR, pT_FSR_corr, ptr_m2mu, ptr_m2mu_corr, tree_pT, tree_m2mu, 
                eta_binedge_ls, pT_binedge_ls, pT_corr_factor_dict)
    print("Filling done.\nPreparing canvas and pads...")
    c1 = r.TCanvas("c1", "c1", 800, 600)
    ptop, pbot = make_pads()
    h_ratio, leg, line = make_ratio_leg_line(h_pT_corr, h_pT_uncorr, x_title="p_{T} [GeV]")
    c1.Print(fullpath_pdf + "[")
    print("...drawing pT distributions.")
    draw_2hists_and_ratio(ptop, pbot, h_pT_uncorr, h_pT_corr, h_ratio, 
                    leg, line, color_before_corr, color_after_corr)
    # line.Draw("same")
    c1.Print(fullpath_pdf)
    # Start new page.
    c1.cd(0)
    c1.Clear()
    print("...drawing pT shift distribution.")
    h_pT_diff.Draw("hist")
    c1.Print(fullpath_pdf)
    print("...drawing BWxCB+exp distributions.")
    c1.Clear()
    # fit_hist_with_BWxCBplusExp(c1, h_m2mu, h_m2mu_corr)
    del ptop, pbot
    ptop, pbot = make_pads()
    h_ratio_m2mu, leg_m2mu, line_m2mu = make_ratio_leg_line(h_m2mu_corr, h_m2mu, x_title="m_{2#mu} [GeV]")
    # tex = r.TLatex()
    draw_2hists_and_ratio(ptop, pbot, h_m2mu, h_m2mu_corr, h_ratio_m2mu, 
                    leg_m2mu, line_m2mu, color_before_corr, color_after_corr, tex=None)
    # Draw the BWxCB
    xframe_uncorr, fit_stats_dict_uncorr = RooFit_CBxBWplusExp_fit_binned(h_m2mu, bin_info_m2mu_dist[:2], 
                                            fit_range=None, show_params=True, params_box=[0.13, 0.4, 0.7], 
                                            linecolor=color_before_corr, markercolor=color_before_corr)
    xframe_corr, fit_stats_dict_corr = RooFit_CBxBWplusExp_fit_binned(h_m2mu_corr, bin_info_m2mu_dist[:2], 
                                            fit_range=None, show_params=True, params_box=[0.63, 0.9, 0.7], 
                                            linecolor=color_after_corr, markercolor=color_after_corr)
    ptop.cd()
    xframe_uncorr.Draw("same")
    xframe_corr.Draw("same")
    tex = r.TLatex()
    tex.SetNDC()
    tex.SetTextColor(r.kRed)
    sigma_uncorr = fit_stats_dict_uncorr["sigma"]
    sigma_corr = fit_stats_dict_corr["sigma"]
    percdiff = abs(sigma_corr - sigma_uncorr) / sigma_uncorr * 100.0
    text_sigma_percdiff = r"#Delta#sigma/#sigma_{uncorr.} = %.1f" % percdiff + "%"
    tex.DrawLatex(0.15, 0.8, text_sigma_percdiff)
    tex.Draw("same")
    c1.Update()
    c1.Print(fullpath_pdf)
    # Close the PDF.
    c1.Print(fullpath_pdf + "]")

    outf.cd()
    print(f"Saving hists to file:\n{fullpath_rootfile}")
    h_pT_uncorr.Write()
    h_pT_corr.Write()
    h_pT_diff.Write()
    h_m2mu.Write()
    h_m2mu_corr.Write()
    tree_pT.Write()
    tree_m2mu.Write()

    outf.Close()
    f.Close()