"""
# Purpose: 
#   This script makes one PDF and one root file.
#   The PDF shows the pT distribution of a m2mu sample, 
#   before and after pT correction.
#   Then the 2017 pT correction factors are applied per muon.
#   Finally the corrected pT distribution is made and compared
#   to the uncorrected dist.
#   The root file contains a TTree ("passed_m2mu_selection")
#   and 2 histograms: one before pT correction and one after. 
# Notes: 
#   Give the file path to the correct pickled dict of pT correction factors.
#   Example key of pT_corr_dict = "2.3eta2.4_100.0pT150.0".
#   Only works on 1 sample_name at a time.
#   FIXME:
# Syntax : python <script>.py
# Author : Jake Rosenzweig
# Edited : 2020-07-19
"""
# NOTE: Possibly make this code ALSO perform the m2mu BWxCB fit.
import ROOT as r
import os
import pickle
from Utils_Python.Utils_Files import makeDirs, check_overwrite
from Utils_ROOT.ROOT_Plotting import Root_Hist_GetLastBinRightEdge
from d0_Studies.d0_Utils.d0_fns import calc_num_bins, correct_muon_pT
from Samples.sample_info import sample_info_dict
#--- User Params ---#
year = "2017"
sample_name = "DY"
# "Jpsi"  
# "DY"    
# "Upsilon
infile_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_{year}/MC_{year}_DY.root"
inpkl_path = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/MC{year}_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0.pkl"

outplot_dir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/pT"
outfile_dir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/pT"

filename_base = f"MC{year}DY_pTdist_applycorr_test09"

n_evts = 10000
overwrite = False
verbose = False

color_before_corr = r.kBlack
color_after_corr = r.kGreen + 2

bin_info_pT_dist = [0.0, 200.0, 1.0]  # [bin_min, bin_max, bin_width]
bin_info_deltapT_dist = [-1.0, 1.0, 0.01]

eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
pT_binedge_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]
#--- Script Functions ---#
def prep_area(outplot_dir, filename, overwrite):
    """Create PDF name, and see if the PDF and output dir exist."""
    fullpath_pdf = os.path.join(outplot_dir, filename + ".pdf")
    fullpath_rootfile = fullpath_pdf.replace(".pdf", ".root")
    check_overwrite(fullpath_pdf, overwrite)  # See if files already exist.
    makeDirs(outplot_dir)  # Make directories, if need be.
    return fullpath_rootfile, fullpath_pdf

def open_pT_corr_factor_dict(inpkl_path):
    """Return the dictionary that holds the pT correction factors."""
    with open(inpkl_path, "rb") as inpkl:
        d = pickle.load(inpkl)
    return d

def make_empty_hists(bin_info_pT_dist, bin_info_deltapT_dist):
    """Create 2 empty histograms: before pT correction, after pT correction."""
    binedge_min = bin_info_pT_dist[0]
    binedge_min_dpT = bin_info_deltapT_dist[0]
    binedge_max = bin_info_pT_dist[1]
    binedge_max_dpT = bin_info_deltapT_dist[1]
    bin_width = bin_info_pT_dist[2]
    bin_width_dpT = bin_info_deltapT_dist[2]
    n_bins = calc_num_bins(*bin_info_pT_dist)
    n_bins_dpT = calc_num_bins(*bin_info_deltapT_dist)

    title = "d_{0} corrections applied to p_{T} of #mu^{#pm} from DY ({%s} MC)" % year
    title_diff = "Shift in p_{T} (#equiv p_{T}^{reco} - p_{T}^{corr.})"
    x_title = "p_{T} [GeV]"
    y_title = f"Events / ({bin_width:.1f} GeV)"
    y_title_diff = f"Events / ({bin_width_dpT:.4f} GeV)"

    h_pT_no_corr = r.TH1F("h_pT_no_corr", f"{title};{x_title};{y_title}", 
                     n_bins, binedge_min, binedge_max)
    h_pT_corr = r.TH1F("h_pT_corr", f"{title};{x_title};{y_title}", 
                     n_bins, binedge_min, binedge_max)
    h_pT_diff = r.TH1F("h_pT_diff", f"{title_diff};{x_title};{y_title_diff}", 
                     n_bins_dpT, binedge_min_dpT, binedge_max_dpT)

    h_pT_no_corr.Sumw2()
    h_pT_corr.Sumw2()
    h_pT_diff.Sumw2()
    return (h_pT_no_corr, h_pT_corr, h_pT_diff)

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
                      h_pT_no_corr, h_pT_corr, h_pT_diff, 
                      eta_binedge_ls, pT_binedge_ls, pT_corr_factor_dict):
    inv_mass_cut_lim = sample_info_dict[sample_name]["inv_mass_cut_lim"]
    for count,evt in enumerate(tree):
        # if count >= n_evts: break
        if (count % 500000) == 0: print(f"Running over event: {count}")
        if passed_2mu_selections(evt, inv_mass_cut_lim):
            mu1 = initialize_muon(evt, "1")
            mu2 = initialize_muon(evt, "2")
            for mu in [mu1, mu2]:
                pT_corr = correct_muon_pT(mu.Eta(), mu.Pt(), mu.charge, mu.d0,
                                          pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls,
                                          verbose=verbose)
                h_pT_no_corr.Fill(mu.Pt())
                h_pT_corr.Fill(pT_corr)
                h_pT_diff.Fill(mu.Pt() - pT_corr)

def passed_2mu_selections(evt, inv_mass_cut_lim):
    """
    Return True, if mu1 and mu2 satisfy m2mu_min < massZ_FSR < m2mu_max GeV.
    m2mu_min/_max come from sample_info_dict[sample_name]["inv_mass_cut_lim"].
    """
    m2mu_min = inv_mass_cut_lim[0]
    m2mu_max = inv_mass_cut_lim[1]
    return True if (m2mu_min < evt.massZ_FSR) and (evt.massZ_FSR < m2mu_max) else False
        
def prep_canvas_and_pads():
    """Create the canvas, an upper pad (hists), and a lower pad (ratio plot)."""
    c1 = r.TCanvas("c1", "c1", 800, 600)
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
    return (c1, ptop, pbot)

def make_ratio_leg_line(h_numer, h_denom, text_size=0.10):
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
    ratio.GetYaxis().SetTitleOffset(0.3);  # Default is 0.005.
    ratio.SetXTitle("p_{T} [GeV]")
    ratio.SetYTitle("corr. / uncorr.")
    ratio.SetMinimum(0.)
    ratio.SetMaximum(2.)
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

def draw_everything(ptop, pbot, h_pT_no_corr, h_pT_corr, h_ratio, 
               leg, line, color_before_corr, color_after_corr):
    """
    Draw 2 histograms to the top pad and ratio plot to the bottom pad.
    Return a ratio line to make sure it doesn't disappear.
    """
    r.gROOT.SetBatch(r.kTRUE)
    ptop.cd()
    # c1.SetTicks(1,1)
    r.gStyle.SetOptStat("iouRMe")
    h_pT_no_corr.Draw("hist")
    h_pT_no_corr.SetLineColor(color_before_corr)
    h_pT_corr.Draw("same hist")
    h_pT_corr.SetLineColor(color_after_corr)
    # Stretch the y-axis to show all bin heights.
    max1 = h_pT_no_corr.GetMaximum()
    max2 = h_pT_corr.GetMaximum()
    y_max = max(max1, max2) * 1.1
    h_pT_corr.GetYaxis().SetRangeUser(0.0, y_max)
    
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
    fullpath_rootfile, fullpath_pdf = prep_area(outplot_dir, filename_base, overwrite)
    outf = r.TFile(outpath_rootfile, "recreate")
    newtree = r.TTree("tree", "m2mu_info")

    ptr_m2mu = array('f', [0.])
    ptr_m2mu_corr = array('f', [0.])
    pT_FSR = array('f', [0.])
    pT_FSR_corr = array('f', [0.])
    newtree.Branch("pT_FSR", pT_FSR, "pT_FSR/F")
    newtree.Branch("pT_FSR_corr", pT_FSR_corr, "pT_FSR_corr/F")
    newtree.Branch("m2mu", ptr_m2mu, "m2mu/F")
    newtree.Branch("m2mu_corr", ptr_m2mu_corr, "m2mu_corr/F")

    print(f"...Opening pickle:\n{inpkl_path}")
    pT_corr_factor_dict = open_pT_corr_factor_dict(inpkl_path)
    print("...Making empty histograms.")
    h_pT_no_corr, h_pT_corr, h_pT_diff = make_empty_hists(bin_info_pT_dist,
                                                          bin_info_deltapT_dist)
    print("...Filling histograms and correcting pTs.")
    fill_hist_pT_vals(sample_info_dict, sample_name, t,
                h_pT_no_corr, h_pT_corr, h_pT_diff, 
                eta_binedge_ls, pT_binedge_ls, pT_corr_factor_dict)

    
    ptr_m2mu[0] = m2mu
    ptr_m2mu_corr[0] = m2mu
    newtree.Fill()
    print("Filling done.\nPreparing canvas and pads.")
    c1, ptop, pbot = prep_canvas_and_pads()
    
    h_ratio, leg, line = make_ratio_leg_line(h_pT_corr, h_pT_no_corr)

    c1.Print(fullpath_pdf + "[")
    draw_everything(ptop, pbot, h_pT_no_corr, h_pT_corr, h_ratio, 
                    leg, line, color_before_corr, color_after_corr)
    # line.Draw("same")
    c1.Print(fullpath_pdf)
    c1.cd(0)
    c1.Clear()
    h_pT_diff.Draw("hist")
    c1.Print(fullpath_pdf)
    c1.Print(fullpath_pdf + "]")

    outf.cd()
    print(f"Saving hists and TTree to file:\n{outpath_rootfile}")
    h_m4mu.Write()
    h_m4mu_corr.Write()
    h_m4mu_diff.Write()
    newtree.Write()

    outf.Close()
    f.Close()


    # all_evts = t.GetEntries()
    # if n_evts == -1:
    #     n_evts = all_evts



    
    # graph = r.TGraph(good_evts_adhoc, array('f', m4mu_ls), array('f', m4mu_corr_ls))
    # graph.GetXaxis().SetTitle("m_{4#mu} [GeV]")
    # graph.GetYaxis().SetTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
    # graph.SetTitle("Effect of p_{T} corrections from d_{0} studies on m_{4#mu}")
    # graph.SetMarkerColor(r.kBlue)

    # print(f"Drawing histograms to:\n{outpath_pdf}")
    # c1 = r.TCanvas()
    # c1.SetTicks(1,1)
    # r.gStyle.SetOptStat("iouRMe")
    # c1.Print(outpath_pdf + "[")
    # h_m4mu.SetXTitle("m_{4#mu} [GeV]")
    # h_m4mu.SetYTitle("Events")
    # fix_stats_box_and_draw(h_m4mu, c1, dim=1)
    # c1.Print(outpath_pdf)
    # c1.Print(outpath_pdf + "]")

    # outf.cd()
    # print(f"Saving hists and TTree to file:\n{outpath_rootfile}")
    # h_m4mu.Write()
    # h_m4mu_corr.Write()
    # h_m4mu_diff.Write()
    # newtree.Write()

    # outf.Close()
    # f.Close()
    # print(f"Found {good_evts_adhoc} good m4mu events after selections.")







# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
# #@@@@@ OLD M2MU PLOTTER @@@@@#
# #@@@@@@@@@@@@@@@@@@@@@@@@@@@@#

# #----------------------#
# #----- Automatons -----#
# #----------------------#
# # Only need to configure this if you introduce new samples,
# # or change file paths. 


# "DY"      : {"LaTeX_name" : r"\mathrm{Z}", 
#                  "inv_mass_cut_lim" : [60.0, 120.0],
#                  "dR_cut" : 0.002}
# # Quick checks.
# if "combined" in sample_ls:
#     # Must also specify at least one sample (like, "Jpsi"). 
#     assert len(sample_ls) > 1
# if (apply_d0_pT_corrections):
#     msg = "You have to have a hist_type with 'corr' in it."
#     assert any(["corr" in x for x in hist_type_ls]), msg
#     # Load up the correction factors. 
#     with open(inpkl_path, "rb") as f:
#         pT_corr_factor_dict = pickle.load(f)

# # Make directories, if need be.
# makeDirs(outplot_dir)
# makeDirs(outpkl_dir)

# # Handle naming of files. 
# full_extra = (
#     f"_{num_sigmas:.2f}sigmas"
#     f"__{min(eta_ls):.1f}_eta_{max(eta_ls):.1f}"
#     f"__{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
# )
# full_extra = make_str_title_friendly(full_extra)

# fullfilename = f"{filename}_{year}_{full_extra}"

# fullpath_pkl = os.path.join(outpkl_dir, fullfilename + ".pkl")
# for samp in sample_ls:
#     # Check if PDF for this sample exists.
#     filename_samp = f"{fullfilename}_{samp}.pdf"
#     fullpath_file_samp = os.path.join(outplot_dir, filename_samp)
#     check_overwrite(fullpath_file_samp, overwrite)

# # See if pickle already exists.
# # Could check for PDFs, but trickier since there are many
# # produced each time this script is run.
# check_overwrite(fullpath_pkl, overwrite)

# # Everything seems good. Start the analysis.

# # Plotting info.
# if not (show_plots_to_screen): gROOT.SetBatch(kTRUE)
# tdrStyle = setTDRStyle()
# tdrGrid(tdrStyle, gridOn=True)

# dpToverpT_str = "(p_{T}^{REC} - p_{T}^{GEN}) / p_{T}^{GEN}"

# #---------------------------#
# #----- Local Functions -----#
# #---------------------------#
# class HistDict:
#     def __init__(self, hist_dict, sample, hist_type):
#         self.hist_dict = hist_dict
#         self.sample = sample
#         self.hist_type = hist_type
        
#     def get_TH1F_from_dict(self, eta_bin_ls, pT_min, pT_max):
#         eta_min = eta_bin_ls[0]
#         eta_max = eta_bin_ls[1]
#         pT_min = pT_bin_ls[0]
#         pT_max = pT_bin_ls[1]

#         key = f"h_{eta_min}eta{eta_max}_{pT_min}pT{pT_max}_{self.sample}_{self.hist_type}"
#         h = self.hist_dict[key]
#         assert isinstance(h, r.TH1F)
#         return h

# class OrgHistDict:
#     def __init__(self):
#         self.dict_of_HistDicts = {}
    
#     def store_HistDict(self, HistDict_):
#         sample = HistDict_.sample
#         hist_type = HistDict_.hist_type
#         key = f"{sample}_{hist_type}"
#         self.dict_of_HistDicts[key] = HistDict_ 

#     def get_HistDict(self, sample, hist_type):
#         key = f"{sample}_{hist_type}"
#         return self.dict_of_HistDicts[key]

#     def get_kinbin_TH1F(self, eta_bin_ls, pT_bin_ls, sample, hist_type):
#         etamin = eta_bin_ls[0]
#         etamax = eta_bin_ls[1]
#         pTmin = pT_bin_ls[0]
#         pTmax = pT_bin_ls[1]

#         hdict = self.get_HistDict(sample, hist_type).hist_dict
#         key = f"h_{etamin}eta{etamax}_{pTmin}pT{pTmax}_{sample}_{hist_type}"
#         return hdict[key]

#     def make_ls_TH1F_samekinbin_diffsample(self, eta_bin_ls, pT_bin_ls, sample_ls, hist_type):
#         assert "combined" not in sample_ls
#         # Get all TH1F that are in this kinbin for all individual samples.
#         th1f_ls = []
#         for sample in sample_ls:
#             th1f_ls.append( self.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, sample, hist_type) )
#             # key = f"{sample}_{hist_type}"
#             # hdict = self.get_HistDict(sample, hist_type)[key]
#         return th1f_ls

#     def add_diffsamples_kinbin_TH1F(self, h_comb, h_ls):
#         """
#         Add all TH1F within one kinbin into single, combined TH1F.
#         """
#         for hist in h_ls:
#             h_comb.Add(hist)

#     def fill_combined_kinbin_hist(self, eta_bin_ls, pT_bin_ls, sample_ls, hist_type):
#         """
#         Based on given parameters, find the corresponding combined TH1F 
#         (from the "combined" hist_dict) in this organizer, and fill the TH1F 
#         using all entries from individual samples.
#         """
#         th1f_ls = self.make_ls_TH1F_samekinbin_diffsample(eta_bin_ls, pT_bin_ls, sample_ls, hist_type)
#         # Get the empty TH1F that corresponds to all the given parameters.
#         h_comb = self.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, "combined", hist_type)
#         assert h_comb.GetEntries() == 0

#         # Fill the empty, combined hist.
#         self.add_diffsamples_kinbin_TH1F(h_comb, th1f_ls)

# def make_iter_binedge_ls(ls): 
#     """
#     Turns a list from this: [1,2,3,4]
#     into this: 
#       [ [1,2],
#         [2,3],
#         [3,4] ]

#     Good for iterative over neighboring elements.

#     ...and I just realized I can simply use a `zip()` while looping:
#     for eta_min, eta_max in zip(myls[:-1], myls[1:]): 
#         print(eta_min, eta_max) 
#     """
#     ls_beg = ls[:-1] 
#     ls_end = ls[1:] 
#     combine = zip(ls_beg, ls_end) 
#     return list(combine) 

# def get_hist_name_parts(h_name):
#     parts = h_name.split("_")
#     eta_piece = parts[1]
#     pT_piece = parts[2]
#     sample = parts[3]
#     h_type = parts[4:]

#     eta_min = float(eta_piece.split("eta")[0])
#     eta_max = float(eta_piece.split("eta")[1])
#     pT_min = float(pT_piece.split("pT")[0])
#     pT_max = float(pT_piece.split("pT")[1])
#     return eta_min, eta_max, pT_min, pT_max, sample, h_type
# #----------------#
# #----- Main -----#
# #----------------#
# # Initialize all hist dicts.
# org_hist_dict = OrgHistDict()
# eta_bin_ls_2D = make_iter_binedge_ls(eta_ls)
# pT_bin_ls_2D = make_iter_binedge_ls(pT_ls)

# for h_type in hist_type_ls:
#     for sample in sample_ls:
#         print(f"#----- Making {sample} hist_dict: {h_type} -----#")
#         h_dict = make_hist_lookuptable(eta_ls=eta_ls,pT_ls=pT_ls,sample=sample,
#                                        hist_type=h_type, bin_info_ls=bin_info_pT_dist)
#         print("hist_dict made.")
        
#         # Store them in an organizer.
#         org_hist_dict.store_HistDict(HistDict(hist_dict=h_dict, 
#                                                 hist_type=h_type, 
#                                                 sample=sample) 
#                                     )

# # Fill each sample hist dict. 
# individ_sample_ls = [x for x in sample_ls if "combined" not in x]
# for sample in individ_sample_ls:
#     hist_dict_ls = []
#     for hist_type in hist_type_ls:
#         hist_dict_ls.append( org_hist_dict.get_HistDict(sample, hist_type).hist_dict )

#     # Now enter each root file ONLY ONCE (maximize efficiency).
#     # Fill the sample hist and corresponding combined hist in parallel.
#     print(f"#----- Filling {sample} hist_dict: {hist_type_ls} -----#")
#     fill_dict_of_dpToverpT_hists(
#         tree=sample_dict[sample]["tree"], 
#         hist_dict=hist_dict_ls, 
#         hist_type=hist_type_ls,
#         sample=sample, 
#         n_evts=n_evts,
#         eta_binedge_ls=eta_ls, pT_binedge_ls=pT_ls,
#         apply_dR_cut=apply_dR_cut, apply_m2l_cut=apply_m2l_cut, verbose=verbose,
#         apply_d0_pT_corrections=apply_d0_pT_corrections, 
#         pT_corr_factor_dict=pT_corr_factor_dict
#     )
#     print("Hist dict filled.")

# # Each (hist_type, sample) has its own dict of TH1F. 
# # Only samples are filled so far, not the combined TH1Fs.
# for h_type in hist_type_ls:
#     for eta_bin_ls in eta_bin_ls_2D:
#         for pT_bin_ls in pT_bin_ls_2D:
#             org_hist_dict.fill_combined_kinbin_hist(eta_bin_ls, pT_bin_ls, individ_sample_ls, h_type)
            
#             # While we're here, check that number of events is consistent
#             # in the combined hists. 
#             num_entries = 0
#             for sample in individ_sample_ls:
#                 num_entries += org_hist_dict.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, sample, h_type).GetEntries()
#             assert num_entries == org_hist_dict.get_kinbin_TH1F(eta_bin_ls, pT_bin_ls, "combined", h_type).GetEntries()

# # Combined hists are all filled. 

# dict_of_fit_stats_dict = {}
# total_entries = 0
# total_entries_corr = 0

# # for sample in sample_ls:
#     # Create PDF name.
# filename_sample = f"{fullfilename}.pdf"  # fullfilename is entire filename, except extension.
# fullpath_pdf_ = os.path.join(outplot_dir, filename_sample)
# # Make a PDF for each sample (including "combined"!)
# HD      = org_hist_dict.get_HistDict("combined", "dpToverpT")
# HD_corr = org_hist_dict.get_HistDict("combined", "dpToverpTcorr")

# canv = TCanvas(filename, filename, 600, 600)
# print(f"Opening up big PDF...")
# canv.Print(fullpath_pdf_ + "[")

# # Plot similar histograms.
# # for eta_bin_ls in eta_bin_ls_2D:
# #     for pT_bin_ls in pT_bin_ls_2D:
# for h,h_corr in zip(HD.hist_dict.values(), HD_corr.hist_dict.values()):
#     # h = hdictobj.get_TH1F_from_dict(eta_bin_ls, pT_bin_ls)
#     # h_corr = hdictobj.get_TH1F_from_dict(eta_bin_ls, pT_bin_ls)

#     # Quick checks.
#     h_name_parts = get_hist_name_parts(h.GetName())
#     h_corr_name_parts = get_hist_name_parts(h_corr.GetName())
#     for part,part_corr in zip(h_name_parts[:-1],h_corr_name_parts[:-1]):
#         # Make sure everything agrees except the hist_type. 
#         assert part == part_corr

#     # Pretty up the plots.
#     latex_str = sample_dict[HD.sample]["fancy_font"]["ROOT"]
#     title_str = (
#         f"MC {year} {latex_str}    binning=("
#         f"{h_name_parts[0]:.1f} < #left|#eta#right| < {h_name_parts[1]:.1f}    "
#         )
#     title_str += r"%.1f < p_{T} < %.1f GeV), " % (h_name_parts[2], h_name_parts[3])
#     title_str += "fit_range=#mu#pm{:.2f}#sigma".format(num_sigmas)
    
#     print(f"Information about these hists:\n"
#             f"  name={h.GetName()},       name_corr={h_corr.GetName()}"
#             f"  entries={h.GetEntries()}, entries_corr={h_corr.GetEntries()}\n")
    
#     print(f"...Performing iterative gaus fit on hists...")
#     canv.cd()

#     x_roofit = r.RooRealVar("x_roofit","x_roofit", -1.0, 1.0)
#     xframe = x_roofit.frame(r.RooFit.Title("Some title here?"))
#     fit_stats_dict      = RooFit_iterative_gaus_fit(h,  x_roofit,    xframe, iters=iters, num_sigmas=num_sigmas, 
#                                                 binned_data=do_binned_fit, draw_stats=draw_stats, 
#                                                 only_draw_last=only_draw_last, 
#                                                 line_color=color_before_corr, marker_color=color_before_corr+2)
#     fit_stats_dict_corr = RooFit_iterative_gaus_fit(h_corr, x_roofit, xframe, iters=iters, num_sigmas=num_sigmas, 
#                                                 binned_data=do_binned_fit, draw_stats=draw_stats, 
#                                                 only_draw_last=only_draw_last, 
#                                                 line_color=color_after_corr, marker_color=color_after_corr+2)
    
#     h_corr.SetTitle(title_str)
#     h_corr.SetXTitle(dpToverpT_str)
#     h_corr.SetYTitle("Events / [%.4f]" % h.GetBinWidth(0) )
    
#     bestfit_mean = fit_stats_dict["mean_ls"][-1]
#     bestfit_mean_err = fit_stats_dict["mean_err_ls"][-1]
#     bestfit_std = fit_stats_dict["std_ls"][-1]
#     bestfit_std_err = fit_stats_dict["std_err_ls"][-1]

#     bestfit_mean_corr = fit_stats_dict_corr["mean_ls"][-1]
#     bestfit_mean_err_corr = fit_stats_dict_corr["mean_err_ls"][-1]
#     bestfit_std_corr = fit_stats_dict_corr["std_ls"][-1]
#     bestfit_std_err_corr = fit_stats_dict_corr["std_err_ls"][-1]

#     latex = r.TLatex()
#     latex_corr = r.TLatex()
#     latex.SetNDC()
#     latex_corr.SetNDC()
#     latex.SetTextSize(0.013) 
#     latex_corr.SetTextSize(0.013) 
#     def shift_y(y, go_down):
#         return y - go_down
#     x_coord = 0.185
#     y_coord = 0.90
#     go_down = 0.03
#     latex.SetTextColor(color_before_corr)
#     latex.DrawText(x_coord, y_coord, "Gaus fit iterations = %d" % iters)
#     new_y = shift_y(y_coord, go_down+0.02)
#     latex.DrawText(x_coord, new_y, "Before pT corrections:")
#     new_y = shift_y(new_y, go_down)
#     latex.DrawText(x_coord, new_y, "mean  = %.5E +- %.5E" % (bestfit_mean, bestfit_mean_err))
#     new_y = shift_y(new_y, go_down)
#     latex.DrawText(x_coord, new_y, "sigma = %.5E +- %.5E" % (bestfit_std, bestfit_std_err))
#     new_y = shift_y(new_y, go_down+0.02)
    
#     latex_corr.SetTextColor(color_after_corr)
#     latex_corr.DrawText(x_coord, new_y, "After pT corrections:")
#     new_y = shift_y(new_y, go_down)
#     latex_corr.DrawText(x_coord, new_y, "mean  = %.5E +- %.5E" % (bestfit_mean_corr, bestfit_mean_err_corr))
#     new_y = shift_y(new_y, go_down)
#     latex_corr.DrawText(x_coord, new_y, "sigma = %.5E +- %.5E" % (bestfit_std_corr, bestfit_std_err_corr))
    

#     dict_of_fit_stats_dict[h.GetName()] = fit_stats_dict
#     dict_of_fit_stats_dict[h_corr.GetName()] = fit_stats_dict_corr

#     # Save some info:
#     total_entries += h.GetEntries()
#     total_entries_corr += h_corr.GetEntries()
#     # Make one page in PDF.
#     canv.Print(fullpath_pdf_)
# # End loop over HistDict obj, analyzed in parallel.
#             # End loop over pT_bin_ls.
#             # End loop over eta_bin_ls.
# print(f"...closing big PDF")
# canv.Print(fullpath_pdf_ + "]")

# tot_mu = 0
# for samp in individ_sample_ls:
#     tot_mu += sample_dict[samp]["tree"].GetEntries() * 2.0
# print(f"  Total muons in files:   {tot_mu}")
# print(f"  Total muons passed sel: {total_entries}")
# print(f"  Total muons corrected:  {total_entries_corr}")
# with open(fullpath_pkl, "wb") as f:
#     pickle.dump(dict_of_fit_stats_dict, f, protocol=2)
# print(f"dict of fit_stats_dict saved at:\n{fullpath_pkl}")