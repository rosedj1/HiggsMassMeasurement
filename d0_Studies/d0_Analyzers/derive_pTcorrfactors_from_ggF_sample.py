"""Ad Hoc pT Correction Analyzer
Purpose:
    Select "good" H->ZZ->4mu events,
    put each muon into User-specified (eta, pT) bins,
    make all kinematic plots into a single PDF.
Syntax:  
    python script.py > output.txt  (recommended)
    python script.py 
Notes:   
    This code runs on a Higgs sample.
    Make sure to check all the parameters in "User Parameters".
    Should be used with Python 3.X.

    This code makes the following distributions, before/after pT corr:
        (1) muon pT dist.
        (2) m4mu dist. (DSCB fit?)
        (3) eta dist., for each eta bin
        (4) muon qd0 dist, per (eta, pT) bin
    The following plots are produced:
Author:  Jake Rosenzweig
Created: 2020-08-24
Updated: 2020-09-02
"""
import ROOT as r
import numpy as np
# Local imports.
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid
from d0_Studies.KinBin_Info.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum)
from d0_Studies.d0_Analyzers.Synch_with_Filippo.save_m4mu_info_andcorr_ggF_synchwithFilippo import (
    make_outfiles, passed_Higgs_selections, get_ndcs_gen, make_muon_ls, verify_all_muons, calc_Hmass
)
#----- User Parameters -----#
year = 2018
inpath_file = f"/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_{year}.root"
# inpath_pkl = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/MC{year}_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0.pkl"
outpath_file_woext = f"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/MC{year}_ggF_derivecorrfromHiggsSample_test03"

eta_binedge_ls = equal_entry_bin_edges_eta_mod1_wholenum
pT_binedge_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum

max_n_evts = 5000
overwrite = False
verbose = True
#----- Functions -----#
def build_muons_from_event(t, evt_num):
    """Look at event in Higgs sample: if muons pass selections, save muons.
    
    Parameters
    ----------
    t : ROOT.TTree
        The TTree which holds all data for each event.
    evt_num : int
        Which event in the TTree.

    Returns
    -------
    If event passes selections: 
        4-tuple of ROOT.Math.LorentzVector objects.
    If event does NOT pass selections: 
        4-tuple of NoneType (None, None, None, None).
    """
    global n_evts_passed
    bad_muons = (None, None, None, None)
    t.GetEntry(evt_num)
    
    if not passed_Higgs_selections(t):
        return bad_muons

    lep_Hindex_ls = list(t.lep_Hindex)  # Elements that correspond to 4 leptons which build Higgs candidate.
    rec_ndcs_ls = lep_Hindex_ls
    lep_genindex_ls = list(t.lep_genindex)
    # if not validate_lep_genindex(lep_genindex_ls):
    #     continue
    # We have a good lep_genindex: at least 4 leptons have been matched. 
    # Could be 5 or more leps.
    gen_ndcs_ls = get_ndcs_gen(rec_ndcs_ls, lep_genindex_ls)
    # Now gen_ndcs_ls should be the same length as rec_ndcs_ls:
    assert len(rec_ndcs_ls) == len(gen_ndcs_ls)

    rec_ID_ls = list(t.lep_id)
    gen_ID_ls = list(t.GENlep_id)
    # if not check_matched_IDs(rec_ndcs_ls, gen_ndcs_ls, rec_ID_ls, gen_ID_ls):
    #     continue
    # if not check_2_OSSF_muon_pairs(gen_ID_ls):
    #     continue
    # if not check_2_OSSF_muon_pairs(rec_ID_ls):
    #     continue

    # Event looks good so far. 
    # Now check kinematics of muons.
    mu_ls = make_muon_ls(t, rec_ndcs_ls, gen_ndcs_ls)
    all_muons_passed = verify_all_muons(mu_ls)
    if not all_muons_passed:
        return bad_muons

    # NOW event is good.
    n_evts_passed += 1

    # Return the 4 good muons.
    good_muons = tuple(mu_ls)
    assert len(good_muons) == 4
    return good_muons

def set_plot_styles(draw=False, gridOn=True):
    """Make your plots consistently pretty."""
    import ROOT as r
    r.gROOT.SetBatch(r.kTRUE)
    # Plotting info.
    if not (draw): r.gROOT.SetBatch(True)
    tdrStyle = setTDRStyle()
    tdrGrid(tdrStyle, gridOn=gridOn)

class MuonCollection:
    """A class to handle a list of ROOT.Math.LorentzVector objects."""
    
    def __init__(self):
        self.muon_ls = []
        self.m4mu_ls = []

    def extract_muons_from_Higgs_file(self, inpath_file, n_evts=-1):
        """Set self.muon_ls of ROOT.Math.LorentzVector muons which pass Higgs selections.
        
        Parameters
        ----------
        inpath_file : str
            Absolute file path to root file.
        n_evts : int
            Number of events to scan over.
            Default of `-1` runs over all events.
        """
        print(f"...Opening root file:\n{inpath_file}")
        f = r.TFile(inpath_file, "read")
        t = f.Get("Ana/passedEvents")

        global max_n_evts
        global n_evts_passed
        all_evts = t.GetEntries()
        if max_n_evts == -1:
            max_n_evts = all_evts
        print(f"...Running over {max_n_evts} events...")

        n_evts_passed = 0
        # Event loop.
        for evt_num in range(max_n_evts):
            if evt_num >= max_n_evts:
                break
            if evt_num % 20000 == 0:
                print(f"Running over evt: {evt_num}")
            # Check if this event passes Higgs selections.
            # Build the 4 muons from the event.
            mu_tup = build_muons_from_event(t, evt_num)
            if None in mu_tup: 
                continue
            # Add the good muons to the final muon list.
            self.muon_ls.extend(mu_tup)
            # Save the m4mu info.
            m4mu = calc_Hmass(*mu_tup)
            self.m4mu_ls.append(m4mu)
            # Correct the pT of each muon. 
            # mu_corr_ls = make_muon_corr_ls(mu_ls, pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls, verbose=verbose)

            # m4mu = calc_Hmass(*mu_tup)
            # m4mu_corr = calc_Hmass(*mu_corr_ls)
            # m4mu_diff = m4mu_corr - m4mu
            
            # m4mu_ls.append(float(m4mu))
            # m4mu_corr_ls.append(float(m4mu_corr))

            # in_bounds = (105 < m4mu) and (m4mu < 140)
            # in_bounds_corr = (105 < m4mu_corr) and (m4mu_corr < 140)
            # if (not in_bounds) or (not in_bounds_corr):
            #     print(f"Event={evt_num}, m4mu={m4mu}, m4mu_corr={m4mu_corr}, nFSRPhotons={t.nFSRPhotons}")
            #     show_diff_array("pT", list(t.lepFSR_pt), list(t.lep_pt))
            #     show_diff_array("eta", list(t.lepFSR_eta), list(t.lep_eta))
            #     show_diff_array("phi", list(t.lepFSR_phi), list(t.lep_phi))
            #     show_diff_array("mass", list(t.lepFSR_mass), list(t.lep_mass))
            #     print_muon_ls_info(mu_tup, mu_corr_ls)

            # rel_diff = m4mu_diff / float(m4mu)
            # if abs(rel_diff) > 0.05:
            #     print(f"event {evt_num}:")
            #     print(f"  m4mu={m4mu}, m4mu_corr={m4mu_corr}, rel_diff={rel_diff}")
            #     print_muon_ls_info(mu_ls, mu_corr_ls)

            # ptr_m4mu[0] = m4mu
            # ptr_m4mu_corr[0] = m4mu_corr
            # newtree.Fill()

            # h_m4mu.Fill(m4mu)
            # h_m4mu_diff.Fill(m4mu_diff)
            # h_m4mu_corr.Fill(m4mu_corr)
            # h_m4muvsm4mucorr.Fill(m4mu, m4mu_corr)
            # for mu_,mucorr_ in zip(mu_ls, mu_corr_ls):
            #     shift = mu_.Pt() - mucorr_.Pt()
            #     h_deltapT_corr.Fill(shift)
        # End evt loop.
        assert len(self.muon_ls) / len(self.m4mu_ls) == 4

    def make_kinematic_plots(self):
        """Run over self.muon_ls and self.m4mu_ls to make and fill all plots."""
        print(f"Creating histograms...")
        h_m4mu = r.TH1F("h_m4mu", "m_{4#mu} distribution, Higgs Production Mode: ggF", 140, 105, 140)
        # h_m4mu.SetXTitle("m_{4#mu} (GeV)")
        h_pT = r.TH1F("h_pT", "muon p_{T}^{reco} distribution", 220, 0, 220)
        h_pT_gen = r.TH1F("h_pT_gen", "muon p_{T}^{gen} distribution", 220, 0, 220)
        h_eta = r.TH1F("h_eta", "muon #eta^{reco} distribution", 100, -2.5, 2.5)
        h_eta_gen = r.TH1F("h_eta_gen", "muon #eta^{gen} distribution", 100, -2.5, 2.5)
        h_phi = r.TH1F("h_phi", "muon #phi^{reco} distribution", 50, -4, 4)
        h_phi_gen = r.TH1F("h_phi_gen", "muon #phi^{gen} distribution", 50, -4, 4)
        h_d0 = r.TH1F("h_d0", "muon d_{0} distribution", 100, -0.01, 0.01)
        h_charge = r.TH1F("h_charge", "muon charge distribution", 100, -10, 10)
        h_dpTOverpT = r.TH1F("h_dpTOverpT", "#Deltap_{T}/p_{T} distribution", 500, -0.6, 0.6)

        print(f"Filling histograms...")
        for m4mu in self.m4mu_ls:
            h_m4mu.Fill(m4mu)
        for muon in self.muon_ls:
            h_pT.Fill(muon.Pt())
            h_eta.Fill(muon.Eta())
            h_phi.Fill(muon.Phi())
            h_d0.Fill(muon.d0)
            h_charge.Fill(muon.charge)
            h_dpTOverpT.Fill(muon.dpTOverpT)
        # h_m4mu_corr = r.TH1F("h_m4mu_corr", "m_{4#mu}^{p_{T},corr.} distribution (Higgs Production Mode: ggF)", 140, 105, 140)
        # h_m4mu_diff = r.TH1F("h_m4mu_diff", "#Deltam_{4#mu} #equiv m_{4#mu}^{corr. p_{T}} - m_{4#mu}", 200, -20, 20)
        # h_m4muvsm4mucorr = r.TH2F("h_m4muvsm4mucorr", "Correlation between m_{4#mu}^{corr. p_{T}} and m_{4#mu}", 
        #                              100, 70, 170, 100, 70, 170)  # n_binX, X_Low,X_Hig, n_binY, Y_low, Y_high
        # h_deltapT_corr = r.TH1F("h_deltapT_corr", "p_{T} correction", 200, -5, 5)
        self.plot_ls = [h_m4mu, h_pT, h_eta, h_phi, h_d0, h_charge, h_dpTOverpT]

    def make_pdf_of_plots(self, outpath_pdf):
        """Draw all plots in self.plot_ls to a TCanvas."""
        c = r.TCanvas()
        c.Print(outpath_pdf + "[")
        for plot in self.plot_ls:
            plot.Draw("hist")
            c.Print(outpath_pdf)
        c.Print(outpath_pdf + "]")

def main():
    """When doing `python this_script.py`, run this code."""
    # Prep your area.
    set_plot_styles()
    outpath_rootfile, outpath_pdf = make_outfiles(outpath_file_woext, overwrite)
    
    # Begin analysis.
    muon_collection = MuonCollection()
    muon_collection.extract_muons_from_Higgs_file(inpath_file, n_evts=max_n_evts)

    # Make and print plots.
    muon_collection.make_kinematic_plots()
    muon_collection.make_pdf_of_plots(outpath_pdf)


    # muon_collection.derive_correction_factors(eta_ls, pT_ls)
if __name__ == "__main__":
    main()
    
    
    




    

    # with open(inpath_pkl, "rb") as pkl:
    #     print(f"...Opening pickle:\n{inpath_pkl}")
    #     pT_corr_factor_dict = pickle.load(pkl)

    # New file, TTree, and TH1Fs.
#     outf = r.TFile(outpath_rootfile, "recreate")
#     newtree = r.TTree("tree", "tree_m4mu_vals")

    
    
#     h_m4mu.Sumw2()
#     h_m4mu_corr.Sumw2()
#     h_m4mu_diff.Sumw2()
#     h_m4muvsm4mucorr.Sumw2()
#     h_deltapT_corr.Sumw2()

#     ptr_m4mu = array('f', [0.])
#     ptr_m4mu_corr = array('f', [0.])
#     ptr_m4mu_diff = array('f', [0.])
#     newtree.Branch("m4mu", ptr_m4mu, "m4mu/F")
#     newtree.Branch("m4mu_corr", ptr_m4mu_corr, "m4mu_corr/F")

    








    
#     graph = r.TGraph(n_evts_passed, array('f', m4mu_ls), array('f', m4mu_corr_ls))
#     graph.GetXaxis().SetTitle("m_{4#mu} [GeV]")
#     graph.GetYaxis().SetTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
#     graph.SetTitle("Effect of p_{T} corrections from d_{0} studies on m_{4#mu}")
#     graph.SetMarkerColor(r.kBlue)

#     def fix_stats_box_and_draw(hist, canv, dim=1):
#         if dim == 1:
#             hist.Draw("hist")
#         elif dim == 2:
#             hist.Draw("colz")
#         r.gPad.Update()
#         statsbox = hist.FindObject("stats")
#         statsbox.SetX1NDC(0.75)
#         statsbox.SetX2NDC(0.90)
#         statsbox.SetY1NDC(0.75)
#         statsbox.SetY2NDC(0.90)
#         if dim == 1:
#             hist.Draw("hist")
#         elif dim == 2:
#             hist.Draw("colz")
#         canv.Update()

#     print(f"Drawing histograms to:\n{outpath_pdf}")
#     c1 = r.TCanvas()
#     c1.SetTicks(1,1)
#     r.gStyle.SetOptStat("iouRMe")
#     c1.Print(outpath_pdf + "[")
#     h_m4mu.SetXTitle("m_{4#mu} [GeV]")
#     h_m4mu.SetYTitle("Events")
#     fix_stats_box_and_draw(h_m4mu, c1, dim=1)
#     # h_m4mu.Draw("hist")
#     c1.Print(outpath_pdf)
#     h_m4mu_corr.SetXTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
#     h_m4mu_corr.SetYTitle("Events")
#     fix_stats_box_and_draw(h_m4mu_corr, c1, dim=1)
#     c1.Print(outpath_pdf)
#     h_m4mu_diff.SetXTitle("#Deltam_{4#mu} [GeV]")
#     h_m4mu_diff.SetYTitle("Events")
#     fix_stats_box_and_draw(h_m4mu_diff, c1, dim=1)
#     c1.Print(outpath_pdf)
#     h_m4muvsm4mucorr.SetXTitle("m_{4#mu} [GeV]")
#     h_m4muvsm4mucorr.SetYTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
#     fix_stats_box_and_draw(h_m4muvsm4mucorr, c1, dim=2)
#     c1.Print(outpath_pdf)
#     h_deltapT_corr.SetXTitle("#Deltap_{T} [GeV]")
#     h_deltapT_corr.SetYTitle("Events")
#     fix_stats_box_and_draw(h_deltapT_corr, c1, dim=1)
#     c1.Print(outpath_pdf)
#     c1.Print(outpath_pdf + "]")

#     outf.cd()
#     print(f"Saving hists and TTree to file:\n{outpath_rootfile}")
#     h_m4mu.Write()
#     h_m4mu_corr.Write()
#     h_m4mu_diff.Write()
#     newtree.Write()

#     outf.Close()
#     f.Close()
#     print(f"Found {n_evts_passed} good m4mu events after selections.")


    











# # OLD d0 equal-entry region code

# from Utils_vaex.vaex_fns import vaex_apply_masks, prepare_vaex_df
# from Samples.sample_info import Sample
# from d0_Utils.d0_fns import find_equal_hist_regions_unbinned
# from d0_Studies.KinBin_Info.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum,
#                                                    bin_edges_pT_sevenfifths_to1000GeV_wholenum)
# from Utils_Python.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite      

# # def ParseOption():
# #     parser = argparse.ArgumentParser(description='submit all')
# #     parser.add_argument('--pklfilename', dest='filename_base', type=str, help='') 
# #     parser.add_argument('--verbose', dest='verbose', type=int, default=1, help='')
# #     parser.add_argument('--overwrite', dest='overwrite', type=int, default=0, help='')  
    
# #     args = parser.parse_args()                                                                                         
# #     return args          
                                                                                                         
# # args = ParseOption()                    
                                                                     
# #---------------------------#
# #----- User Parameters -----#
# #---------------------------#
# # Where to save pkl and csv files.
# outdir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/"

# # filename_base = args.filename_base
# # overwrite = args.overwrite
# # verbose = args.verbose
# makeDirs(outdir)

# if (verbose):
#     print(f"|eta| regions:\n  {np.round(eta_ls, decimals=2)}")
#     print(f"pT regions:\n  {np.round(pT_ls, decimals=2)}\n")
# #----------------#
# #----- Main -----#
# #----------------#
# extra = (
#     f"_{r}reg"
#     f"with{algo[1]}perreg"
#     f"__{min(eta_ls):.1f}_eta_{max(eta_ls):.1f}"
#     f"__{min(pT_ls):.1f}_pT_{max(pT_ls):.1f}_GeV"
# )
# extra = make_str_title_friendly(extra)
# extra += ".csv"

# fullpath_csv = os.path.join(outdir, filename_base + extra)
# fullpath_pickle = fullpath_csv.replace(".csv", ".pkl")


# equal_entry_binedge_dict = {
#     "all_eta_bins" : eta_ls,
#     "all_pT_bins" : pT_ls,
#     "sample_name_ls" : ["DY", "Jpsi", "DY+Jpsi"],
# }

# with open(fullpath_csv, "w") as myfile:
#     myfile.write(f"all_eta_bins : {eta_ls}\n")
#     myfile.write(f"all_pT_bins  : {pT_ls}\n\n")

#     # Loop over eta regions.
#     for k in range(len(eta_ls)-1):
#         eta_min = eta_ls[k]
#         eta_max = eta_ls[k+1]
#         eta_range = [eta_min, eta_max]
#         eta_key = f"eta_bin={eta_range}"

#         equal_entry_binedge_dict[eta_key] = {}

#         # Column names.
#         col_str  = (
#             f"{'sample':<7}, {'bins_found':<10}, {'bins_wanted':<11}, {'muons_per_reg':<13}, "
#             f"\teta_range,\t\tpT_range,\t\tq*d0_bins\n"
#         )
#         myfile.write(col_str)

        

            

#             # Append to dict.
#             equal_entry_binedge_dict[eta_key][pT_key]["equalentry_qd0ls_DY"] = equalentry_binedge_ls_DY
#             equal_entry_binedge_dict[eta_key][pT_key]["equalentry_qd0ls_Jpsi"] = equalentry_binedge_ls_Jpsi
#             equal_entry_binedge_dict[eta_key][pT_key]["equalentry_qd0ls_DY+Jpsi"] = equalentry_binedge_ls_comb

#             n_muons_DY = len(qd0_arr_sel_DY)
#             n_muons_Jpsi = len(qd0_arr_sel_Jpsi)
#             n_muons_comb = len(qd0_arr_sel_comb)
#             info  = (
#                 f"{'DY':<7}, {r_updated_DY:<10}, {r:<11}, {math.ceil(n_muons_DY / r_updated_DY):<13}, "
#                 f"\t{eta_range},\t\t{pT_range},\t\t{equalentry_binedge_ls_DY}\n"
#                 f"{'Jpsi':<7}, {r_updated_Jpsi:<10}, {r:<11}, {math.ceil(n_muons_Jpsi / r_updated_Jpsi):<13}, "
#                 f"\t{eta_range},\t\t{pT_range},\t\t{equalentry_binedge_ls_Jpsi}\n"
#                 f"{'DY+Jpsi':<7}, {r_updated_comb:<10}, {r:<11}, {math.ceil(n_muons_comb / r_updated_comb):<13}, "
#                 f"\t{eta_range},\t\t{pT_range},\t\t{equalentry_binedge_ls_comb}\n"
#             )
#             myfile.write(info)

#             total_entries += n_muons_DY + n_muons_Jpsi
#             print(f"total_entries up to this point: {total_entries}")
#         # End pT loop. Go to next eta range
#         myfile.write('\n')
#         print(f"Finished eta_range = {eta_range}\n\n")
#     # End eta loop.
    
# print(f"[INFO] q*d0 bin edge info written to csv file:\n{fullpath_csv}")

# with open(fullpath_pickle,'wb') as output:
#     pickle.dump(equal_entry_binedge_dict, output, protocol=2)
# print(f"[INFO] eta, pT, q*d0 bin dict written to pickle file:\n{fullpath_pickle}\n")

# total_muons_original = MC_2018_DY_hdf5.vdf_prepped.count() + MC_2018_Jpsi_hdf5.vdf_prepped.count()
# perc = total_entries / float(total_muons_original) * 100.
# print(
#     f"Total muons found before selections: {total_muons_original}\n"
#     f"Total muons found after selections (inclusive all regions): {total_entries} ({perc:.2f}%)"
# )