"""Make and Fill KinBin2Ds for RC Studies.

Given a list of eta and pT bin edges, create a KB2D (eta, pT) for each
combination. Then a ggH->4mu sample is opened and muons are selected, based on
eta, pT, etc. Inclusive muon kinematic distributions are made.

The plots are saved in a root file.
The KB2Ds are saved in a pkl file.
 
The produced pkl file can then be used to do iterative Gaussian fits.
Use:
- roch_vs_noroch_itergausfit_template.py
- submit_to_slurm_inbatch.py
- roch_vs_noroch_slurm.sbatch

Author: Jake Rosenzweig
Updated: 2021-03-04
"""
import sys
import ROOT
import gc  # Garbage Collector. Good for cleaning up after loops.
import os
# Local imports.
from ParticleCollections import MyMuonCollection
from Utils_Python.Selections import build_muons_from_HZZ4mu_event
from Utils_Python.Utils_Files import check_overwrite
from Utils_ROOT.ROOT_classes import make_TH1F
from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid
from d0_Studies.d0_Utils.d0_cls import KinBin2D
from Utils_Python.Utils_Files import save_to_pkl
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum
from Utils_Python.Utils_Files import make_dirs
#----- User Switches -----#
overwrite = 1
verbose = 1
max_evts = -1
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum #[0.0, 0.9, 1.7, 2.4]
pT_ls = [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 75.0, 200.0]
# pT_ls = [5.0, 10.0, 15.0, 20.0, 30.0, 40.0, 50.0, 75.0, 100.0, 150.0, 200.0]
#pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum[:-1] #[5, 25, 50, 100]
inpath_file_RC = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2018.root"
inpath_file_noRC = "/cmsuf/data/store/user/t2/users/ferrico/SingleBS_studies/After/NoRoch/ggF_2018.root"

filename = "MC2018ggH_KB2D_onlymuoninfo_fullstats_pTbinsroundnumbers_75then200GeV"
outpath_file = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/plots/tests/{filename}.pdf"
outpath_pkl = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/tests/{filename}.pkl"
outpath_root = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/rootfiles/tests/{filename}_inclusivekinemhists.root"

def draw_Xhists_to_canvas(c, outpath_file, h_tup):
    """Draw 4 hists in h_ls to canvas divided into 2x2 pads."""
    X = len(h_tup)
    if X in (3,4):
        c.Divide(2,2)
    for h in range(X):
        c.cd(h + 1)
        h_tup[h].Draw("hist same")
    c.Print(outpath_file)
    c.Clear()

if __name__ == "__main__":
    tdrStyle = setTDRStyle(show_statsbox=True)
    tdrGrid(tdrStyle, gridOn=False)
    ROOT.gROOT.SetBatch(True)
    def prep_area(file_tup):
        for f in (outpath_file, outpath_root, outpath_pkl):
            dir_ = os.path.dirname(f)
            make_dirs(dir_)
            check_overwrite(f, overwrite)
    
    eta_min = min(eta_ls)
    eta_max = max(eta_ls)
    pT_min = min(pT_ls)
    pT_max = max(pT_ls)

    h_pT1 = make_TH1F("h_pT1", title=r"p_{T1,#mu} (leading)", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT2 = make_TH1F("h_pT2", title=r"p_{T2,#mu} (subleading)", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT3 = make_TH1F("h_pT3", title=r"p_{T3,#mu}", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT4 = make_TH1F("h_pT4", title=r"p_{T4,#mu}", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT1_RC = make_TH1F("h_pT1_RC", title=r"p^{RC}_{T1,#mu} (leading)", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT2_RC = make_TH1F("h_pT2_RC", title=r"p^{RC}_{T2,#mu} (subleading)", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT3_RC = make_TH1F("h_pT3_RC", title=r"p^{RC}_{T3,#mu}", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_pT4_RC = make_TH1F("h_pT4_RC", title=r"p^{RC}_{T4,#mu}", xlabel=r"p_{T}", n_bins=100, x_min=0, x_max=100, units=r"GeV")
    h_dpT_RC_reco = make_TH1F("h_dpT_RC_reco", title=r"p_{T}^{reco,RC} - p_{T}^{reco}", xlabel=r"#Deltap_{T}", n_bins=100, x_min=-2, x_max=2, units=r"GeV")
    h_dpT_RC_gen  = make_TH1F("h_dpT_RC_gen", title=r"p_{T}^{reco,RC} - p_{T}^{gen}", xlabel=r"#Deltap_{T}", n_bins=100, x_min=-5, x_max=5, units=r"GeV")
    h_dpT_reco_gen  = make_TH1F("h_dpT_reco_gen", title=r"p_{T}^{reco} - p_{T}^{gen}", xlabel=r"#Deltap_{T}", n_bins=100, x_min=-5, x_max=5, units=r"GeV")
    h_dpTOpT_RC_reco = make_TH1F("h_dpTOpT_RC_reco", title=r"(p_{T}^{reco,RC} - p_{T}^{reco})/p_{T}^{reco}", xlabel=r"#frac{#Deltap_{T}}{p_{T}}", n_bins=100, x_min=-0.02, x_max=0.02, units=None)
    h_dpTOpT_RC_gen = make_TH1F("h_dpTOpT_RC_gen", title=r"(p_{T}^{reco,RC} - p_{T}^{gen})/p_{T}^{gen}", xlabel=r"#frac{#Deltap_{T}}{p_{T}}", n_bins=100, x_min=-0.15, x_max=0.15, units=None)
    h_dpTOpT_reco_gen = make_TH1F("h_dpTOpT_reco_gen", title=r"(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", xlabel=r"#frac{#Deltap_{T}}{p_{T}}", n_bins=100, x_min=-0.15, x_max=0.15, units=None)

    hist_tup = (h_pT1, h_pT2, h_pT3, h_pT4,
        h_pT1_RC, h_pT2_RC, h_pT3_RC, h_pT4_RC,
        h_dpT_RC_reco, h_dpT_RC_gen, h_dpT_reco_gen,
        h_dpTOpT_RC_reco, h_dpTOpT_RC_gen, h_dpTOpT_reco_gen)
        
    for h in hist_tup:
        h.Sumw2()

    # RC = Rochester Corrections
    f = ROOT.TFile(inpath_file_noRC)
    f_RC = ROOT.TFile(inpath_file_RC)
    tree = f.Get("Ana/passedEvents")
    tree_RC = f_RC.Get("Ana/passedEvents")
    n_evts = tree.GetEntries()
    n_evts_RC = tree_RC.GetEntries()
    assert n_evts == n_evts_RC, "Files don't contain equal events!"

    run_over_n = n_evts if max_evts == -1 else max_evts
    print(f"Running over {run_over_n} events:")
    muon_coll = MyMuonCollection(prod_mode_ls=["ggH"])
    muon_coll.create_KinBins(eta_ls, pT_ls)
    for ctr, (evt, evt_RC) in enumerate(zip(tree, tree_RC)):
        if (ctr % 50000) == 0:
            print(f"...Processing event {ctr}")
        if ctr >= run_over_n:
            break
        # Make sure event number is the same.
        assert evt.Event == evt_RC.Event
        if (evt.mass4mu < 100) or (140 < evt.mass4mu):
            continue

        mu_ls = build_muons_from_HZZ4mu_event(tree, ctr, eta_bin=[eta_min, eta_max], pT_bin=[pT_min, pT_max], d0_max=1, use_FSR=False)
        mu_RC_ls = build_muons_from_HZZ4mu_event(tree_RC, ctr, eta_bin=[eta_min, eta_max], pT_bin=[pT_min, pT_max], d0_max=1, use_FSR=False)
        if (None in mu_ls) or (None in mu_RC_ls):
            continue
        # Associate the RC pT to each muon:
        for mu, muRC in zip(mu_ls, mu_RC_ls):
            mu.pT_RC = muRC.pT

        # Sort from highest to lowest pT:
        mu_ls = sorted(mu_ls, key=lambda mu : mu.pT)
        mu_ls.reverse()
        del mu_RC_ls  # For good measure.
        h_pT1.Fill(mu_ls[0].pT)
        h_pT2.Fill(mu_ls[1].pT)
        h_pT3.Fill(mu_ls[2].pT)
        h_pT4.Fill(mu_ls[3].pT)

        h_pT1_RC.Fill(mu_ls[0].pT_RC)
        h_pT2_RC.Fill(mu_ls[1].pT_RC)
        h_pT3_RC.Fill(mu_ls[2].pT_RC)
        h_pT4_RC.Fill(mu_ls[3].pT_RC)

        for mu in mu_ls:
            # Put into the correct kb2d.
            # if kb2d.check_muon_belongs(mu):
            #     kb2d.add_muon(mu)
            muon_coll.place_single_muon_into_KinBin(mu, eta_ls, pT_ls, verbose=verbose)
            h_dpT_RC_reco.Fill(mu.pT_RC - mu.pT)
            h_dpT_RC_gen.Fill(mu.pT_RC - mu.gen_pT)
            h_dpT_reco_gen.Fill(mu.pT - mu.gen_pT)

            h_dpTOpT_RC_reco.Fill((mu.pT_RC - mu.pT) / mu.pT)
            h_dpTOpT_RC_gen.Fill((mu.pT_RC - mu.gen_pT) / mu.gen_pT)
            h_dpTOpT_reco_gen.Fill((mu.pT - mu.gen_pT) / mu.gen_pT)
        muon_coll.muon_ls.extend(mu_ls)
    gc.collect()  # Can improve performance after a long loop.
    muon_coll.make_inclusive_kinematic_plots()  # Stores self.hist_inclusive_ls.
    print("...Finished looping over events.")
    print("...Putting muons into KinBin2Ds and deleting muons in MyMuonCollection...")
    muon_coll.muon_ls = "overwritten"
    print("...Plotting inclusive distributions.")
    c = ROOT.TCanvas()
    c.Print(outpath_file + "[")
    draw_Xhists_to_canvas(c, outpath_file, (h_pT1, h_pT2, h_pT3, h_pT4))
    draw_Xhists_to_canvas(c, outpath_file, (h_pT1_RC, h_pT2_RC, h_pT3_RC, h_pT4_RC))
    draw_Xhists_to_canvas(c, outpath_file, (h_dpT_RC_reco, h_dpT_RC_gen, h_dpT_reco_gen))
    draw_Xhists_to_canvas(c, outpath_file, (h_dpTOpT_RC_reco, h_dpTOpT_RC_gen, h_dpTOpT_reco_gen))
    for h in muon_coll.hist_inclusive_ls:
        h.Draw("hist")
        c.Print(outpath_file)
    c.Print(outpath_file + "]")
    print(f"...Saving distributions to root file:\n{outpath_root}")
    muon_coll.plot_ls.extend(hist_tup)

    outf = ROOT.TFile(outpath_root, "recreate")
    for h in muon_coll.plot_ls + muon_coll.hist_inclusive_ls:
        h.Write()
    outf.Close()

    save_to_pkl(muon_coll, outpath_pkl)
    print("Done.")