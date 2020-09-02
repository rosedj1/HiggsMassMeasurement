"""
PURPOSE:
    This code opens up a gluon fusion Higgs sample
    from Filippo's area on HPG T2.
    It selects m4mu events which pass all specified selection criteria,
    stores the m4mu and m4mu_corr values in a TTree, and makes 
    distributions of these 2 variables.
    m4mu_corr is the m4mu reevaluated after each muon has
    had its pT corrected from the d0 studies. 
Syntax: 
Notes: 
Author: Jake Rosenzweig
Created:
Updated: 
"""
import pickle
import ROOT
import numpy as np
from array import array

# Local imports.
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid#, fixOverlay
from Utils_Python.Utils_Files import check_overwrite
from d0_Studies.d0_Utils.d0_fns import correct_muon_pT
from array import array
#----- User Parameters -----#
year = 2017
inpath_file = f"/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_{year}.root"
inpath_pkl = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/MC{year}_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0.pkl"
outpath_file_woext = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/MC{year}_ggF_synchwithFilippo_basiccuts_usingFSR_absd0cut0p010"

eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
pT_binedge_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]

n_evts = -1
overwrite = False
verbose = False
#----- Functions -----#
def make_outfiles(outpath_file_woext, overwrite):
    """Make outgoing '.root' and '.pdf' files. Also throw error if they already exist."""
    outpath_rootfile = outpath_file_woext + ".root"
    outpath_pdf = outpath_file_woext + ".pdf"

    check_overwrite(outpath_rootfile, overwrite)
    check_overwrite(outpath_pdf, overwrite)
    return outpath_rootfile, outpath_pdf

def calc_Hmass(mu1, mu2, mu3, mu4):
    """
    NOTE: This method relies on calls like: mu.Pt(), mu.Eta(), etc. 
    So you must change the muon kinematics by doing: mu.SetPt(val), e.g.
    """
    H = mu1 + mu2 + mu3 + mu4
    return H.M()

def passed_Higgs_selections(evt):
    """
    Returns True, if this event passes the given selections.
    """
    selec_ls = []
    selec_ls.append(evt.passedFullSelection)
    selec_ls.append(evt.finalState == 1)
    selec_ls.append((105 < evt.mass4l) & (evt.mass4l < 140))
#     selec_ls.append(evt.passedFiducialSelection)
    return all(selec_ls)

def check_2_OSSF_muon_pairs(id_ls):
    """
    Returns True if there are at least 2 OSSF pairs of muons.
    """
    gteq2mupos = id_ls.count(13) >= 2
    gteq2muneg = id_ls.count(-13) >= 2
    return all([gteq2mupos, gteq2muneg])

def get_ndcs_rec(lep_genindex_ls):
    """
    Returns a list of only the indices in lep_kinem list
    which have a value >= 0.

    WARNING: 
    - Can return a list with length >= 4.
      Example: 
        lep_genindex_ls = [3, 2, -1, 0, 1, -1, 3]
        returns: [0, 1, 3, 4, 6]
      This is problematic because the 0th reco and the 6th reco leps
      were matched to the same gen. Skip this event!
    
    The index of lep_genindex tells which rec lepton matches which gen.
    
    Example:
        GENlep_pt = [27.4, 9.0, 28.7, 18.1]  # gen pTs
        lep_pt = [27.4, 27.3, 18.2, 8.8]  # rec pTs
        
        Use lep_genindex = [0,2,3,1] to find which gen matches which rec. 
            0th element is 0 : lep_pt[0] matches GENlep_pt[0]
            1st element is 2 : lep_pt[1] matches GENlep_pt[2]
            2nd element is 3 : lep_pt[2] matches GENlep_pt[3]
            3rd element is 1 : lep_pt[3] matches GENlep_pt[1]"""
    return [ct for ct,x in enumerate(lep_genindex_ls) if x >= 0]

def get_ndcs_gen(rec_ndcs_ls, lep_genindex):
    """ 
    Use the good reco indices to retrieve the corresponding gen indices.
    These gen indices can be used to slice GEN vectors. 

    rec_ndcs_ls = [2, 1, 0, 3]
    lep_genindex = [2, 3, 1, 0, -1]
    returns: [1, 3, 2, 0]  # good gen indices to slices GEN vectors.
    """
    return [lep_genindex[ndx] for ndx in rec_ndcs_ls]

def validate_lep_genindex(lep_genindex_ls):
    """
    Return True, if lep_genindex_ls satisfies the following:
        - Make sure it has at least 4 positive elements:
            - [2, 3, 0, 1]
            - [2, 0, 1, -1, -1, 3]
        - For indices >= 0, make sure there are no duplicate indices. 
    """
    # Find all positive elements.
    pruned_ls = [x for x in lep_genindex_ls if x >= 0]
    if len(pruned_ls) < 4:
        return False
    elif len(pruned_ls) != len(set(pruned_ls)):
        # Duplicates detected.
        return False 
    else:
        return True

def initialize_muon(evt, reco_ndx, gen_ndx):
    """
    For a given event (evt), build a muon (mu) using: 
        - the reco values at a given index (reco_ndx) in the lep_kinem vectors
        - the gen values at a given gen index (gen_ndx) in the GENlep_kinem vectors. 

    NOTE:
        - The reco values are retrieved using ROOT.Math.LorentzVector methods (mu.Pt(), mu.Phi(), etc.)
        - If you use this muon in ROOT calculations like: mu1 + mu2,
          then be sure that the kinematics are stored using ROOT.Math.SetPtEtaPhiM(pt,eta,phi,m).

    Additionally: 
    - checks ID to make sure it is a muon.
    - Stores: pT, eta, phi, mass, charge, d0_BS.

    Returns
    -------
    mu : ROOT.Math.LorentzVector
    """
    # Record reco info. 
    mu = ROOT.Math.PtEtaPhiMVector(evt.lepFSR_pt[reco_ndx], evt.lepFSR_eta[reco_ndx], evt.lepFSR_phi[reco_ndx], evt.lepFSR_mass[reco_ndx])
    mu.charge = evt.lep_id[reco_ndx] / -13.
    if abs(mu.charge) != 1:
        print(f"mu.charge ({mu.charge}]) != +-1")
        raise ValueError
    mu.d0 = evt.lep_d0BS[reco_ndx]
    # Gen kinematics.
    mu.gen_Pt = evt.GENlep_pt[gen_ndx]
    mu.gen_Eta = evt.GENlep_eta[gen_ndx]
    mu.gen_Phi = evt.GENlep_phi[gen_ndx]
    mu.gen_Mass = evt.GENlep_mass[gen_ndx]
    # Kinematical combinations.
    mu.dpTOverpT = (mu.Pt() - mu.gen_Pt) / mu.gen_Pt
    return mu

def clone_muon(mu):
    cloned_mu = ROOT.Math.PtEtaPhiMVector(mu.Pt(), mu.Eta(), mu.Phi(), mu.M())
    cloned_mu.charge = mu.charge
    cloned_mu.d0 = mu.d0
    return cloned_mu

def make_muon_ls(evt, rec_ndcs_ls, gen_ndcs_ls):
    """Create a list of muons with reco and gen kinematic info.
    
    Parameters
    ----------
    evt : TTree event (e.g. tree.GetEntry(2))
    rec_ndcs_ls : list
        The reco indices of the muons in a "lep vector" (e.g. lep_pt, lep_id).
    gen_ndcs_ls : list
        The gen indices of a lep_gen vector that correspond to the matched reco muons. 
        Should be the same length as rec_ndcs_ls.
    """
    assert len(rec_ndcs_ls) == len(gen_ndcs_ls)
    mu_ls = []
    # Run over the reco muons.
    for reco_ndx,gen_ndx in zip(rec_ndcs_ls, gen_ndcs_ls):
        mu = initialize_muon(evt, reco_ndx, gen_ndx)
        mu_ls.append(mu)
    return mu_ls

def make_muon_corr_ls(mu_ls, pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls, verbose=False):
    mu_corr_ls = []
    # Run over the 4 rec muons.
    for mu in mu_ls:
        # Make corr_reco mu.
        # mu_corr = clone_muon(mu)
        # Apply pT correction only works if muon is found 
        # Find corrected pT for original muon.
        corr_pT = correct_muon_pT(mu.Eta(), mu.Pt(), mu.charge, mu.d0,
                            pT_corr_factor_dict=pT_corr_factor_dict,
                            eta_binedge_ls=eta_binedge_ls, pT_binedge_ls=pT_binedge_ls, verbose=verbose)
        # Assign it the corrected pT. 
        mu_corr = ROOT.Math.PtEtaPhiMVector(corr_pT, mu.Eta(), mu.Phi(), mu.M())
        mu_corr.charge = mu.charge
        mu_corr.d0 = mu.d0
        # mu_corr.SetPt(corr_pT)
        # Save muons for this event.
        mu_corr_ls.append(mu_corr)
    return mu_corr_ls

def check_muon_kinem(mu):
    """Make sure that this muon has reco info:
        0.0 < abs(eta) < 2.4
        5 < pT < 200
    """
    passed_eta = True
    passed_pT =  True
    # passed_eta = (0 < abs(mu.Eta())) and (abs(mu.Eta()) < 2.4)
    # passed_pT = (5 < mu.Pt()) and (mu.Pt() < 200)
    passed_d0 = True #(abs(mu.d0) < 0.010)
    return all([passed_eta, passed_pT, passed_d0])

def verify_all_muons(mu_ls):
    """
    Make sure all muons passed selections and do final pT checks.
    Returns True if all muons passed, False otherwise.
    """
    # count_pT_gt10_for2mu = 0
    # count_pT_gt20_for1mu = 0
    for mu in mu_ls:
        pass_kinem = check_muon_kinem(mu)
        if not pass_kinem:
            return False
        # if mu.Pt() > 10:
        #     count_pT_gt10_for2mu += 1
        # if mu.Pt() > 20:
        #     count_pT_gt20_for1mu += 1
    # End loop over muons.
    # if count_pT_gt10_for2mu < 2:
    #     return False
    # elif count_pT_gt20_for1mu < 1:
    #     return False
    # else:
        # return True
    return True

def check_matched_IDs(rec_ndcs_ls, gen_ndcs_ls, rec_ID_ls, gen_ID_ls):
    """
    Use the reco indices and the gen indices to slice the corresponding
    ID list and see if the IDs match. Return True if all IDs match. 

    OLD DOC STRING:
    Use the indices and elements of lep_genindex to see if the IDs
    between matched reco and gen leptons match. Returns True if 
    all IDs match.
    
    The elements of lep_genindex are the indices of the gen particles.
    The indices of lep_genindex are the indices of the reco particles.

    Example:
        GENlep_pt = [27.4, 9.0, 28.7, 18.1]  # gen pTs
        lep_pt = [27.4, 27.3, 18.2, 8.8]  # rec pTs
        
        Use lep_genindex = [0,2,3,1] to find which gen matches which rec. 
            0th element is 0 : lep_pt[0] matches GENlep_pt[0]
            1st element is 2 : lep_pt[1] matches GENlep_pt[2]
            2nd element is 3 : lep_pt[2] matches GENlep_pt[3]
            3rd element is 1 : lep_pt[3] matches GENlep_pt[1]
    """
    check_ID_ls = []
    for rec_elem,gen_elem in zip(rec_ndcs_ls, gen_ndcs_ls):
        match_bool = (rec_ID_ls[rec_elem] == gen_ID_ls[gen_elem])
        check_ID_ls.append(match_bool)
    return all(check_ID_ls)

def print_muon_info(mu, mu_corr):
    print(f"  mu.charge={mu.charge:18.8f}    mu_corr.charge={mu_corr.charge:18.8f}")
    print(f"  mu.d0=    {mu.d0:18.8f}    mu_corr.d0=    {mu_corr.d0:18.8f}")
    print(f"  mu.Pt()=  {mu.Pt():18.8f}    mu_corr.Pt()=  {mu_corr.Pt():18.8f}")
    print(f"  mu.Eta()= {mu.Eta():18.8f}    mu_corr.Eta()= {mu_corr.Eta():18.8f}")
    print(f"  mu.Phi()= {mu.Phi():18.8f}    mu_corr.Phi()= {mu_corr.Phi():18.8f}")
    print(f"  mu.M()=   {mu.M():18.8f}    mu_corr.M()=   {mu_corr.M():18.8f}\n")

def print_muon_ls_info(mu_ls, mu_corr_ls):
    for mu, mu_corr in zip(mu_ls, mu_corr_ls):
        print_muon_info(mu, mu_corr)
#--------------------#
#----- Analysis -----#
#--------------------#
global_atleast2mu_pTgt10 = 0 # Delete this later
if __name__ == "__main__":
    ROOT.gROOT.SetBatch(ROOT.kTRUE)
    outpath_rootfile, outpath_pdf = make_outfiles(outpath_file_woext, overwrite)

    print(f"...Opening root file:\n{inpath_file}")
    f = ROOT.TFile(inpath_file, "read")
    t = f.Get("Ana/passedEvents")

    with open(inpath_pkl, "rb") as pkl:
        print(f"...Opening pickle:\n{inpath_pkl}")
        pT_corr_factor_dict = pickle.load(pkl)

    # New file, TTree, and TH1Fs.
    outf = ROOT.TFile(outpath_rootfile, "recreate")
    newtree = ROOT.TTree("tree", "tree_m4mu_vals")

    h_m4mu = ROOT.TH1F("h_m4mu", "m_{4#mu} distribution, Higgs Production Mode: ggF", 140, 105, 140)
    h_m4mu_corr = ROOT.TH1F("h_m4mu_corr", "m_{4#mu}^{p_{T},corr.} distribution (Higgs Production Mode: ggF)", 140, 105, 140)
    h_m4mu_diff = ROOT.TH1F("h_m4mu_diff", "#Deltam_{4#mu} #equiv m_{4#mu}^{corr. p_{T}} - m_{4#mu}", 200, -20, 20)
    h_m4muvsm4mucorr = ROOT.TH2F("h_m4muvsm4mucorr", "Correlation between m_{4#mu}^{corr. p_{T}} and m_{4#mu}", 
                                 100, 70, 170, 100, 70, 170)  # n_binX, X_Low,X_Hig, n_binY, Y_low, Y_high
    h_deltapT_corr = ROOT.TH1F("h_deltapT_corr", "p_{T} correction", 200, -5, 5)
    h_m4mu.Sumw2()
    h_m4mu_corr.Sumw2()
    h_m4mu_diff.Sumw2()
    h_m4muvsm4mucorr.Sumw2()
    h_deltapT_corr.Sumw2()

    ptr_m4mu = array('f', [0.])
    ptr_m4mu_corr = array('f', [0.])
    ptr_m4mu_diff = array('f', [0.])
    newtree.Branch("m4mu", ptr_m4mu, "m4mu/F")
    newtree.Branch("m4mu_corr", ptr_m4mu_corr, "m4mu_corr/F")

    all_evts = t.GetEntries()
    if n_evts == -1:
        n_evts = all_evts

    good_evts_adhoc = 0
    m4mu_ls = []
    m4mu_corr_ls = []
    # Event loop.
    for ct in range(n_evts):
        if ct % 50000 == 0:
            print("Running over evt:",ct)
        t.GetEntry(ct)
        
        if not passed_Higgs_selections(t):
            continue

        lep_Hindex_ls = list(t.lep_Hindex)  # Elements that correspond to 4 leptons which build Higgs candidate.
        rec_ndcs_ls = lep_Hindex_ls
        lep_genindex_ls = list(t.lep_genindex)
        # if not validate_lep_genindex(lep_genindex_ls):
        #     continue

        # We have a good lep_genindex: at least 4 leptons have been matched. 
        # Could be 5 or more leps.
        # rec_ndcs_ls = get_ndcs_rec(lep_genindex_ls)
        gen_ndcs_ls = get_ndcs_gen(rec_ndcs_ls, lep_genindex_ls)

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
            continue

        # NOW event is good.
        good_evts_adhoc += 1
        # Correct the pT of each muon. 
        mu_corr_ls = make_muon_corr_ls(mu_ls, pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls, verbose=verbose)

        if len(rec_ndcs_ls) != 4:
            print(f"[WARNING] Event {ct}: len(rec_ndcs_ls) != 4, (it is {len(rec_ndcs_ls)})")
            print_muon_ls_info(mu_ls, mu_corr_ls)
            # continue
        if (-1 in rec_ndcs_ls):
            print(f"[WARNING] Event {ct}: -1 in rec_ndcs_ls, (rec_ndcs_ls={rec_ndcs_ls})")
        m4mu = calc_Hmass(*mu_ls)
        m4mu_corr = calc_Hmass(*mu_corr_ls)
        m4mu_diff = m4mu_corr - m4mu
        
        m4mu_ls.append(float(m4mu))
        m4mu_corr_ls.append(float(m4mu_corr))

        # Seeing strange m4mu calculation. Investigation:
        def show_diff_array(kinem, fsr_ls, ls):
            fsr_arr = np.array(fsr_ls)
            arr = np.array(ls)
            print(f"{kinem}: fsr - wofsr = {fsr_arr - arr}")

        in_bounds = (105 < m4mu) and (m4mu < 140)
        in_bounds_corr = (105 < m4mu_corr) and (m4mu_corr < 140)
        if (not in_bounds) or (not in_bounds_corr):
            print(f"Event={ct}, m4mu={m4mu}, m4mu_corr={m4mu_corr}, nFSRPhotons={t.nFSRPhotons}")
            show_diff_array("pT", list(t.lepFSR_pt), list(t.lep_pt))
            show_diff_array("eta", list(t.lepFSR_eta), list(t.lep_eta))
            show_diff_array("phi", list(t.lepFSR_phi), list(t.lep_phi))
            show_diff_array("mass", list(t.lepFSR_mass), list(t.lep_mass))
            print_muon_ls_info(mu_ls, mu_corr_ls)

        rel_diff = m4mu_diff / m4mu
        if abs(rel_diff) > 0.05:
            print(f"event {ct}:")
            print(f"  m4mu={m4mu}, m4mu_corr={m4mu_corr}, rel_diff={rel_diff}")
            print_muon_ls_info(mu_ls, mu_corr_ls)

        ptr_m4mu[0] = m4mu
        ptr_m4mu_corr[0] = m4mu_corr
        newtree.Fill()

        h_m4mu.Fill(m4mu)
        h_m4mu_diff.Fill(m4mu_diff)
        h_m4mu_corr.Fill(m4mu_corr)
        h_m4muvsm4mucorr.Fill(m4mu, m4mu_corr)
        for mu_,mucorr_ in zip(mu_ls, mu_corr_ls):
            shift = mu_.Pt() - mucorr_.Pt()
            h_deltapT_corr.Fill(shift)
    # End evt loop.
    
    graph = ROOT.TGraph(good_evts_adhoc, array('f', m4mu_ls), array('f', m4mu_corr_ls))
    graph.GetXaxis().SetTitle("m_{4#mu} [GeV]")
    graph.GetYaxis().SetTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
    graph.SetTitle("Effect of p_{T} corrections from d_{0} studies on m_{4#mu}")
    graph.SetMarkerColor(ROOT.kBlue)

    def fix_stats_box_and_draw(hist, canv, dim=1):
        if dim == 1:
            hist.Draw("hist")
        elif dim == 2:
            hist.Draw("colz")
        ROOT.gPad.Update()
        statsbox = hist.FindObject("stats")
        statsbox.SetX1NDC(0.75)
        statsbox.SetX2NDC(0.90)
        statsbox.SetY1NDC(0.75)
        statsbox.SetY2NDC(0.90)
        if dim == 1:
            hist.Draw("hist")
        elif dim == 2:
            hist.Draw("colz")
        canv.Update()

    print(f"Drawing histograms to:\n{outpath_pdf}")
    c1 = ROOT.TCanvas()
    c1.SetTicks(1,1)
    ROOT.gStyle.SetOptStat("iouRMe")
    c1.Print(outpath_pdf + "[")
    h_m4mu.SetXTitle("m_{4#mu} [GeV]")
    h_m4mu.SetYTitle("Events")
    fix_stats_box_and_draw(h_m4mu, c1, dim=1)
    # h_m4mu.Draw("hist")
    c1.Print(outpath_pdf)
    h_m4mu_corr.SetXTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
    h_m4mu_corr.SetYTitle("Events")
    fix_stats_box_and_draw(h_m4mu_corr, c1, dim=1)
    c1.Print(outpath_pdf)
    h_m4mu_diff.SetXTitle("#Deltam_{4#mu} [GeV]")
    h_m4mu_diff.SetYTitle("Events")
    fix_stats_box_and_draw(h_m4mu_diff, c1, dim=1)
    c1.Print(outpath_pdf)
    h_m4muvsm4mucorr.SetXTitle("m_{4#mu} [GeV]")
    h_m4muvsm4mucorr.SetYTitle("m_{4#mu}^{corr. p_{T}} [GeV]")
    fix_stats_box_and_draw(h_m4muvsm4mucorr, c1, dim=2)
    c1.Print(outpath_pdf)
    h_deltapT_corr.SetXTitle("#Deltap_{T} [GeV]")
    h_deltapT_corr.SetYTitle("Events")
    fix_stats_box_and_draw(h_deltapT_corr, c1, dim=1)
    c1.Print(outpath_pdf)
    c1.Print(outpath_pdf + "]")

    outf.cd()
    print(f"Saving hists and TTree to file:\n{outpath_rootfile}")
    h_m4mu.Write()
    h_m4mu_corr.Write()
    h_m4mu_diff.Write()
    newtree.Write()

    outf.Close()
    f.Close()
    print(f"Found {good_evts_adhoc} good m4mu events after selections.")
