"""
FIXME: This code looks like it was superceded by MyMuonCollection methods.
    - Investigate!
FIXME: This file needs to have its Higgs mass cuts figured out. 
    Particularly for lep_genindex! 
    [ ] Do gen indexing
    [ ] Do reco indexing
    [ ] Make sure 4 muons in event
    [ ] No electrons!
Purpose:
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
Updated: 2021-02-01
"""
import pickle
import ROOT
import numpy as np
from array import array

# Local imports.
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid#, fixOverlay
from Utils_Python.Utils_Files import check_overwrite
from d0_Studies.d0_Utils.d0_fns import correct_muon_pT
#----- User Parameters -----#
inpath_file = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2018.root"
inpath_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/MC2018_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0_fixedkeys.pkl"
outpath_file = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/m4mu_withcorr_MC2018_ggF_fullstats_pTlt200_mayhaveelec_d0cut.root"

eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
pT_binedge_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]

n_evts = -1

overwrite = True
verbose = True
#----- Functions -----#
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
    selec_ls.append((105 < evt.mass4mu) & (evt.mass4mu < 140))
#     selec_ls.append(evt.passedFiducialSelection)
    return all(selec_ls)

def check_2_OSSF_muon_pairs(id_ls):
    """
    Returns True if there are at least 2 OSSF pairs of muons.
    """
    gteq2mupos = id_ls.count(13) >= 2
    gteq2muneg = id_ls.count(-13) >= 2
    return all([gteq2mupos, gteq2muneg])

def check_atleast_4matchedmuons(evt):
    """Return True if there are at least 4 matched muons in this event."""
    # At this point, sometimes particles will sneak into lep_genindex
    # but not really be muons, even though they were reconstructed as muons. 
    # MC 2017 ggF sample, event 21266 showed a ubar faking a muon. 
    # Sometimes we think we might have 4 or 5 muons, but 2 objects 
    # are actually fakes. Upon removing -1's from lep_genindex, we only
    # have 3 good muons. Definitely skip these events.
    lep_genindex = list(evt.lep_genindex)
    n_possible_mu = len([x for x in lep_genindex if x != -1])
    return True if n_possible_mu >= 4 else False

# def taus_faking_muons(evt):
#     return any(abs(ID) == 15 for ID in evt.GENlep_id)

def get_skimmed_rec_ID_ls(evt):
    ID_rec_ls = list(evt.lep_id)
    indices_rec = get_4mu_ndcs_rec(evt)
    return [ID_rec_ls[ndx] for ndx in indices_rec]

def get_4mu_ndcs_rec(evt):
    """Returns a list of 
    
    The index of lep_genindex tells which rec lepton matches which gen.
    
    Example:
        GENlep_pt = [27.4, 9.0, 28.7, 18.1]  # gen pTs
        lep_pt = [27.4, 27.3, 18.2, 8.8]  # rec pTs
        
        Use lep_genindex = [0,2,3,1] to find which gen matches which rec. 
            0th element is 0 : lep_pt[0] matches GENlep_pt[0]
            1st element is 2 : lep_pt[1] matches GENlep_pt[2]
            2nd element is 3 : lep_pt[2] matches GENlep_pt[3]
            3rd element is 1 : lep_pt[3] matches GENlep_pt[1]"""
    indices_gen = get_4mu_ndcs_gen(evt)
    lep_genindex_ls = list(evt.lep_genindex)
    return [lep_genindex_ls.index(ndx) for ndx in indices_gen]

def get_4mu_ndcs_gen(evt):
    """ 
    Use these indices to slice GEN vectors. 
    Effectively gets rid of all -1. Maintains the original order.
    """
    return [ndx for ndx in evt.lep_genindex if ndx != -1]

def sanity_check(evt):
    """
    Takes a 4mu event, which has passed_Higgs_selections:
      * 105 < m4mu < 140 GeV
      * passedFullSelection == True
      * finalState == True
      
    Will raise an error if the kinematics of reco or gen look strange.
    """
    # Make sure most kinematic vectors all have the same length. 
    # All rec vectors are of same length; all gen vectors are of same length.
    equal_len_kinem_vec(evt)
    
    # Check that gen 4mu and rec 4mu match in IDs. 
    skimmed_ID_gen_ls = get_skimmed_gen_ID_ls(evt)
    skimmed_ID_rec_ls = get_skimmed_rec_ID_ls(evt)
    check_IDs(skimmed_ID_gen_ls, skimmed_ID_rec_ls)
    
def initialize_muon(evt, k):
    mu = ROOT.Math.PtEtaPhiMVector(evt.lep_pt[k], evt.lep_eta[k], evt.lep_phi[k], evt.lep_mass[k])
    mu.charge = evt.lep_id[k] / -13.
    mu.d0 = evt.lep_d0BS[k]
    return mu

def clone_muon(mu):
    cloned_mu = ROOT.Math.PtEtaPhiMVector(mu.Pt(), mu.Eta(), mu.Phi(), mu.M())
    cloned_mu.charge = mu.charge
    cloned_mu.d0 = mu.d0
    return cloned_mu

def make_muon_ls(evt):
    mu_ls = []
    # Run over the 4 rec muons.
    for k in range(4):
        mu = initialize_muon(evt, k)
        mu_ls.append(mu)
    return mu_ls

def make_muon_corr_ls(mu_ls, pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls, verbose=False):
    mu_corr_ls = []
    # Run over the 4 rec muons.
    for mu in mu_ls:
        # Find corrected pT for original muon.
        corr_pT = correct_muon_pT(mu.Eta(), mu.Pt(), mu.charge, mu.d0,
                              pT_corr_factor_dict=pT_corr_factor_dict,
                              eta_binedge_ls=eta_binedge_ls, pT_binedge_ls=pT_binedge_ls, verbose=verbose)
        # Make corr_reco mu.
        mu_corr = clone_muon(mu)
        # Assign it the corrected pT. 
        mu_corr.SetPt(corr_pT)
        # Save muons for this event.
        mu_corr_ls.append(mu_corr)

        assert mu.Pt() != mu_corr.Pt()
        
    return mu_corr_ls

def check_muon_kinem(mu):
    """Make sure that this muon has:
        0.0 < abs(eta) < 2.4
        5 < pT < 200
    """
    passed_eta = (0 < abs(mu.Eta())) and (abs(mu.Eta()) < 2.4)
    passed_pT = (5 < mu.Pt()) and (mu.Pt() < 200)
    passed_d0 = (abs(mu.d0) < 0.015)
    return all([passed_eta, passed_pT, passed_d0])

def verify_all_muons(mu_ls):
    """
    Make sure all muons passed selections and do final pT checks.
    Returns True if all muons passed, False otherwise.
    """
    count_pT_gt10_for2mu = 0
    count_pT_gt20_for1mu = 0
    for mu in mu_ls:
        pass_kinem = check_muon_kinem(mu)
        if not pass_kinem:
            return False
        if mu.Pt() > 10:
            count_pT_gt10_for2mu += 1
        if mu.Pt() > 20:
            count_pT_gt20_for1mu += 1
    # End loop over muons.
    if count_pT_gt10_for2mu < 2:
        return False
    elif count_pT_gt20_for1mu < 1:
        return False
    else:
        return True

#--------------------#
#----- Analysis -----#
#--------------------#
if __name__ == "__main__":
    check_overwrite(outpath_file, overwrite)

    tdrStyle = setTDRStyle()
    tdrGrid(tdrStyle, gridOn=True)

    print(f"...Opening root file:\n{inpath_file}")
    f = ROOT.TFile(inpath_file, "read")
    t = f.Get("Ana/passedEvents")

    with open(inpath_pkl, "rb") as pkl:
        print(f"...Opening pickle:\n{inpath_pkl}")
        pT_corr_factor_dict = pickle.load(pkl)

    # New file, TTree, and TH1F.
    outf = ROOT.TFile(outpath_file, "recreate")
    newtree = ROOT.TTree("tree", "tree_m4mu_vals")

    h_m4mu = ROOT.TH1F("h_m4mu", "m_{4#mu} distribution, Higgs Production Mode: ggF", 140, 105, 140)
    h_m4mu_corr = ROOT.TH1F("h_m4mu_corr", "m_{4#mu}^{p_{T},corr.} distribution (Higgs Production Mode: ggF)", 140, 105, 140)
    h_m4mu.Sumw2()
    h_m4mu_corr.Sumw2()

    ptr_m4mu = array('f', [0.])
    ptr_m4mu_corr = array('f', [0.])
    newtree.Branch("m4mu", ptr_m4mu, "m4mu/F")
    newtree.Branch("m4mu_corr", ptr_m4mu_corr, "m4mu_corr/F")

    all_evts = t.GetEntries()
    if n_evts == -1:
        n_evts = all_evts

    good_evts_adhoc = 0
    # Event loop.
    for ct,x in enumerate(range(n_evts)):
        if ct % 50000 == 0: 
            print("Running over evt:",ct)
        t.GetEntry(x)
        
        if not passed_Higgs_selections(t):
            continue
        if not check_2_OSSF_muon_pairs(list(t.GENlep_id)):
            continue
        if not check_2_OSSF_muon_pairs(list(t.lep_id)):
            continue
        if not validate_lep_genindex(t):
            continue

        def check_lep_genindex(evt):
            """
            Make sure lep_genindex has 4 unique leptons from Higgs.
            
            Explanation: lep_genindex can look like: [-1, 0, -1, 3, 1, 2]
                The elements 0,3,1,2 say that the 
            Corresponding lep_id may look like: [11, -13, -15, 13, 13, -13]
            """
            lep_genindex = list(evt.lep_genindex)
            occurrences = [lep_genindex.count(ndx) for ndx in range(4)]
            if sum()

            # for elem in lep_genindex:
            #     if elem in [0,1,2,3]:
            #         ndcs_4mu_fromHiggs_rec.append(elem.index())

            def get_ndcs_4lepsfromHiggs_rec(evt):
                """
                Return the list of indices of reco leptons which came from Higgs. 

                Suppose lep_genindex looks like: [-1, 0, -1, 3, 1, 2]
                Then this function will return the list of indices of the 
                4 leptons which came from the Higgs: [1, 3, 4, 5]
                Use this list to slice reco vectors, like: lep_id, lep_pt, lep_eta
                """
                lep_genindex = list(evt.lep_genindex)
                return [lep_genindex.index(elem) for elem in lep_genindex if elem in [0,1,2,3]]

            def get_ndcs_4lepsfromHiggs_gen(evt):
                """
                Return the list of indices of gen leptons which came from Higgs. 

                Suppose lep_genindex looks like: [-1, 0, -1, 3, 1, 2]
                Then this function will return the list of indices of the 
                4 leptons which came from the Higgs: [0, 3, 1, 2]
                Use this list to slice gen vectors, like: GENlep_id, GENlep_pt, GENlep_eta
                """
                lep_genindex = list(evt.lep_genindex)
                return [elem for elem in lep_genindex if elem in [0,1,2,3]]


            zeroes = lep_genindex.count(0)
            ones = lep_genindex.count(1)
            twos = 
            three = 


        if not check_2_OSSF_pairs(t):
            continue
        if not check_atleast_4matchedmuons(t):
            continue
        
        # Now, lep_genindex definitely has at least 4 muons.
        # if not check_4recomuons(t):
        #     continue
        print(f"t.lep_id looks like: {list(t.lep_id)}")
        skimmed_ID_rec_ls = get_skimmed_rec_ID_ls(t)
        
        # FIXME: Code could be improved here! 
    #     if len(skimmed_ID_rec_ls) > 4:
    #         ndcs_to_kill = get_ndcs_of_fake_muons(t)
            # FIXME: MUST PRUNE ALL lep_kinem vectors using ndcs_to_kill, like below:
    #         lep_pt_good = remove_fake_muons_from_recvector(t.lep_pt, ndcs_to_kill)
        if len(skimmed_ID_rec_ls) != 4:
            continue
        
    #     sanity_check(t)
        # Event looks good so far. 
        # Now check kinematics of muons.
        mu_ls = make_muon_ls(t)
        all_muons_passed = verify_all_muons(mu_ls)
        if not all_muons_passed:
            continue

        # NOW event is good.
        good_evts_adhoc += 1
        # Correct the pT of each muon. 
        mu_corr_ls = make_muon_corr_ls(mu_ls, pT_corr_factor_dict, eta_binedge_ls, pT_binedge_ls, verbose=verbose)

        assert all([mu.Pt() != mucorr.Pt() for mu,mucorr in zip(mu_ls,mu_corr_ls)])

        m4mu = calc_Hmass(*mu_ls)
        m4mu_corr = calc_Hmass(*mu_corr_ls)
        
        ptr_m4mu[0] = m4mu
        ptr_m4mu_corr[0] = m4mu_corr
        newtree.Fill()
        h_m4mu.Fill(m4mu)
        h_m4mu_corr.Fill(m4mu_corr)
    # End evt loop.

    print(f"Found {good_evts_adhoc} good m4mu events after selections.")
    outf.cd()
    print(f"Saving m4mu histograms and TTree to file:\n{outpath_file}")
    h_m4mu.Write()
    h_m4mu_corr.Write()
    newtree.Write()
    outf.Close()
    f.Close()