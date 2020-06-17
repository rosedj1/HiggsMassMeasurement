"""
PURPOSE:
  This code opens up one of Filippo's ggF files on HPG T2,
  selects m4mu events which satisfy the desired selections, 
  and makes a m4mu histogram. 
  This hist is saved to a '.root' file along with the 
  individual data in a TTree.

NOTES: This code will rewrite over outfile_path.
"""

import ROOT
import numpy as np
from array import array
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid#, fixOverlay
from Utils_Python.Utils_Files import check_overwrite
#----- User Parameters -----#
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2017.root"
outfile_path = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/MC_2017_ggF_m4mu.root"

n_evts = -1

overwrite = False
#----- Automatons -----#
check_overwrite(outfile_path, overwrite)

tdrStyle = setTDRStyle()
tdrGrid(tdrStyle, gridOn=True)

f = ROOT.TFile(infile_path, "read")
t = f.Get("Ana/passedEvents")

# New file, TTree, and TH1F.
outf = ROOT.TFile(outfile_path, "recreate")
newtree = ROOT.TTree("t1", "m4mu_uncorrected_pT")

h_m4mu = ROOT.TH1F("h_m4mu", "m_{4#mu} distribution, Higgs Production Mode: ggF", 140, 105, 140)
h_m4mu.Sumw2()

ptr = array('f', [0.])
newtree.Branch("m4mu", ptr, "m4mu/F")

#----- Functions -----#
def calc_Hmass(mu1, mu2, mu3, mu4):
    H = mu1 + mu2 + mu3 + mu4
    return H.M()

def passed_Higgs_selections(evt):
    """
    Check to see if this is a 4mu event. 
    
    Returns True, if this event passes the given selections.
    """
    selec_ls = []
    
    selec_ls.append(evt.passedFullSelection)
    selec_ls.append(evt.finalState == 1)
    selec_ls.append((105 < evt.mass4mu) & (evt.mass4mu < 140))
#     selec_ls.append(evt.passedFiducialSelection)
    
    # Also make sure to do: 
    # abs(eta) < 2.4
    # all muon pT > 5
    return all(selec_ls)

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
#     print("skimmed_ID_gen_ls:", skimmed_ID_gen_ls)
#     print("skimmed_ID_rec_ls:", skimmed_ID_rec_ls)
    check_IDs(skimmed_ID_gen_ls, skimmed_ID_rec_ls)
    
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
        
    # Make sure there are at least 2 OSSF pairs of muons.
    # t.GENlep_id examples:
    # [-13,13,-15,15] or [13,-13,-13,13,13,-13] or [13,-13,-13,13,11,-11]
    GENlep_id = list(t.GENlep_id)
    gteq2mupos = GENlep_id.count(13) >= 2
    gteq2muneg = GENlep_id.count(-13) >= 2
    if not (gteq2mupos and gteq2muneg):
        continue
        
    # At this point, sometimes particles will sneak into lep_genindex
    # but not really be muons, even though they were reconstructed as muons. 
    # MC 2017 ggF sample, event 21266 showed a ubar faking a muon. 
    # Sometimes we think we might have 4 or 5 muons, but 2 objects 
    # are actually fakes. Upon removing -1s from lep_genindex, we only
    # have 3 good muons. Definitely skip these events.
    lep_genindex = list(t.lep_genindex)
    n_possible_mu = len([x for x in lep_genindex if x != -1])
    if n_possible_mu < 4:
        continue
        
    skimmed_ID_rec_ls = get_skimmed_rec_ID_ls(t)
    
    # FIXME: Code could be improved here! 
#     if len(skimmed_ID_rec_ls) > 4:
#         ndcs_to_kill = get_ndcs_of_fake_muons(t)
        # FIXME: MUST PRUNE ALL lep_kinem vectors using ndcs_to_kill, like below:
#         lep_pt_good = remove_fake_muons_from_recvector(t.lep_pt, ndcs_to_kill)
    if len(skimmed_ID_rec_ls) != 4:
        continue
    
#     sanity_check(t)
    # Event is good. Do analysis.
    good_evts_adhoc += 1
    
    # Get indices
    rec_ncds_ls = get_4mu_ndcs_rec(t)
    gen_ncds_ls = get_4mu_ndcs_gen(t)
    


    index_rec_mu1 = rec_ncds_ls[0]
    index_rec_mu2 = rec_ncds_ls[1]
    index_rec_mu3 = rec_ncds_ls[2]
    index_rec_mu4 = rec_ncds_ls[3]
    
    mu1 = ROOT.Math.PtEtaPhiMVector(t.lep_pt[index_rec_mu1], t.lep_eta[index_rec_mu1], t.lep_phi[index_rec_mu1], t.lep_mass[index_rec_mu1])
    mu2 = ROOT.Math.PtEtaPhiMVector(t.lep_pt[index_rec_mu2], t.lep_eta[index_rec_mu2], t.lep_phi[index_rec_mu2], t.lep_mass[index_rec_mu2])
    mu3 = ROOT.Math.PtEtaPhiMVector(t.lep_pt[index_rec_mu3], t.lep_eta[index_rec_mu3], t.lep_phi[index_rec_mu3], t.lep_mass[index_rec_mu3])
    mu4 = ROOT.Math.PtEtaPhiMVector(t.lep_pt[index_rec_mu4], t.lep_eta[index_rec_mu4], t.lep_phi[index_rec_mu4], t.lep_mass[index_rec_mu4])
    
    m4mu = calc_Hmass(mu1, mu2, mu3, mu4)
    ptr[0] = m4mu
    newtree.Fill()
    h_m4mu.Fill(m4mu)
# End evt loop.

outf.cd()
h_m4mu.Write()
newtree.Write()
outf.Close()
f.Close()