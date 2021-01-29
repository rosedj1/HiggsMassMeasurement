from Particles import MyMuon

class Selector():
    """
    A class to organize the selections for different samples. 
    """
    def __init__(self, sample_name, verbose=False):
        """
        sample_name : str
            "Jpsi", "DY", "Upsilon", "Higgs"
        verbose : bool, optional
            Show debug info.
        """
#         pass_overall=None
#         pass_overall : bool  
#             If False, failed at least 1 selection. 
#             If True, has passed all selections.
#             If None, has not had any selections applied.
#         self.pass_overall = pass_overall
            
        self.sample_name = sample_name
        self.verbose = verbose
        
#         self.passed_dR_cut = False
#         self.passed_m2l_cut = False
#         self.passed_pT_cut = False
#         self.passed_eta_cut = False
        
        self.cut_dict = {
            "Jpsi"    : {"m2l_cut" : [2.9, 3.3], "dR_cut" : 0.005},  # Filippo: 0.005
            "DY"      : {"m2l_cut" : [60, 120],  "dR_cut" : 0.002},  # Filippo: 0.002
            "Upsilon" : {"m2l_cut" : [8.5, 11],  "dR_cut" : 0.005},  # Filippo: 0.005
            }
        
#     def check_passed_all_selections(self):
#         self.passed_all = all( (self.passed_dR_cut, 
#                                 self.passed_m2l_cut,
#                                 self.passed_pT_cut,
#                                 self.passed_eta_cut)
#                              )
    
    def pass_selection_val(self, val, val_min=None, val_max=None):
        """
        Return True if val is within val_min, val_max.
        Otherwise, returns False.
        
        Useful for applying pT cuts but is used in other methods. 
        """
        if (val_min is not None) and (val_max is not None):
            return True if (val > val_min) and (val < val_max) else False
        elif (val_min is None) and (val_max is not None):
            # Check against val_max only. 
            return True if (val < val_max) else False
        elif (val_min is not None) and (val_max is None):
            # Check against val_min only. 
            return True if (val > val_min) else False
        else:
            msg = "[WARNING] You performed a cut, but didn't specify any bounds."
            raise RuntimeWarning(msg)
        
    def pass_selection_abseta(self, eta, eta_min=None, eta_max=None):
        """
        Return True if abs(eta) is within eta_min, eta_max bounds.
        Otherwise, returns False.
        """
        return self.pass_selection_val(abs(eta), eta_min, eta_max)

    def pass_selection_m2l(self, m2l):
        """Return True if m2l passes selections, based on sample_name."""
        m2l_cut_ls = self.cut_dict[self.sample_name]["m2l_cut"]
        return self.pass_selection_val(val=m2l, val_min=m2l_cut_ls[0], val_max=m2l_cut_ls[1])
    
    def pass_selection_dR(self, dR):
        """Return True if dR passes selections, based on sample_name."""
        dR_max = self.cut_dict[self.sample_name]["dR_cut"]
        return self.pass_selection_val(val=dR, val_max=dR_max)

#----------------------------#

def passed_Higgs_evt_selection(evt):
    """
    Returns True, if this event passes the given selections.

    evt : TTree
        NOTE: evt must have attributes:
        - passedFullSelection
        - finalState
        - mass4l
    """
    selec_ls = []
    selec_ls.append(evt.passedFullSelection)
    selec_ls.append(evt.finalState == 1)
    selec_ls.append((105 < evt.mass4l) & (evt.mass4l < 140))
#     selec_ls.append(evt.passedFiducialSelection)
    return all(selec_ls)

def passed_Hmumu_evt_selection(evt):
    """
    TODO: implement Hmumu selections here, if any. 

    Returns True, if this event passes the given selections.

    evt : TTree
        NOTE: evt must have attributes:
    """
    selec_ls = [True]
    return all(selec_ls)

def initialize_muon(evt, reco_ndx, gen_ndx, use_FSR=False):
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
    mu : MyMuon object
    """
    # Record reco info. 
    # mu = ROOT.Math.PtEtaPhiMVector(evt.lepFSR_pt[reco_ndx], evt.lepFSR_eta[reco_ndx], evt.lepFSR_phi[reco_ndx], evt.lepFSR_mass[reco_ndx])
    charge = evt.lep_id[reco_ndx] / -13.
    mu = MyMuon(charge)
    # Set reco and gen kinematics.
    mu.set_GENPtEtaPhiMass(evt.GENlep_pt[gen_ndx], evt.GENlep_eta[gen_ndx], evt.GENlep_phi[gen_ndx], evt.GENlep_mass[gen_ndx])
    if use_FSR:
        this_pT = evt.lepFSR_pt[reco_ndx]
    else:
        # On 2020-12-16 Filippo says we should use WITHOUT FSR quantities.
        this_pT = evt.lep_pt[reco_ndx]
    mu.set_PtEtaPhiMass(this_pT, evt.lep_eta[reco_ndx], evt.lep_phi[reco_ndx], evt.lep_mass[reco_ndx])
    # Set other kinematic quantities.
    mu.d0 = evt.lep_d0BS[reco_ndx]
    mu.dpTOverpT = (mu.pT - mu.gen_pT) / mu.gen_pT
    return mu

def initialize_Htomumu_muons_fromliteskim(evt, num):
    """
    For a given event (evt), build a muon (mu) using: 
        - the reco values in the lep_kinem vectors
        - the gen values in the GENlep_kinem vectors. 

    NOTE:
        - The reco values are retrieved using ROOT.Math.LorentzVector methods (mu.Pt(), mu.Phi(), etc.)
        - If you use this muon in ROOT calculations like: mu1 + mu2,
          then be sure that the kinematics are stored using ROOT.Math.SetPtEtaPhiM(pt,eta,phi,m).

    Additionally: 
    - checks ID to make sure it is a muon.
    - Stores: pT, eta, phi, mass, charge, d0_BS.

    Parameters
    ----------
    evt : TTree event
    
    num : int
        Which muon in this 2 mu event. Either `1` or `2`.
    Returns
    -------
    mu : MyMuon object
    """
    assert num in [1, 2]
    # Record reco info. 
    # mu = ROOT.Math.PtEtaPhiMVector(evt.lepFSR_pt[reco_ndx], evt.lepFSR_eta[reco_ndx], evt.lepFSR_phi[reco_ndx], evt.lepFSR_mass[reco_ndx])
    charge = getattr(evt, f"Id{num}") / -13.0
    mu = MyMuon(charge)
    # Set reco and gen kinematics.
    # mu.set_PtEtaPhiMass(getattr(evt, f"pT_FSR{num}"),
    #                     getattr(evt, f"eta_FSR{num}"),
    #                     getattr(evt, f"phi_FSR{num}"),
    #                     getattr(evt, f"m_FSR{num}")
    #                     )
    # On 2020-12-16 Filippo says we should use WITHOUT FSR quantities.
    mu.set_PtEtaPhiMass(getattr(evt, f"pT{num}"),
                        getattr(evt, f"eta{num}"),
                        getattr(evt, f"phi{num}"),
                        getattr(evt, f"m{num}")
                        )
    mu.set_GENPtEtaPhiMass(getattr(evt, f"genLep_pt{num}"),
                           getattr(evt, f"genLep_eta{num}"),
                           getattr(evt, f"genLep_phi{num}"),
                           getattr(evt, f"genLep_mass{num}")
                           )
    # Set other kinematic quantities.
    mu.d0 = getattr(evt, f"d0BS{num}")
    mu.dpTOverpT = (mu.pT - mu.gen_pT) / mu.gen_pT
    return mu

def make_muon_ls(evt, rec_ndcs_ls, gen_ndcs_ls, use_FSR=False):
    """Return a list of MyMuon objects with stored reco and gen kinematic info.
    
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
    # mu_ls = []
    # Run over the reco muons.
    # for reco_ndx,gen_ndx in zip(rec_ndcs_ls, gen_ndcs_ls):
    #     mu = initialize_muon(evt, reco_ndx, gen_ndx)
    #     mu_ls.append(mu)
    # return mu_ls
    return [initialize_muon(evt, reco_ndx, gen_ndx, use_FSR) for reco_ndx,gen_ndx in zip(rec_ndcs_ls, gen_ndcs_ls)]

def make_muon_ls_fromliteskim(evt):
    """Return a 2-elem list of MyMuon objects with stored reco and gen info.

    This function assumes that reco-gen index matching (via lep_genindex)
    was already done in a lite skim (like MioSkim_2L.C).
    
    Parameters
    ----------
    evt : TTree event (e.g. tree.GetEntry(2))
    """
    return [initialize_Htomumu_muons_fromliteskim(evt, num) for num in [1,2]]

# def check_muon_kinem(mu):
#     """Make sure that this muon has reco info:
#         0.0 < abs(eta) < 2.4
#         5 < pT < 200
#     """
    
#     return 

def apply_kinem_selections(mu_ls, eta_bin=[0,2.4], pT_bin=[5,200], d0_max=1):
    """
    Make sure all muons passed selections and do final pT checks.
    Returns True if all muons passed, False otherwise.

    Parameters
    ----------
    eta_bin : 2-elem list
        Keep muons with reco abs(eta) in range: [eta_min, eta_max]
    pT_bin : 2-elem list
        Keep muons with reco pT in range: [pT_min, pT_max] (GeV)
    d0_max : float
        Keep muons with abs(d0) < d0_max (cm).
    """
    eta_min, eta_max = eta_bin[0], eta_bin[1]
    pT_min, pT_max = pT_bin[0], pT_bin[1]
    for mu in mu_ls:
        passed_eta = (eta_min < abs(mu.eta)) and (abs(mu.eta) < eta_max)
        passed_pT = (pT_min < mu.pT) and (mu.pT < pT_max)
        passed_d0 = abs(mu.d0) < d0_max
        pass_kinem = all([passed_eta, passed_pT, passed_d0])
        if not pass_kinem:
            return False
    return True

def get_ndcs_gen(rec_ndcs_ls, lep_genindex):
    """ 
    Use the good reco indices to retrieve the corresponding gen indices.
    These gen indices can be used to slice GEN vectors. 

    rec_ndcs_ls = [2, 1, 0, 3]
    lep_genindex = [2, 3, 1, 0, -1]
    returns: [1, 3, 2, 0]  # good gen indices to slices GEN vectors.
    """
    return [lep_genindex[ndx] for ndx in rec_ndcs_ls]

def build_muons_from_HZZ4mu_event(t, evt_num, eta_bin=[0,2.4], pT_bin=[5,200], d0_max=1, use_FSR=False):
    """Return a 4-tuple of MyMuon objects which pass muon selections in H->ZZ->4mu sample.

    NOTE: This function works for post-UFHZZ4L analyzer, not lite skim.
    
    Parameters
    ----------
    t : ROOT.TTree
        The TTree which holds all data for each event.
    evt_num : int
        Which event in the TTree.
    eta_bin : 2-elem list
        Keep muons with reco eta in range: [eta_min, eta_max]
    pT_bin : 2-elem list
        Keep muons with reco pT in range: [pT_min, pT_max] (GeV)
    d0_max : float
        Keep muons with d0 < d0_max (cm).
    use_FSR : bool
        If True, then store reco pT accounting for FSR.

    Returns
    -------
    If event passes selections: 
        4-tuple of MyMuon objects.
    If event does NOT pass selections: 
        4-tuple of NoneType (None, None, None, None).
    """
    # global n_evts_passed
    bad_muons = (None, None, None, None)
    t.GetEntry(evt_num)
    
    if not passed_Higgs_evt_selection(t):
        return bad_muons

    rec_ndcs_ls = list(t.lep_Hindex)  # Elements that correspond to 4 leptons which build Higgs candidate.
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
    mu_ls = make_muon_ls(t, rec_ndcs_ls, gen_ndcs_ls, use_FSR)
    all_muons_passed = apply_kinem_selections(mu_ls, eta_bin=eta_bin, pT_bin=pT_bin, d0_max=d0_max)
    if not all_muons_passed:
        return bad_muons

    # NOW event is good.
    # n_evts_passed += 1

    # Return the 4 good muons.
    good_muons = tuple(mu_ls)
    assert len(good_muons) == 4
    return good_muons

def build_muons_from_Hmumu_event(t, evt_num, eta_bin=[0,2.4], pT_bin=[20,200], d0_max=1):
    """Return a 2-tuple of MyMuon objects which pass muon selections in H->2mu sample.
    
    Parameters
    ----------
    t : ROOT.TTree
        The TTree which holds all data for each event.
    evt_num : int
        Which event in the TTree.
    eta_bin : 2-elem list
        Keep muons with reco eta in range: [eta_min, eta_max]
    pT_bin : 2-elem list
        Keep muons with reco pT in range: [pT_min, pT_max] (GeV)
    d0_max : float
        Keep muons with d0 < d0_max (cm).

    Returns
    -------
    If event passes selections: 
        2-tuple of MyMuon objects.
    If event does NOT pass selections: 
        2-tuple of NoneType (None, None).
    """
    # global n_evts_passed
    bad_muons = (None, None)
    t.GetEntry(evt_num)
    
    if not passed_Hmumu_evt_selection(t):
        return bad_muons

    # Filippo's MioSkim_2L.C already accounts for lep_genindex.
    # gen_ndcs_ls = get_ndcs_gen(rec_ndcs_ls, lep_genindex_ls)

    # rec_ID_ls = list(t.lep_id)
    # gen_ID_ls = list(t.GENlep_id)
    # if not check_matched_IDs(rec_ndcs_ls, gen_ndcs_ls, rec_ID_ls, gen_ID_ls):
    #     continue
    # if not check_2_OSSF_muon_pairs(gen_ID_ls):
    #     continue
    # if not check_2_OSSF_muon_pairs(rec_ID_ls):
    #     continue

    # Event looks good so far. 
    # Now check kinematics of muons.
    mu_ls = make_muon_ls_fromliteskim(t)
    all_muons_passed = apply_kinem_selections(mu_ls, eta_bin=eta_bin, pT_bin=pT_bin, d0_max=d0_max)
    if not all_muons_passed:
        return bad_muons

    # NOW event is good.
    # n_evts_passed += 1

    # Return the 2 good muons.
    good_muons = tuple(mu_ls)
    assert len(good_muons) == 2
    return good_muons