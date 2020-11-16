import ROOT
import numpy as np

from Utils_Python.Selections import Selector
from Utils_Python.Utils_Physics import calc_dphi, calc_dR  
from d0_Studies.d0_Utils.d0_fns import find_bin_edges_of_value, correct_muon_pT, calc_num_bins, print_header_message

def get_bin_centers_from_TH1(hist):
    """Return a list of the centers of each bin."""
    n_bins = hist.GetNbinsX()
    return [hist.GetBinCenter(x) for x in range(1, n_bins+1)]

def get_bin_vals_from_TH1(hist): 
    """Return a list of the values in each bin."""
    n_bins = hist.GetNbinsX()
    return [hist.GetBinContent(x) for x in range(1,n_bins+1)]
    
def get_bin_edges_from_TH1(hist):
    """
    Return a list of bin edges from a TH1. 
    If the TH1 has N bins, then len(edge_ls) == N + 1.
    """
    n_bins = hist.GetNbinsX()
    edge_ls = [hist.GetBinLowEdge(x) for x in range(1,n_bins+1)]
    rightmost_edge = hist.GetBinLowEdge(n_bins) + hist.GetBinWidth(n_bins)
    edge_ls.append(rightmost_edge)
    return edge_ls

def skip_black_yellow_fit_line_colors(count):
        """
        Skip black and yellow fit lines:
          - Avoid black since data points are usually black
          - Avoid yellow since it is difficult to see in ROOT. 
        This should be used in a for/while loop to avoid black and yellow.
        Use the counter of the loop to shift the colors.

        Returns
        -------
        color : int
            Color of fit line based on which iterative fit you're doing (count).
            Look here for which ints refer to which colors:
            https://root.cern.ch/doc/master/classTColor.html
        """
        if (1 <= count) and (count <= 3):
            # First, second, third fit should shift colors by 1 to avoid black fit line.
            color = count + 1
        else:
            # Fourth fit and beyond should avoid yellow.
            color = count + 2
        return color

def fill_m2mu_hist(tree, hist_m2mu, 
                    sample_name, n_evts,
                    bounds_ls_eta=None, bounds_ls_pT=None, 
                    apply_dR_cut=True, apply_m2l_cut=True):
    """
    Fills a histogram with dimuon invariant mass events.
    Checks to make sure that BOTH muons pass selections specified. 

    Parameters
    ----------
    tree : ROOT.TTree
        Contains the data to be histogrammed.
    hist_m2mu : ROOT.TH1
        Histogram to store dimuon invariant mass.
    sample_name : str
        Given a sample_name, apply all relevent cuts to each event. 
        References a "cuts dictionary".
        "Jpsi", "Upsilon", "DY", "Higgs"
    n_evts : int
        The number of entries to loop over. 
        Not necessarily the number of entries to be found. 
        Use n_evts = -1 to run over all events in TTree.
    bounds_ls_eta : 2-element list
        Specifies the eta range that BOTH muons must satisfy. 
        eta_min < abs(mu.eta) < eta_max
    bounds_ls_pT : 2-element list
        Specifies the pT range that BOTH muons must satisfy. 
        eta_min < mu.pT < eta_max
    apply_dR_cut : bool
        If True, make sure for both muons: mu.dR < dR_cut.
        dR_cut is determined by sample_name and specified in Utils_Selections.
    apply_m2l_cut : bool
        If True, make sure that: m2l_min < m2l_this_event < m2l_max
        m2l_cut is determined by sample_name and specified in Utils_Selections.
    """
    print(f"Scanning {sample_name} events. Filling hists.")

    if bounds_ls_eta is not None:
        check_eta = True
        eta_min = bounds_ls_eta[0]
        eta_max = bounds_ls_eta[1]
    if bounds_ls_pT is not None:
        check_pT = True
        pT_min = bounds_ls_pT[0]
        pT_max = bounds_ls_pT[1]

    if n_evts == -1:
        # Run over all entries.
        n_evts = tree.GetEntries()
    
    # Handles selection criteria, based on sample_name.
    event_selector = Selector(sample_name)

    for evt in range(n_evts):
        tree.GetEntry(evt)

        # Use RECO info to make muons.
        muon1 = ROOT.Math.PtEtaPhiMVector(tree.pT1, tree.eta1, tree.phi1, tree.m1)
        muon2 = ROOT.Math.PtEtaPhiMVector(tree.pT2, tree.eta2, tree.phi2, tree.m2)

        # Form resonance.
        res = muon1 + muon2
        mll = res.M()

        # Get GEN data from event.
        muon1.genLep_pt   = tree.genLep_pt1
        muon2.genLep_pt   = tree.genLep_pt2
        muon1.genLep_eta  = tree.genLep_eta1
        muon2.genLep_eta  = tree.genLep_eta2
        muon1.genLep_phi  = tree.genLep_phi1
        muon2.genLep_phi  = tree.genLep_phi2
        muon1.charge = -1 if tree.Id1 == 13 else +1
        muon2.charge = -1 if tree.Id2 == 13 else +1

        for mu in [muon1, muon2]:
            mu.pT = mu.Pt()
            mu.eta = mu.Eta()
            mu.phi = mu.Phi()
            mu.mass = mu.M()
            # Kinematic differences. 
            mu.deta = mu.eta - mu.genLep_eta
            mu.dphi = calc_dphi(mu.phi, mu.genLep_phi)
            mu.dR   = calc_dR(mu.deta, mu.dphi)
            mu.dpToverpT = (mu.pT - mu.genLep_pt) / mu.genLep_pt

            mu.selection_ls = []
            if apply_dR_cut:
                mu.passed_dR = event_selector.pass_selection_dR(mu.dR)
                mu.selection_ls.append(mu.passed_dR)

            if apply_m2l_cut:
                # Both muons should have the same passed_m2l bool.
                mu.passed_m2l = event_selector.pass_selection_m2l(mll)
                mu.selection_ls.append(mu.passed_m2l)

            if check_eta:
                mu.passed_eta = event_selector.pass_selection_abseta(eta=mu.eta, eta_min=eta_min, eta_max=eta_max)
                mu.selection_ls.append(mu.passed_eta)

            if check_pT:
                mu.passed_pT = event_selector.pass_selection_val(val=mu.pT, val_min=pT_min, val_max=pT_max)
                mu.selection_ls.append(mu.passed_pT)

            mu.passed_all = all(mu.selection_ls)

        # End muon loop. 
        # Fill m2mu hist - but only once. 
        # Only requires that both muons pass dR cut and m2mu cut.
        if (muon1.passed_all) and (muon2.passed_all):
            hist_m2mu.Fill(mll)
    # End evt loop.

def fill_dpToverpT_hist(tree, hist_dpToverpT, 
                    sample_name, n_evts,
                    bounds_ls_eta=None, bounds_ls_pT=None, 
                    apply_dR_cut=True, apply_m2l_cut=True):
    """
    Fills a histogram with the (pT_rec - pT_gen) / pT_gen values 
    of each muon that passes the given selections. 
    Works for dimuon events. 
    Each muon from each event is treated individually.

    Parameters
    ----------
    tree : ROOT.TTree
        Contains the data to be histogrammed.
    hist_dpToverpT : ROOT.TH1
        Histogram to store (pT_rec - pT_gen) / pT_gen values.
    sample_name : str
        Given a sample_name, apply all relevent cuts to each event. 
        References a "cuts dictionary".
        "Jpsi", "Upsilon", "DY", "Higgs"
    n_evts : int
        The number of entries to loop over. 
        Not necessarily the number of entries to be found. 
        Use n_evts = -1 to run over all events in TTree.
    bounds_ls_eta : 2-element list
        Specifies the eta range that BOTH muons must satisfy. 
        eta_min < abs(mu.eta) < eta_max
    bounds_ls_pT : 2-element list
        Specifies the pT range that BOTH muons must satisfy. 
        eta_min < mu.pT < eta_max
    apply_dR_cut : bool
        If True, make sure for both muons: mu.dR < dR_cut.
        dR_cut is determined by sample_name and specified in Utils_Selections.
    apply_m2l_cut : bool
        If True, make sure that: m2l_min < m2l_this_event < m2l_max
        m2l_cut is determined by sample_name and specified in Utils_Selections.
    """
    print(f"Scanning {sample_name} events. Filling hists.")

    if bounds_ls_eta is not None:
        check_eta = True
        eta_min = bounds_ls_eta[0]
        eta_max = bounds_ls_eta[1]
    if bounds_ls_pT is not None:
        check_pT = True
        pT_min = bounds_ls_pT[0]
        pT_max = bounds_ls_pT[1]

    if n_evts == -1:
        # Run over all entries.
        n_evts = tree.GetEntries()
    
    # Handles selection criteria, based on sample_name.
    event_selector = Selector(sample_name)

    for evt in range(n_evts):
        tree.GetEntry(evt)

        # Use RECO info to make muons.
        muon1 = ROOT.Math.PtEtaPhiMVector(tree.pT1, tree.eta1, tree.phi1, tree.m1)
        muon2 = ROOT.Math.PtEtaPhiMVector(tree.pT2, tree.eta2, tree.phi2, tree.m2)

        # Form resonance.
        res = muon1 + muon2
        mll = res.M()

        # Get GEN data from event.
        muon1.genLep_pt   = tree.genLep_pt1
        muon2.genLep_pt   = tree.genLep_pt2
        muon1.genLep_eta  = tree.genLep_eta1
        muon2.genLep_eta  = tree.genLep_eta2
        muon1.genLep_phi  = tree.genLep_phi1
        muon2.genLep_phi  = tree.genLep_phi2
        muon1.charge = -1 if tree.Id1 == 13 else +1
        muon2.charge = -1 if tree.Id2 == 13 else +1

        for mu in [muon1, muon2]:
            mu.pT = mu.Pt()
            mu.eta = mu.Eta()
            mu.phi = mu.Phi()
            mu.mass = mu.M()
            # Kinematic differences. 
            mu.deta = mu.eta - mu.genLep_eta
            mu.dphi = calc_dphi(mu.phi, mu.genLep_phi)
            mu.dR   = calc_dR(mu.deta, mu.dphi)
            mu.dpToverpT = (mu.pT - mu.genLep_pt) / mu.genLep_pt

            mu.selection_ls = []
            if apply_dR_cut:
                mu.passed_dR = event_selector.pass_selection_dR(mu.dR)
                mu.selection_ls.append(mu.passed_dR)

            if apply_m2l_cut:
                # Both muons should have the same passed_m2l bool.
                mu.passed_m2l = event_selector.pass_selection_m2l(mll)
                mu.selection_ls.append(mu.passed_m2l)

            mu.passed_dR_and_m2l = all(mu.selection_ls)
            if not mu.passed_dR_and_m2l:
                # This muon is bad. Skip it.
                continue

            if check_eta:
                mu.passed_eta = event_selector.pass_selection_abseta(eta=mu.eta, eta_min=eta_min, eta_max=eta_max)
                mu.selection_ls.append(mu.passed_eta)

            if check_pT:
                mu.passed_pT = event_selector.pass_selection_val(val=mu.pT, val_min=pT_min, val_max=pT_max)
                mu.selection_ls.append(mu.passed_pT)

            mu.passed_all = all(mu.selection_ls)

            if (mu.passed_all):
                # Fill dpT/pT hist.
                # Put this muon in this hist, only if it satisfies the KinBin criteria.
                # (i.e. correct eta, pT, and dR, m2l.)
                # Treat all muons independently.
                # Could have up to 2x number of events in TTree.
                hist_dpToverpT.Fill(mu.dpToverpT)
        # End muon loop. 
    # End evt loop.

def fill_m2mu_and_dpToverpT_hists(tree, hist_m2mu, hist_dpToverpT, 
                                    sample_name, n_evts,
                                    bounds_ls_eta=None, bounds_ls_pT=None, 
                                    apply_dR_cut=True, apply_m2l_cut=True):
    print(f"Scanning {sample_name} events. Filling hists.")

    if bounds_ls_eta is not None:
        check_eta = True
        eta_min = bounds_ls_eta[0]
        eta_max = bounds_ls_eta[1]
    if bounds_ls_pT is not None:
        check_pT = True
        pT_min = bounds_ls_pT[0]
        pT_max = bounds_ls_pT[1]

    if n_evts == -1:
        # Run over all entries.
        n_evts = tree.GetEntries()
    
    # Handles selection criteria, based on sample_name.
    event_selector = Selector(sample_name)

    for evt in range(n_evts):
        tree.GetEntry(evt)

        # Use RECO info to make muons.
        muon1 = ROOT.Math.PtEtaPhiMVector(tree.pT1, tree.eta1, tree.phi1, tree.m1)
        muon2 = ROOT.Math.PtEtaPhiMVector(tree.pT2, tree.eta2, tree.phi2, tree.m2)

        # Form resonance.
        res = muon1 + muon2
        mll = res.M()

        # Get GEN data from event.
        muon1.genLep_pt   = tree.genLep_pt1
        muon2.genLep_pt   = tree.genLep_pt2
        muon1.genLep_eta  = tree.genLep_eta1
        muon2.genLep_eta  = tree.genLep_eta2
        muon1.genLep_phi  = tree.genLep_phi1
        muon2.genLep_phi  = tree.genLep_phi2
        muon1.charge = -1 if tree.Id1 == 13 else +1
        muon2.charge = -1 if tree.Id2 == 13 else +1

        for mu in [muon1, muon2]:
            mu.pT = mu.Pt()
            mu.eta = mu.Eta()
            mu.phi = mu.Phi()
            mu.mass = mu.M()
            # Kinematic differences. 
            mu.deta = mu.eta - mu.genLep_eta
            mu.dphi = calc_dphi(mu.phi, mu.genLep_phi)
            mu.dR   = calc_dR(mu.deta, mu.dphi)
            mu.dpToverpT = (mu.pT - mu.genLep_pt) / mu.genLep_pt

            mu.selection_ls = []
            if apply_dR_cut:
                mu.passed_dR = event_selector.pass_selection_dR(mu.dR)
                mu.selection_ls.append(mu.passed_dR)

            if apply_m2l_cut:
                # Both muons should have the same passed_m2l bool.
                mu.passed_m2l = event_selector.pass_selection_m2l(mll)
                mu.selection_ls.append(mu.passed_m2l)

            if (apply_dR_cut) and (apply_m2l_cut):
                mu.passed_dR_and_m2l = all(mu.selection_ls)
            else: 
                mu.passed_dR_and_m2l = False

            if check_eta:
                mu.passed_eta = event_selector.pass_selection_abseta(eta=mu.eta, eta_min=eta_min, eta_max=eta_max)
                mu.selection_ls.append(mu.passed_eta)

            if check_pT:
                mu.passed_pT = event_selector.pass_selection_val(val=mu.pT, val_min=pT_min, val_max=pT_max)
                mu.selection_ls.append(mu.passed_pT)

            # See if this muon passed all selections.
            if all(mu.selection_ls):
                # Fill dpT/pT hist.
                # Put this muon in this hist, only if it satisfies the KinBin criteria.
                # (i.e. correct eta, pT, and dR, m2l.)
                # Treat all muons independently.
                # Could have up to 2x number of events in TTree.
                hist_dpToverpT.Fill(mu.dpToverpT)

        # End muon loop. 
        # Fill m2mu hist - but only once. Regardless of KinBin.
        if (muon1.passed_dR_and_m2l) and (muon2.passed_dR_and_m2l):
            hist_m2mu.Fill(mll)
    # End evt loop.

def fill_hists_m2mu_dpToverpT(tree, 
               hist_m2mu_plus, hist_m2mu_minus, hist_dpToverpT, 
               N, sample_name, 
               bounds_ls_eta=None, bounds_ls_pT=None, 
               apply_dR_cut=True, apply_m2l_cut=True):
    """
    FIXME: This function is just like fill_m2mu_and_dpToverpT_hists().
           Possibly merge them into single function?

    Fill three hists with m_mumu events:
      (1) hist_m2mu_plus: dimuon invariant mass hist (good for J/psi, DY, Upsilon, etc.)
        - Only when the mu+ from the event passes selections does it enter this hist.
      (2) hist_m2mu_minus: same as above, but with plus -> minus
      (3) hist_dpToverpT: the (pT_rec - pT_gen) / pT_gen hist

    Parameters
    ----------
    tree : ROOT.TTree
        Contains the data to be histogrammed.
    hist_m2mu_plus, hist_m2mu_minus, hist_dpToverpT : ROOT.TH1
        Explained above.
    N : int
        The number of entries to loop over. 
        Not necessarily the number of entries to be found. 
        Use N = -1 to run over all events in TTree.
    sample_name : str
        Given a sample_name, apply all relevent cuts to each event. 
        References a "cuts dictionary".
    """
    if bounds_ls_eta is not None:
        check_eta = True
        eta_min = bounds_ls_eta[0]
        eta_max = bounds_ls_eta[1]
    if bounds_ls_pT is not None:
        check_pT = True
        pT_min = bounds_ls_pT[0]
        pT_max = bounds_ls_pT[1]
    
    event_selector = Selector(sample_name)

    if N == -1:
        # Run over all entries.
        N = tree.GetEntries()
        
    evts_found = 0
    for evt in range(N):
        tree.GetEntry(evt)
        
        # Use RECO info to make muons.
        muon1 = ROOT.Math.PtEtaPhiMVector(tree.pT1, tree.eta1, tree.phi1, tree.m1)
        muon2 = ROOT.Math.PtEtaPhiMVector(tree.pT2, tree.eta2, tree.phi2, tree.m2)
        
        # Form resonance.
        res = muon1 + muon2
        mll = res.M()
        
        # Get GEN data from event.
        muon1.genLep_pt   = tree.genLep_pt1
        muon2.genLep_pt   = tree.genLep_pt2
        muon1.genLep_eta  = tree.genLep_eta1
        muon2.genLep_eta  = tree.genLep_eta2
        muon1.genLep_phi  = tree.genLep_phi1
        muon2.genLep_phi  = tree.genLep_phi2
        muon1.charge = -1 if tree.Id1 == 13 else +1
        muon2.charge = -1 if tree.Id2 == 13 else +1
        
        for mu in [muon1, muon2]:
            mu.pT = mu.Pt()
            mu.eta = mu.Eta()
            mu.phi = mu.Phi()
            mu.mass = mu.M()
            # Kinematic differences. 
            mu.deta = mu.eta - mu.genLep_eta
            mu.dphi = calc_dphi(mu.phi, mu.genLep_phi)
            mu.dR   = calc_dR(mu.deta, mu.dphi)
            mu.dpToverpT = (mu.pT - mu.genLep_pt) / mu.genLep_pt

            mu.selection_ls = []
            if check_eta:
                mu.passed_eta = event_selector.pass_selection_abseta(eta=mu.eta, eta_min=eta_min, eta_max=eta_max)
                mu.selection_ls.append(mu.passed_eta)

            if check_pT:
                mu.passed_pT = event_selector.pass_selection_val(val=mu.pT, val_min=pT_min, val_max=pT_max)
                mu.selection_ls.append(mu.passed_pT)

            if apply_dR_cut:
                mu.passed_dR = event_selector.pass_selection_dR(mu.dR)
                mu.selection_ls.append(mu.passed_dR)

            if apply_m2l_cut:
                # Both muons should have the same passed_m2l bool.
                mu.passed_m2l = event_selector.pass_selection_m2l(mll)
                mu.selection_ls.append(mu.passed_m2l)
                
            # See if this muon passed all selections.
            mu.passedall = all(mu.selection_ls)

        for mu in [muon1, muon2]:
            if (mu.passedall):
                # Fill the first hist. Treat all muons independently.
                # I.e. this hist could have up to 2x number of events in TTree.
                hist_dpToverpT.Fill(mu.dpToverpT)
                # Fill the other two hists - maybe.
                if mu.charge == -1:
                    hist_m2mu_minus.Fill(mll)
                else:
                    hist_m2mu_plus.Fill(mll)

        # if (muon1.passedall) or (muon2.passedall): evts_found += 1

def make_hist_lookuptable(eta_ls, pT_ls, sample=None, hist_type=None, bin_info_ls=None):#sample_ls, hist_type_ls, ):
    """
    Return a dictionary of ROOT.TH1 objects, one for each eta, pT bin
    specified by User. 

    Parameters
    ----------
    eta_ls : list
        The eta bin edges. 
    pT_ls : list
        The pT bin edges. 
    sample : str
        The sample name, like: "Jpsi" or "Upsilon"
        Used for identification of different histograms. 
    hist_type : str
        Any str to differentiate between types of hists. 
        E.g., "dpToverpT" or "m2mu"

    #--- Below is old.
    sample_ls : list of str
        The sample names, like: ["Jpsi", "Upsilon"]
        Used for identification of different histograms. 
    hist_type_ls : list of str
        The type of histogram you want per eta, per pT, per sample. 
    bin_dict : dict
        Organizes the binning info for each hist. 
        Structure: 
        bin_dict = {
            "Sample_name" (like "Jpsi") : {"hist_type" : [x_min, x_max, bin_width]}
        }
    #--- Above is old.

    Note: All these input parameters are just used to generate keys 
          of the hist_dict.

    Returns
    -------
    hist_dict : dict
        Dictionary of ROOT.TH1 objects, all empty. 

    Notes
    -----
    It seems better to use a single dictionary of histograms instead of 
    nested dictionaries:
      - fewer for loops
      - more readable
      - easier to acccess keys and values
      
    Testing results:
      - 4 nested dicts, 936 histos total, each filled with 20K samples from normal dist:
        Took: 1min 32s
      - 1 dict, 936 histos (all at same level of dict), also filled with 20K samples:
        Took: 1min 13s
        
    """
    from ROOT import TH1F
    
    hist_dict = {}
    for k in range(len(eta_ls)-1):
        eta_min = eta_ls[k]
        eta_max = eta_ls[k+1]
        eta_key = "{}eta{}".format(eta_min, eta_max)

        for j in range(len(pT_ls)-1):
            pT_min = pT_ls[j]
            pT_max = pT_ls[j+1]
            pT_key = "{}pT{}".format(pT_min, pT_max)

            # for sample in sample_ls:
                # for h_type in hist_type_ls:
            if sample is None:
                h_name = "h_{}_{}_{}".format(eta_key, pT_key, hist_type)
            else:
                # Sample specific.
                h_name = "h_{}_{}_{}_{}".format(eta_key, pT_key, sample, hist_type)
            print(f"Making h_name: {h_name}")
            
            # x_min = bin_dict[sample][h_type][0]
            # x_max = bin_dict[sample][h_type][1]
            # bwidth = bin_dict[sample][h_type][2]

            x_min = bin_info_ls[0]
            x_max = bin_info_ls[1]
            bwidth = bin_info_ls[2]
            n_bins = calc_num_bins(*bin_info_ls)
            # n_bins = int(round( (x_max - x_min)/float(bwidth) ))
        
            hist_dict[h_name] = TH1F(h_name, h_name, n_bins, x_min, x_max)
            hist_dict[h_name].Sumw2()
                    
    return hist_dict

def fill_dict_of_dpToverpT_hists(tree, hist_dict, hist_type,
                                 sample, n_evts,
                                 eta_binedge_ls=None, pT_binedge_ls=None, 
                                 apply_dR_cut=True, apply_m2l_cut=True, verbose=False,
                                 apply_d0_pT_corrections=False, pT_corr_factor_dict=None):
    """
    Fills a dictionary of (pT_rec - pT_gen) / pT_gen histograms.
    The key 
    The dpT/pT value is only stored if the muon passes the given selections. 
    Works for dimuon events. 
    Each muon from each event is treated individually.

    Parameters
    ----------
    tree : ROOT.TTree
        Contains the data to be histogrammed.
    hist_dict : dict, or list of dict
        Dicionary of ROOT.TH1 objects. 
        Key : str, 'h_<eta_min>eta<eta_max>_<pT_min>pT<pT_max>_<sample>_dpToverpT'
        Val : The ROOT.TH1 hist to store the dpT/pT data for this (eta, pT) bin.
          Example of a key = 'h_0.2eta0.4_10.0pT14.0_Jpsi_dpToverpT'
        Must be the same length as hist_type.
    hist_type : str, or list of str
        Just a label of what kind of hist(s) to add muons to:
            e.g., "dpToverpT", or ["dpToverpT", "dpToverpTcorr"]
        Must be the same length as hist_dict.
    sample : str
        Given a sample, apply all relevent cuts to each event. 
        References a "cuts dictionary".
        "Jpsi", "Upsilon", "DY", "Higgs"
    n_evts : int
        The number of entries to loop over. 
        Not necessarily the number of entries to be found. 
        Use n_evts = -1 to run over all events in TTree.
    eta_binedge_ls : list
        Specifies the abs(eta) bin edges.
        Only muons that fall within one of these eta bins will be kept.
    pT_binedge_ls : list
        Specifies the pT bin edges.
        Only muons that fall within one of these pT bins will be kept.
    apply_dR_cut : bool
        If True, make sure for both muons: mu.dR < dR_cut.
        dR_cut is determined by sample and specified in Utils_Selections.
    apply_m2l_cut : bool
        If True, make sure that: m2l_min < m2l_this_event < m2l_max
        m2l_cut is determined by sample and specified in Utils_Selections.
    verbose : bool
        Print debug info.
    """
    # Quick checks.
    if (apply_d0_pT_corrections):
        assert isinstance(hist_dict, list)
        assert isinstance(hist_type, list)
        assert len(hist_dict) == len(hist_type)
        print_header_message("Applying correction factors, behbeh")

        someone_is_bad = any([x is None for x in (eta_binedge_ls, pT_binedge_ls, pT_corr_factor_dict)])
        assert not someone_is_bad
        
    event_selector = Selector(sample)
    
    if n_evts == -1:
        n_evts = tree.GetEntries()
    
    for evt in range(n_evts):
        tree.GetEntry(evt)
        if (evt % 500000) == 0: 
            print(f"...Running over event {evt}")
        # Use RECO info to make muons.
        muon1 = ROOT.Math.PtEtaPhiMVector(tree.pT1, tree.eta1, tree.phi1, tree.m1)
        muon2 = ROOT.Math.PtEtaPhiMVector(tree.pT2, tree.eta2, tree.phi2, tree.m2)

        # Form resonance.
        res = muon1 + muon2
        mll = res.M()

        # Get GEN data from event.
        muon1.genLep_pt   = tree.genLep_pt1
        muon2.genLep_pt   = tree.genLep_pt2
        muon1.genLep_eta  = tree.genLep_eta1
        muon2.genLep_eta  = tree.genLep_eta2
        muon1.genLep_phi  = tree.genLep_phi1
        muon2.genLep_phi  = tree.genLep_phi2
        muon1.d0BS        = tree.d0BS1
        muon2.d0BS        = tree.d0BS2
        muon1.charge = -1 if tree.Id1 == 13 else +1
        muon2.charge = -1 if tree.Id2 == 13 else +1

        if (verbose): 
            print("evt:", evt)
        for count,mu in enumerate([muon1, muon2], 1):
            mu.pT = mu.Pt()
            mu.eta = mu.Eta()
            mu.phi = mu.Phi()
            mu.mass = mu.M()
            # Kinematic differences. 
            mu.deta = mu.eta - mu.genLep_eta
            mu.dphi = calc_dphi(mu.phi, mu.genLep_phi)
            mu.dR   = calc_dR(mu.deta, mu.dphi)
            mu.dpToverpT = (mu.pT - mu.genLep_pt) / mu.genLep_pt

            # In the very unlikely case that eta,pT of mu EXACTLY falls
            # on one of the bin edges (it's happened twice now)...,

            mu.selection_ls = []
            # First check m2l and dR cuts. If these fail, then the muon is useless.
            if apply_dR_cut:
                mu.passed_dR = event_selector.pass_selection_dR(mu.dR)
                mu.selection_ls.append(mu.passed_dR)

            if apply_m2l_cut:
                # Both muons should have the same passed_m2l bool.
                mu.passed_m2l = event_selector.pass_selection_m2l(mll)
                mu.selection_ls.append(mu.passed_m2l)

            mu.passed_dR_and_m2l = all(mu.selection_ls)
            if not mu.passed_dR_and_m2l:
                # This muon is bad. Skip it.
                continue

            # Check this muon's eta and pT and see where it should fall in the given lists.
            binedge_min_eta, binedge_max_eta = find_bin_edges_of_value(abs(mu.eta), np.array(eta_binedge_ls))
            binedge_min_pT,  binedge_max_pT  = find_bin_edges_of_value(mu.pT, np.array(pT_binedge_ls))

            # Are any of them bad? (i.e. the pT or eta was outside the bin range)
            binls = [binedge_min_eta, binedge_max_eta, binedge_min_pT, binedge_max_pT]
            outside = any(x is None for x in binls)
            if (verbose):
                info = "  eta_bin=[{}, {}], pT_bin=[{}, {}]".format(binedge_min_eta, binedge_max_eta, binedge_min_pT, binedge_max_pT)
                print("  muon {} belongs to: {}".format(count, info))
                print("  passed_dR_and_m2l = {}".format(mu.passed_dR_and_m2l))
                print("            outside = {}".format(outside))

            if (outside): 
                continue
                
            # If you made it here, you have a good muon!
            # Put this muon in the hist where it belongs. 
            eta_key = "{}eta{}".format(binedge_min_eta, binedge_max_eta)
            pT_key = "{}pT{}".format(binedge_min_pT, binedge_max_pT)
            
            # See if we should apply corrections from d0 studies. 
            hist_key_ls = []
            if (apply_d0_pT_corrections):
                #Some ad hoc shit right here.
                for ht in hist_type:
                    dude = "h_{}_{}_{}_{}".format(eta_key, pT_key, sample, ht)
                    hist_key_ls.append(dude)

                mu.pT_corr = correct_muon_pT(mu.eta, mu.pT, mu.charge, mu.d0BS, 
                                            pT_corr_factor_dict,
                                            eta_binedge_ls=eta_binedge_ls, 
                                            pT_binedge_ls=pT_binedge_ls,
                                            verbose=verbose)
                mu.dpToverpT_corr = (mu.pT_corr - mu.genLep_pt) / mu.genLep_pt
                
                dptdivpt_ls = [mu.dpToverpT, mu.dpToverpT_corr]

            # Finally, fill the hist(s).
            try:
                if (apply_d0_pT_corrections):
                    # if (evt % 2000) == 0: print("Doing weird triple zipper...")
                    for hdict,hkey,dptdivpt in zip(hist_dict,hist_key_ls,dptdivpt_ls):
                        hdict[hkey].Fill(dptdivpt)
                else:
                    hist_dict[hist_key].Fill(this_dpToverpT)
            except KeyError:
                print("mu.pT", mu.pT)
                print("mu.eta ", mu.eta )
                print("eta_key ", eta_key )
                print("pT_key", pT_key)
                print("hist_key,", hist_key)
                print("hist_dict:\n",hist_dict)
                raise KeyError
        # End muon loop for this event. 
        assert muon1.passed_m2l == muon2.passed_m2l
    # End evt loop.