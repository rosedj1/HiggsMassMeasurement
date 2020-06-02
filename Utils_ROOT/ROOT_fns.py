import ROOT
from Utils_Python.Utils_Selections import Selector
from Utils_Python.Utils_Physics import calc_dphi, calc_dR  
from d0_Studies.d0_Utils.d0_fns import find_bin_edges_of_value

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
    Fills a histogram with dimuon invariant mass events.
    Checks to make sure that BOTH muons pass selections specified. 

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