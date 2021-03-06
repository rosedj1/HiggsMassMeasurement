import os
import time
import pickle
import threading
import ROOT
from array import array
from pprint import pprint
# Local imports.
from Utils_Python.Selections import (build_muons_from_HZZ4mu_event,
                                     build_muons_from_Hmumu_event,
                                     build_muons_from_DY_event)
from Utils_Python.Utils_Physics import calc_Hmass, perc_diff
from Utils_Python.Utils_Files import check_overwrite, save_to_pkl, make_dirs
from Utils_Python.printing import print_header_message
from Utils_ROOT.ROOT_classes import make_TH1F
from d0_Studies.d0_Utils.d0_cls import KinBin2D
from d0_Studies.d0_Utils.d0_fns import (correct_muon_pT, parse_etapT_key,
                                       get_binedges_from_keys, make_key_from_binedges,
                                       find_bin_edges_of_value)
from d0_Studies.d0_Utils.d0_dicts import process_dct, color_dict_RooFit

def check_n_muons(prod_mode, muon_ls, per_event=True):
    """Ensure that muon_ls has correct number of muons, based on prod_mode.
    
    prod_mode : str
        "DY2mu", "DY2e", "H2mu", "H2e", "H4mu", "H4e"
    muon_ls : list
        List of MyMuons.
    per_event : bool
        If True, make sure final state has exactly expected number of muons
        on a per-event basis. Example: DY2e should have exactly 2 muons.
        If False, make sure that the len(muon_ls) has an integer number
        of events, based on the final state of the prod_mode. 
    """
    n = len(muon_ls)
    assert n > 0, "[ERROR] No muons were found!"
    if prod_mode in ["DY2mu", "DY2e", "H2mu", "H2e"]:
        assert n == 2 if per_event else (n % 2) == 0, f"n_muons_total={n}"
    elif prod_mode in ["H4mu", "H4e"]:
        assert n == 4 if per_event else (n % 4 == 0), f"n_muons_total={n}"
    else:
        raise ValueError("Unknown prod_mode.")

def build_muons_from_process(prod_mode, evt, evt_num,
                             eta_bin=[0, 2.4],
                             pT_bin=[5, 200],
                             d0_bin=[0, 1],
                             inv_m_bin=[105, 140],
                             dR_max=0.002,
                             verbose=False):
    """Return a tuple of MyMuons for 1 event based on prod_mode and selections.
    
    prod_mode : str
        "DY2mu", "DY2e", "H2mu", "H2e", "H4mu", "H4e"
    FIXME:
    - Update docstring.
    - Associate the m2l, m4l, event number to each MyMuon.
    """
    if prod_mode in ["DY2mu", "DY2e"]:
        mu_tup = build_muons_from_DY_event(evt, evt_num,
                                        eta_bin=eta_bin,
                                        pT_bin=pT_bin,
                                        d0_bin=d0_bin,
                                        inv_m_bin=inv_m_bin,
                                        dR_max=dR_max,
                                        verbose=verbose)
    elif prod_mode in ["H4mu", "H4e"]:
        # Check if this event passes Higgs selections.
        # Build the 4 MyMuons from the event.
        mu_tup = build_muons_from_HZZ4mu_event(evt, evt_num,
                                        eta_bin=eta_bin,
                                        pT_bin=pT_bin,
                                        d0_bin=d0_bin,
                                        inv_m_bin=inv_m_bin,
                                        dR_max=dR_max,
                                        verbose=verbose)

    elif prod_mode in ["H2mu", "H2e"]:
        # Check if this event passes H->2mu selections.
        # Build the 2 MyMuons from the event.
        # mu_tup = build_muons_from_Hmumu_event(t, evt_num,
        mu_tup = build_muons_from_Hmumu_event(evt, evt_num,
                                        eta_bin=eta_bin,
                                        pT_bin=pT_bin,
                                        d0_bin=d0_bin,
                                        inv_m_bin=inv_m_bin,
                                        dR_max=dR_max,
                                        verbose=verbose)
    else:
        msg = f"""`prod_mode` ({prod_mode}) must be in:\n["DY2mu", "DY2e", "H4mu", "H4e", "H2mu", "H2e"]"""
        raise ValueError(msg)
    return mu_tup

class MyMuonCollection:
    """A class to handle a list of MyMuon objects."""
    
    def __init__(self, prod_mode):
        self.prod_mode = prod_mode  # Process from which muons were produced.
        # Works for: "DY2mu", "DY2e", "H4mu", "H4e", "H2mu", "H2e"
        self.prod_mode_ls = [prod_mode]  # List of strings. Currently not used.
        self.muon_ls = []
        self.m4mu_ls = []
        self.m4mu_corr_ls = []
        # self.m4mu_corr_geofit_ls = []

        self.KinBin2D_dict = {}
        self.pT_corr_factor_dict = None
        self.do_mu_pT_corr = False   # Will turn True when pT corr is performed.

        self.hist_inclusive_ls = []  # All muons in same plot.
        self.kinbin2d_graph_ls = []  # dpTOverpT vs. qd0 plots.
        self.kinbin2d_hist_ls = []  # qd0, dpT/pT hists.
        self.kinbin3d_iterfitplot_ls = []   
        self.kinbin3d_hist_ls = []   # h_qd0
        # self.plot_ls = []            # Specific KinBin2D plots.
        # self.kinbin2d_plot_ls = []   #
        # self.kinbin3d_plot_ls = []   # KinBin3D distributions, like qd0 and dpTOverpT.

        self.multigraph_ls = []
        self.multigraph_1divpT_ls = []
        self.multigraph_leg_ls = []
        self.multigraph_1divpT_leg_ls = []
        self.multigraph_avgOf1divpT_ls = []
        self.multigraph_avgOf1divpT_leg_ls = []
        self.multigraph_muOf1divpT_ls = []
        self.multigraph_muOf1divpT_leg_ls = []

    def get_prod_modes_str(self):
        """Return a title string of all the production modes of the muons."""
        return ', '.join(self.prod_mode_ls)

    def extract_muons(self, infile_path, prod_mode,
                      n_evts=-1, n_evt_beg=None, n_evt_end=None,
                      print_out_every=10000,
                      eta_lim=[0.0, 2.4],
                      pT_lim=[5, 1000],
                      d0_lim=[0, 1],
                      inv_m_lim=[60.0, 120.0], 
                      dR_max=None,
                      do_mu_pT_corr=False,
                      force_zero_intercept=False,
                      pT_corr_factor_dict=None,
                      correction_type=None,
                      use_GeoFit_algo=False,
                      verbose=False):
        """Fill self.muon_ls with all MyMuon muons from `infile_path`
        which pass selections in `prod_mode`. Also fills self.m4mu_ls per event.

        TODO: Update docstring.
        
        NOTE: 
        - self.m4mu_ls will also become self.m2mu_ls if prod_mode has 2 lep
        in final state.
        - If do_mu_pT_corr, then use pT_corr_factor_dict to
        correct the per muon pT.

        Parameters
        ----------
        infile_path : str
            Absolute file path to root file.
        prod_mode : str
            "DY2mu", "DY2e", "H2mu", "H2e", "H4mu", "H4e"
        n_evts : int
            Number of events to scan over.
            Default of `-1` runs over all events.
        n_evt_beg : int
            Initial INCLUSIVE event number to run over.
        n_evt_end : int
            Final INCLUSIVE event number to run over.
        print_out_every : int
            Show which event number is being processed in batches
            of `print_out_every`.
        do_mu_pT_corr : bool
            If True, perform muon pT correction based on ad hoc d0 studies.
            If no pT_corr_factor_dict is given, will default to
            self.pT_corr_factor_dict.
        pT_corr_factor_dict : dict, optional
            {
                'dpTOverpT_vs_qd0' : {
                    '0.0eta0.2_10.0pT14.0' : {
                        'intercept'     : -0.0004421981787913626,
                        'intercept_err' : 4.7910538089331156e-05,
                        'slope'         : 1.1863059538824643,
                        'slope_err'     : 0.02705983833855872
                    },
                    '0.0eta0.2_100.0pT150.0' : {
                        'intercept'     : ...,
                    },
                    ...,
                },
                'dpTOverpTscaled_vs_qd0' : {
                    '0.0eta0.2_10.0pT14.0' : {
                        ...
                    },
                    ...,
                }
                'dpTOverpTtimesavgOf1divpT_vs_qd0' : {
                    '0.0eta0.2_10.0pT14.0' : {
                        ...
                    },
                    ...,
                }
                'dpTOverpTtimesmuOf1divpT_vs_qd0' : {
                    '0.0eta0.2_10.0pT14.0' : {
                        ...
                    },
                    ...,
                }
            }
        correction_type : str
            "dpTOverpT_vs_qd0"
            "dpTOverpTscaled_vs_qd0"
            "dpTOverpTtimesavgOf1divpT_vs_qd0"
            "dpTOverpTtimesmuOf1divpT_vs_qd0"
        verbose : bool
            Print juicy debug info.
        """
        print(f"[INFO] Extracting muons from:\n{infile_path}")

        specified_max_evts = (n_evts == -1)
        specified_range = (n_evt_beg is not None and n_evt_end is not None)
        assert not (specified_max_evts and specified_range), "Specify all events or range of events."

        if do_mu_pT_corr:
            print(f"...Correcting muon pT using supplied correction factors.")
            self.do_mu_pT_corr = True
            # Make sure only one pT corr factor dict is known.
            if isinstance(pT_corr_factor_dict, dict) and isinstance(self.pT_corr_factor_dict, dict):
                msg = "MyMuonCollection already has a pT_corr_factor_dict, but you are supplying another."
                raise ValueError(msg)
            # If one was not supplied, use the one that already exists.
            if pT_corr_factor_dict is None:
                pT_corr_factor_dict = self.pT_corr_factor_dict[correction_type]
            if use_GeoFit_algo:
                msg = "Using GeoFit Correction Algorithm"
                print_header_message(msg, pad_char="@", n_center_pad_chars=5)
                
        print(f"[INFO] ROOT file opened.")
        f = ROOT.TFile(infile_path, "read")
        dirname = list(f.GetListOfKeys())[0].GetName()
        if "Ana" in dirname:
            t = f.Get("Ana/passedEvents")
        elif "passedEvents" in dirname:
            t = f.Get("passedEvents")
        else:
            raise ValueError("Could not attach to TTree.")

        all_evts = t.GetEntries()
        if n_evts == -1:
            n_evt_beg = 0
            n_evt_end = all_evts - 1
            n_requested_evts = all_evts
            msg = f"Running over ALL ({all_evts}) events."
            print_header_message(msg, pad_char="-", n_center_pad_chars=5)
        elif n_evt_beg is None and n_evt_end is None:
            # User specified n_evts < max.
            n_evt_beg = 0
            n_evt_end = n_evts - 1
            n_requested_evts = n_evts
            msg = f"Running over the first {n_requested_evts} events."
            print_header_message(msg, pad_char="-", n_center_pad_chars=5)
        else:
            # User specified range of events.
            assert n_evt_beg is not None and n_evt_end is not None
            n_requested_evts = n_evt_end - n_evt_beg + 1
            print(f"Running over the specified range of event numbers: {n_evt_beg} -> {n_evt_end}")
            
        # print(f"...Finding {n_requested_evts} events that pass selections...")
        print(f"...Running over events: ({n_evt_beg} -> {n_evt_end})")

        # Event loop.
        n_good_evts = 0
        time_start = time.perf_counter()
        # for evt_num, evt in enumerate(t):
        for evt_num in range(n_evt_beg, n_evt_end + 1):
            t.GetEntry(evt_num)
            if evt_num % print_out_every == 0:
                time_end = time.perf_counter()
                dt = time_end - time_start
                print(f"Running over evt: {evt_num}.")
                print(f"  Time to process last {print_out_every} events: {dt:.6f} s")
                time_start = time.perf_counter()
            mu_tup = build_muons_from_process(prod_mode, t, evt_num,
                                              eta_bin=eta_lim,
                                              pT_bin=pT_lim,
                                              d0_bin=d0_lim,
                                              inv_m_bin=inv_m_lim,
                                              dR_max=dR_max, verbose=verbose)
            if None in mu_tup: 
                continue

            if do_mu_pT_corr:
                # Correct each muon's pT according to the
                # pT correction factors given.
                # Save the muons with old and new pTs.
                # corr_mu_ls = []
                lorentzvec_corr_mu_withFSR_ls = []
                # corr_mu_geofit_ls = []
                for mu in mu_tup:
                    # For all muons in this event, correct the muon pT.
                    # Then evaluate the m4mu_corr for this event.
                    # NOTE: Depending on prod_mode
                    # Correct muon pT WITHOUT accounting for FSR,
                    # then add in FSR later!
                    mu.pT_corr = correct_muon_pT(
                        mu.eta, mu.pT, mu.charge, mu.d0, 
                        # mu.eta, mu.pT_withFSR, mu.charge, mu.d0, 
                        pT_corr_factor_dict, detection="auto",
                        force_zero_intercept=force_zero_intercept,
                        use_GeoFit_algo=use_GeoFit_algo,
                        print_all_muon_info=False,
                        verbose=verbose)
                    # lorentzvec_mu_corr = ROOT.Math.PtEtaPhiMVector(mu.pT_corr, mu.eta, mu.phi, mu.mass)
                    # corr_mu_ls.append(lorentzvec_mu_corr)
                # Now combine any FSR photon to its muon with corrected pT.
                    
                    #--- Clever way to rebuild the FSR photon:
                    # Idea: mu_reco + FSR_photon = mu_withFSR
                    # ==> FSR_photon = mu_withFSR - mu_reco
                    lvec_photon = mu.get_LorentzVector(kind="withFSR") - mu.get_LorentzVector(kind="reco")
                    # Then: FSR_photon + mu_recopTcorr = mu_pTcorrwithFSR
                    lvec_mu_corrpT_FSR = lvec_photon + ROOT.Math.PtEtaPhiMVector(mu.pT_corr, mu.eta, mu.phi, mu.mass)
                    lorentzvec_corr_mu_withFSR_ls.append(lvec_mu_corrpT_FSR)
                #--- End loop over muons.

                # check_n_muons(prod_mode, corr_mu_ls, per_event=True)
                check_n_muons(prod_mode, mu_tup, per_event=True)

                # The below is not tested, but should work!
                # If you have problems, uncomment some of the lines below.
                inv_m_corr = calc_Hmass(lorentzvec_corr_mu_withFSR_ls)  # Can accommodate different num fs muons.

                # for mu in mu_tup:
                #     # An attempt to get rid of self.m4mu_corr_ls.
                #     mu.inv_m_event_corr = inv_m_corr
                # m4mu_corr = calc_Hmass(corr_mu_ls)  # Can accommodate different num fs muons.
                self.m4mu_corr_ls.append(inv_m_corr)

            # Save the inv_mass(muons) info.
            lorentzvec_mu_withFSR_ls = [mu.get_LorentzVector("withFSR") for mu in mu_tup]
            inv_m = calc_Hmass(lorentzvec_mu_withFSR_ls)

            for mu in mu_tup:
                # An attempt to get rid of self.m4mu_ls and self.m4mu_corr_ls.
                mu.inv_m_event = inv_m
                if do_mu_pT_corr:
                    mu.inv_m_event_corr = inv_m_corr
            # lorentzvector_mu_ls = [mu.get_LorentzVector("reco") for mu in mu_tup]
            # m4mu = calc_Hmass(lorentzvector_mu_ls)
            self.m4mu_ls.append(inv_m)
            # Add the good muons to the final muon list.
            self.muon_ls.extend(mu_tup)

            # Some checks to make sure FSR was computed correctly.
            if verbose and "4" in prod_mode:
                percent_diff = perc_diff(inv_m, t.mass4l)
                if (percent_diff > 2):
                    print(
                        f"[WARNING] (calc_inv_m - t.mass4l) / t.mass4l * 100.0 > 2% spotted!\n"
                        f"  Event #{evt_num}\n"
                        f"  nFSRPhotons   = {t.nFSRPhotons}\n"
                        f"  t.mass4l      = {t.mass4l}\n"
                        f"  calc_inv_m    = {inv_m}\n"
                        f"  percent_diff  = {percent_diff}\n"
                        )
                if do_mu_pT_corr:
                    percent_diff = perc_diff(inv_m_corr, t.mass4l)
                    if (percent_diff > 2):
                        print(
                            f"[WARNING] (inv_m - t.mass4l) / t.mass4l * 100.0 > 2% spotted!\n"
                            f"  Event #{evt_num}\n"
                            f"  nFSRPhotons   = {t.nFSRPhotons}\n"
                            f"  t.mass4l      = {t.mass4l}\n"
                            f"  inv_m_corr    = {inv_m_corr}\n"
                            f"  percent_diff  = {percent_diff}\n"
                            )
            n_good_evts += 1
            # if n_good_evts >= n_requested_evts:
            # if n_evt + 1 >= n_requested_evts:
            #     break
        # End evt loop.
        check_n_muons(prod_mode, self.muon_ls, per_event=False)
        if prod_mode in ["H2mu", "H2e"]:
            self.m2mu_ls = self.m4mu_ls
            self.m2mu_corr_ls = self.m4mu_corr_ls
        if do_mu_pT_corr:
            assert len(self.m4mu_ls) == len(self.m4mu_corr_ls)

    def sort_muons(self, eta_ls=None, pT_ls=None,
                   pT_corr_factor_dict=None,
                   n_bins_dpTOverpT=100, x_lim_dpTOverpT=[-0.4, 0.4],
                   n_bins_qd0=100, x_lim_qd0=[-0.01, 0.01], verbose=False):
        """Go through list of muons and put each muon
        into its corresponding KinBin2D.

        A check is performed to see if self.KinBin2D_dict exists.
        If it does then put muons from self.muon_ls into that dict.
        Otherwise make the dict using specified eta_ls and pT_ls.
        """
        if (eta_ls is None) or (pT_ls is None):
            print(
                f"[INFO] No eta_ls or pT_ls detected.\n"
                f"[INFO] Making sure self.KinBin2D_dict has been prepared..."
                )
            assert len(self.KinBin2D_dict) > 0, "No KinBin2D_dict found."
        else:
            print(
                f"[INFO] Creating a KinBin2D for each (eta, pT) bin\n"
                f"       using these provided bin lists:\n"
                f"  eta bins: {eta_ls}\n"
                f"  pT bins:  {pT_ls}\n"
                )

            self.create_KinBins(eta_ls, pT_ls)
        print("[INFO] ...Placing muons into KinBin2Ds...")
        self.place_muons_into_KinBins(eta_ls=eta_ls, pT_ls=pT_ls,
                                      pT_corr_factor_dict=pT_corr_factor_dict,
                                      verbose=verbose)
        print("[INFO] ...Making and filling hists for all KinBin2Ds...")
        # self.make_empty_KinBin_hists()
        self.make_KinBin_hists(n_bins_dpTOverpT=100, x_lim_dpTOverpT=[-0.3, 0.3],
                                n_bins_qd0=100, x_lim_qd0=[-0.01, 0.01])

    def create_KinBins(self, eta_ls, pT_ls):
        """Create a dict of KinBin2Ds, one for each (eta, pT) bin specified in eta_ls, pT_ls.
        
        Example: 
            eta_ls = [0.0, 0.2, 0.4]
            pT_ls  = [5.0, 7.0, 10.0]
            - NOTE: these lists do not need to be the same length!
            
            => Will make 4 KinBin2Ds each with an (eta, pT) bin: 
                {
                    "0.0eta0.2_5.0pT7.0"  : KinBin2D(0.0-0.2, 5.0-7.0),
                    "0.0eta0.2_7.0pT10.0" : KinBin2D(0.0-0.2, 7.0-10.0), 
                    "0.2eta0.4_5.0pT7.0"  : KinBin2D(0.2-0.4, 5.0-7.0),
                    "0.2eta0.4_7.0pT10.0" : KinBin2D(0.2-0.4, 7.0-10.0),
                }
        """
        for eta_min, eta_max in zip(eta_ls[:-1], eta_ls[1:]):
            for pT_min, pT_max in zip(pT_ls[:-1], pT_ls[1:]):
                key = self.make_bin_key(eta_min, eta_max, pT_min, pT_max, title_friendly=False)
                this_KinBin = KinBin2D(eta_range=[eta_min,eta_max], pT_range=[pT_min,pT_max])
                self.KinBin2D_dict[key] = this_KinBin
        print(f"[INFO] All KinBin2Ds have been made from given eta_ls and pT_ls:")
        print(f"  eta : {eta_ls}")
        print(f"   pT : {pT_ls}")

    def create_KinBins_from_pT_corr_dict(self, pT_corr_factor_dict):
        """Create a dict of KinBin2Ds, one for each key in pT_corr_dict.
        
        Example key: '0.0eta0.2_5.0pT7.0'
        Produces KinBin2D(eta_range=[0.0, 0.2], pT_range=[5.0, 7.0])
        """
        err_msg = "You are trying to make a KinBin2D_dict, when one already exists."
        assert len(self.KinBin2D_dict) == 0, err_msg

        for key in pT_corr_factor_dict.keys():
            eta_min, eta_max, pT_min, pT_max = parse_etapT_key(key)
            this_KinBin = KinBin2D(eta_range=[eta_min,eta_max],
                                   pT_range=[pT_min,pT_max])
            self.KinBin2D_dict[key] = this_KinBin
        print(f"[INFO] All KinBin2Ds have been made from pT_corr_dict keys.")

    def make_KinBin3Ds(self, regions, min_muons_per_qd0_bin=1000, verbose=False):
        """Each KinBin2D splits its muons up among a number of equal-entry KinBin3Ds equal to `regions`.
        
        regions : int
            The number of equal-entry KinBin3Ds to sort muons into.
        """
        print(f"...Splitting each KinBin2D up into {regions} KinBin3Ds...")
        for kb2d in self.KinBin2D_dict.values():
            kb2d.make_empty_equalentry_KinBin3Ds(regions, algo=("at_least", min_muons_per_qd0_bin), 
                                                 verbose=verbose, title_friendly=False)
            kb2d.store_muon_info_in_KinBin3Ds(title_friendly=False, verbose=verbose)

    def make_KinBin_hists(self, n_bins_dpTOverpT=100, x_lim_dpTOverpT=[-0.12, 0.12],
                                n_bins_qd0=100, x_lim_qd0=[-0.01, 0.01]):
        """Have each KinBin2D fill and store its own self.h_dpTOverpT and
        self.h_qd0 hists.
        Append the hists to muon_collection.kinbin2d_hist_ls.
        
        NOTE:
        - Plots get stored like, e.g.: self.h_dpTOverpT, self.h_qd0
        - Used to be called: make_empty_KinBin_hists()
        """
        for kb in self.KinBin2D_dict.values():
            kb.make_dpTOverpT_hist(n_bins=n_bins_dpTOverpT, x_lim=x_lim_dpTOverpT)
            kb.make_qd0_hist(n_bins=n_bins_qd0, x_lim=x_lim_qd0)
        self.kinbin2d_hist_ls.extend([kb.h_dpTOverpT for kb in self.KinBin2D_dict.values()])
        self.kinbin2d_hist_ls.extend([kb.h_qd0 for kb in self.KinBin2D_dict.values()])
        print("Done making KinBin hists")

    def place_single_muon_into_KinBin(self, muon, eta_ls, pT_ls, verbose=False):
        """Put MyMuon object into the correct KinBin2D, based on muon's (eta, pT) value.
        
        Parameters
        ----------
        muon : MyMuon obj.
        eta_ls : list of eta bin edges.
        pT_ls : list of pT bin edges.
        verbose : bool for debugging.
        """
        eta_min, eta_max = find_bin_edges_of_value(abs(muon.eta), eta_ls, verbose)
        pT_min, pT_max = find_bin_edges_of_value(muon.pT, pT_ls, verbose)
        if any([x is None for x in (eta_min, eta_max, pT_min, pT_max)]):
            if verbose:
                msg = f"[WARNING] Muon found with strange (eta, pT) values: eta={[eta_min, eta_max]}, pT={[pT_min, pT_max]}"
                print(msg)
                print("Skipping this muon, since your bins cannot hold it!")
            # n_muons_skipped += 1
        else:
            # Put this muon into correct KinBin2D.
            key = self.make_bin_key(eta_min, eta_max, pT_min, pT_max, title_friendly=False)
            self.KinBin2D_dict[key].add_muon(muon)

    def place_muons_into_KinBins(self, eta_ls=None, pT_ls=None, 
                                 pT_corr_factor_dict=None, auto_detect=False, verbose=False):
        """Put muons from muon_ls into correct KinBin2D, based on muon's (eta, pT) value.
        
        After muons have been assigned to KB2D, self.muon_ls is cleared.

        Parameters
        ----------
        auto_detect : bool
            If True, then will auto-detect the KinBin2Ds according to the keys
            of pT_corr_factor_dict.
        """
        if verbose: 
            print(f"...Putting {len(self.muon_ls)} muons into their respective KinBin2Ds...")
        # Store muons in KinBins.
        if (eta_ls is None) or (pT_ls is None) or auto_detect:
            # Automatically decide into which KinBin2D to put muons
            # based on the pT_corr_factor_dict keys.
            assert pT_corr_factor_dict is not None
            for mu in self.muon_ls:
                binedge_tup = get_binedges_from_keys(mu.eta, mu.pT, pT_corr_factor_dict)
                key = make_key_from_binedges(binedge_tup)
                self.KinBin2D_dict[key].add_muon(mu)
        else:
            # Put muons into KB2Ds using provided eta and pT bin edges.
            for ct, muon in enumerate(self.muon_ls):
                if (ct % 1E6) == 0:
                    print(f"  Placed muon #{ct}.")
                self.place_single_muon_into_KinBin(muon, eta_ls, pT_ls, verbose=verbose)
            # End muon loop.
            if verbose:
                n_mu = len(self.muon_ls)
                n_muons_in_KB2Ds = sum([len(kb2d.muon_ls) for kb2d in self.KinBin2D_dict.values()])
                n_muons_skipped = n_mu - n_muons_in_KB2Ds
                print(f"Skipped {n_muons_skipped}/{n_mu} muons ({n_muons_skipped/n_mu * 100.0:.3f})%\n")
        print("[INFO] All muons assigned to KinBin2Ds.\n")
        print_header_message("Overwriting MyMuonCollection.muon_ls to save space.", pad_char="*")
        self.muon_ls = "overwritten"

    def make_bin_key(self, eta_min, eta_max, pT_min, pT_max, title_friendly=False):
        """Return a string which identifies a 2D (eta, pT) bin.
        
        Optional parameter: title_friendly 
            Useful for making file names and title names. 
            If True, replaces '.' with 'p'.
        """
        key = f"{eta_min}eta{eta_max}_{pT_min}pT{pT_max}"
        return key.replace(".", "p") if title_friendly else key

    def do_2D_iter_gaus_fits(self, bins_dpTOverpT=100, bins_qd0=100,
                       x_lim_dpTOverpT=[-0.4,0.4], x_lim_qd0=[-0.01,0.01],
                       fit_whole_range_first_iter=True,
                       iters=1, num_sigmas=2,
                       switch_to_binned_fit=2000, verbose=False, alarm_level="warning",
                       use_mu_pT_corr=False, only_draw_last=False):
        """Perform the iterated Gaussian fit algorithm on each KB2D in
        self.KinBin2D_dict. Record the statistics of the fits.

        If use_mu_pT_corr is True, then another set of iter Gauss fits
        is performed on the muon.pT_corr data.
        """
        if verbose:
            # Make a dict to record the KinBin2Ds which have been processed.
            kb2d_checkdct = {key : "" for key in self.KinBin2D_dict.keys()}
            print("...Performing iterated fits on the following KinBin2Ds:")
            pprint(kb2d_checkdct)
        
        # arg_tup = (
        #     bins_dpTOverpT, bins_qd0,
        #     x_lim_dpTOverpT, x_lim_qd0,
        #     fit_whole_range_first_iter,
        #     iters, num_sigmas,
        #     switch_to_binned_fit,
        #     verbose, alarm_level,
        #     use_mu_pT_corr, only_draw_last
        #     )
        # thread_ls = []

        #--- Multithreading not yet implemented!
        #--- Look into the 'multiprocessing' module.
        ## Option 1: Start all threads at once. Join all at once.
        # for kb2d in self.KinBin2D_dict.values():
        #     thread_ls.append(threading.Thread(target=kb2d.do_itergausfit, args=arg_tup))
        #     thread_ls[-1].start()
        # for thread in thread_ls:
        #     thread_ls.join()

        ## Option 2: Start then join each thread.
        # thread_ls = []
        # for kb2d in self.KinBin2D_dict.values():
            # thread_ls.append(threading.Thread(target=kb2d.do_itergausfit, args=arg_tup))
            # thread_ls[-1].start()
            # thread_ls[-1].join()

        for kb2d in self.KinBin2D_dict.values():
            kb2d.do_itergausfit(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0,
                       x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0,
                       fit_whole_range_first_iter=fit_whole_range_first_iter,
                       iters=iters, num_sigmas=num_sigmas,
                       switch_to_binned_fit=switch_to_binned_fit,
                       verbose=verbose, alarm_level=alarm_level,
                       use_mu_pT_corr=use_mu_pT_corr, only_draw_last=only_draw_last)
            if verbose:
                # Show which KB2Ds are completed.
                kb2d_checkdct[kb2d.get_bin_key()] = "FITS COMPLETE"
                pprint(kb2d_checkdct)

    def do_3D_iter_gaus_fits(self, binned_fit=False, bins_dpTOverpT=100, bins_qd0=100, 
                             x_lim_dpTOverpT=[-0.4,0.4], x_lim_qd0=[-0.01,0.01],
                             fit_whole_range_first_iter=False,
                             iters=1, num_sigmas=2.5, switch_to_binned_fit=2000,
                             alarm_level="warning",
                             verbose=False, use_data_in_xlim=False):
        """
        Store a dict of iterated Gaussian fit statistics to each KinBin3D.
        Append KB3D plots to self.kinbin3d_iterfitplot_ls and self.kinbin3d_hist_ls.

        Plots:
        - Iterated Gaussian Fit of dpT/pT dist.
        - Binned dist of q*d0.
        
        NOTE: For each KinBin3D, performs a fit on the: 
            - dpT/pT distribution
            - qd0 distribution (Is this true?)
        """
        for kb2d in self.KinBin2D_dict.values():
            kb2d.is_qd0_binned = True
            for kb3d in kb2d.KinBin3D_dict.values():
                kb3d.analyze_KinBin3D(bins_dpTOverpT, bins_qd0, x_lim_dpTOverpT, x_lim_qd0,
                                      binned_fit=binned_fit, fit_whole_range_first_iter=fit_whole_range_first_iter,
                                      iters=iters, num_sigmas=num_sigmas,
                                      switch_to_binned_fit=switch_to_binned_fit, 
                                      verbose=verbose, alarm_level=alarm_level,
                                      use_data_in_xlim=use_data_in_xlim)
        self.kinbin3d_iterfitplot_ls.extend([kb3d.frame_1OverpT for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
        self.kinbin3d_iterfitplot_ls.extend([kb3d.frame_dpTOverpT for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
        self.kinbin3d_hist_ls.extend([kb3d.h_qd0 for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])

    def fill_KinBin_hists(self):
        """Fill the hists of each KinBin2D with its own muon kinematic values."""
        print("...Filling KinBin2D histograms...")
        for kinbin in self.KinBin2D_dict.values():
            kinbin.fill_hists()

    def make_KinBin2D_graphs(self, fit_with_zero_interc=False):
        """
        For each KinBin2D, use its stored muons to make all dpT/pT vs. qd0
        graphs:
        - dpT/pT            vs. q*d0
        - dpT/pT * 1/<pT>   vs. q*d0
        - dpT/pT * <1/pT>   vs. q*d0
        - dpT/pT * mu(1/pT) vs. q*d0
        Store the graph in self.kinbin2d_graph_ls.
        """
        print("...Making TGraphs for all KinBin2Ds...")
        for kb2d in self.KinBin2D_dict.values():
            kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_1divpT=False, fit_with_zero_interc=fit_with_zero_interc)
            kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_1divpT=True, fit_with_zero_interc=fit_with_zero_interc)
            kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_avgOf1divpT=True, fit_with_zero_interc=fit_with_zero_interc)
            # kb2d.make_dpTOverpT_graph(color=4, do_fit=True, scale_by_muOf1divpT=True, fit_with_zero_interc=fit_with_zero_interc)
            self.kinbin2d_graph_ls.extend([kb2d.gr_dpTOverpT_vs_qd0])
            self.kinbin2d_graph_ls.extend([kb2d.gr_dpTOverpTscaled_vs_qd0])
            self.kinbin2d_graph_ls.extend([kb2d.gr_dpTOverpTtimesavgOf1divpT])
            self.kinbin2d_graph_ls.extend([kb2d.gr_dpTOverpTtimesmuOf1divpT])

    def make_inclusive_kinematic_plots(self):
        """Store self.hist_inclusive_ls as a list of histograms.
        
        Run over self.muon_ls and self.m4mu_ls to make
        and fill all inclusive plots.
        """
        print(f"Creating histograms of inclusive muon kinematics...")
        title_prefix = "inclusive "
        
        label_m_inv = process_dct[self.prod_mode]["m_inv_label"]
        label_process = process_dct[self.prod_mode]["process"]
        m_inv_lim = process_dct[self.prod_mode]["m_inv_lim"]
        m_inv_nbins = process_dct[self.prod_mode]["m_inv_nbins"]

        h_invm = make_TH1F("h_invm", title=r"%s%s from %s" % (title_prefix, label_m_inv, label_process), n_bins=m_inv_nbins, xlabel=label_m_inv, x_min=m_inv_lim[0], x_max=m_inv_lim[1], units="GeV")
        h_pT = make_TH1F("h_pT", title=title_prefix+"p_{T,#mu}^{reco}", n_bins=220, xlabel=r"p_{T,#mu}^{reco}", x_min=0, x_max=220, units="GeV")
        h_pT_gen = make_TH1F("h_pT_gen", title=title_prefix+"p_{T,#mu}^{gen}", n_bins=220, xlabel=r"p_{T,#mu}^{gen}", x_min=0, x_max=220, units="GeV")
        h_eta = make_TH1F("h_eta", title=title_prefix+"#eta_{#mu}^{reco}", n_bins=100, xlabel=r"#eta_{#mu}^{reco}", x_min=-2.5, x_max=2.5, units=None)
        h_eta_gen = make_TH1F("h_eta_gen", title=title_prefix+"#eta_{#mu}^{gen}", n_bins=100, xlabel=r"#eta_{#mu}^{gen}", x_min=-2.5, x_max=2.5, units=None)
        h_phi = make_TH1F("h_phi", title=title_prefix+"#phi_{#mu}^{reco}", n_bins=50, xlabel=r"#phi_{#mu}^{reco}", x_min=-4, x_max=4, units=None)
        h_phi_gen = make_TH1F("h_phi_gen", title=title_prefix+"#phi_{#mu}^{gen}", n_bins=50, xlabel=r"#phi_{#mu}^{gen}", x_min=-4, x_max=4, units=None)
        h_qd0 = make_TH1F("h_qd0", title=title_prefix+"muon qd_{0}", n_bins=100, xlabel=r"qd_{0}", x_min=-0.01, x_max=0.01, units="cm")
        h_charge = make_TH1F("h_charge", title=title_prefix+"muon charge", n_bins=8, xlabel=r"Charge of muon", x_min=-4, x_max=4, units=None)
        h_charge.SetYTitle(r"Number of muons")
        h_dpTOverpT = make_TH1F("h_dpTOverpT", title=title_prefix+"(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", n_bins=100, xlabel=r"#Deltap_{T}^{reco}/p_{T}^{gen}", x_min=-0.3, x_max=0.3, units=None)

        if self.do_mu_pT_corr:
            h_pT_corr = make_TH1F("h_pT_corr", title=title_prefix+"p_{T,#mu}^{reco, corr.}", n_bins=220, xlabel=r"#Deltap_{T}^{reco}/p_{T}^{gen}", x_min=0, x_max=220, units=None)
            h_invm_corr = make_TH1F("h_invm_corr", title=title_prefix+"%s^{p_{T}, corr.} from %s" % (label_m_inv, label_process), n_bins=m_inv_nbins, xlabel=r"%s^{p_{T},corr.}" % label_process, x_min=105, x_max=140, units="GeV")
            h_invm_diff = make_TH1F("h_invm_diff", title=title_prefix+"#Delta%s #equiv %s^{corr. p_{T}} - %s" % (label_m_inv, label_m_inv, label_m_inv), n_bins=100, xlabel=r"#Delta%s" % label_m_inv, x_min=-10, x_max=10, units="GeV")
            h_invm_vs_invmcorr = make_TH2F("h_invm_vs_invmcorr", title=title_prefix+"Correlation between %s^{corr. p_{T}} and %s", 
                                n_binsx=100, x_label=r"%s" % label_m_inv, x_units="GeV", x_min=70, x_max=170,
                                n_binsy=100, y_label=r"%s^{p_{T},corr.}" % label_m_inv, y_units="GeV", y_min=70, y_max=170,
                                z_min=None, z_max=None, z_label_size=None,
                                n_contour=100)
            h_rel_dpT_corr2gen = make_TH1F("h_rel_dpT_corr2gen", title=title_prefix+"(p_{T}^{corr} - p_{T}^{gen})/p_{T}^{gen}", n_bins=100, xlabel=r"#Deltap_{T}^{corr}/p_{T}^{gen}", x_min=-0.3, x_max=0.3, units=None)
            h_rel_dpT_corr2rec = make_TH1F("h_rel_dpT_corr2rec", title=title_prefix+"(p_{T}^{corr} - p_{T}^{reco})/p_{T}^{reco}", n_bins=100, xlabel=r"#Deltap_{T}^{corr}/p_{T}^{reco}", x_min=-0.05, x_max=0.05, units=None)

            hist_inclusive_ls = [
                h_invm, h_invm_corr, h_invm_diff, h_invm_vs_invmcorr,
                h_dpTOverpT, h_rel_dpT_corr2gen, h_rel_dpT_corr2rec,
                h_pT_gen, h_pT, h_eta_gen, h_eta, h_phi_gen, h_phi,
                h_qd0, h_charge
            ]
        else:
            hist_inclusive_ls = [
                h_invm, h_dpTOverpT,
                h_pT_gen, h_pT, h_eta_gen, h_eta, h_phi_gen, h_phi,
                h_qd0, h_charge
            ]

        for h in hist_inclusive_ls:
            h.Sumw2()
        
        print(f"Filling histograms...")
        # Loop over stored values inside MyMuonCollection.
        # start = time.perf_counter()
        if self.do_mu_pT_corr:
            print(f"...Filling m_inv and m_inv_corr values...")
            # From here to below, "m4mu" is just m_inv (invariant mass of whatever process).
            for ct,(m4mu,m4mu_corr) in enumerate(zip(self.m4mu_ls, self.m4mu_corr_ls)):
                m4mu_diff = m4mu_corr - m4mu
                rel_diff = m4mu_diff / m4mu

                if abs(rel_diff) > 0.05:
                    msg = f"[INFO] Event {ct}: m_inv={m4mu}, m_inv_corr={m4mu_corr}, rel_diff={rel_diff}"
                    print(msg)

                h_invm.Fill(m4mu)
                h_invm_corr.Fill(m4mu_corr)
                h_invm_diff.Fill(m4mu_diff)
                h_invm_vs_invmcorr.Fill(m4mu, m4mu_corr)

            for mu in self.muon_ls:
                rel_dpT_corr2gen = (mu.pT_corr - mu.gen_pT) / mu.gen_pT
                rel_dpT_corr2rec = (mu.pT_corr - mu.pT) / mu.pT
                
                h_rel_dpT_corr2gen.Fill(rel_dpT_corr2gen)
                h_rel_dpT_corr2rec.Fill(rel_dpT_corr2rec)

                h_pT_corr.Fill(mu.pT_corr)

                h_pT.Fill(mu.pT)
                h_pT_gen.Fill(mu.gen_pT)
                h_eta.Fill(mu.eta)
                h_eta_gen.Fill(mu.gen_eta)
                h_phi.Fill(mu.phi)
                h_phi_gen.Fill(mu.gen_phi)
                h_qd0.Fill(mu.charge * mu.d0)
                h_charge.Fill(mu.charge)
                h_dpTOverpT.Fill(mu.dpTOverpT)
        else:
            if len(self.m4mu_ls) > 0:
                for m4mu in self.m4mu_ls:
                    h_invm.Fill(m4mu)
            for mu in self.muon_ls:
                h_pT.Fill(mu.pT)
                h_pT_gen.Fill(mu.gen_pT)
                h_eta.Fill(mu.eta)
                h_eta_gen.Fill(mu.gen_eta)
                h_phi.Fill(mu.phi)
                h_phi_gen.Fill(mu.gen_phi)
                h_qd0.Fill(mu.charge * mu.d0)
                h_charge.Fill(mu.charge)
                h_dpTOverpT.Fill(mu.dpTOverpT)
        # end = time.perf_counter()
        print(f"Filling histograms complete.")  # Took {end - start:.3f} seconds.")
        self.hist_inclusive_ls.extend(hist_inclusive_ls)

    def overwrite_kb2d_muon_info(self):
        """Overwrite each KinBin2D's long-listed attributes to save memory."""
        for kb2d in self.KinBin2D_dict.values():
            kb2d.overwrite_muon_info()

    def overwrite_kb3d_muon_info(self):
        """Overwrite each KinBin3D's long-listed attributes to save memory."""
        for kb3d in self.KinBin3D_dict.values():
            kb3d.overwrite_muon_info()
        
    def make_pT_corr_dict(self):
        """Store each KinBin2D's pT correction factors in a nested dict.

        NOTE: 
            - Each KinBin2D is associated with an (eta, pT) bin.
            - Each KinBin2D has its own set of pT corr factors (intecept, slope).
        
        Structure:
        {
            'dpTOverpT_vs_qd0' : {
                '0.0eta0.2_10.0pT14.0' : {
                    'intercept'     : -0.0004421981787913626,
                    'intercept_err' : 4.7910538089331156e-05,
                    'slope'         : 1.1863059538824643,
                    'slope_err'     : 0.02705983833855872
                },
                '0.0eta0.2_100.0pT150.0' : {
                    'intercept'     : ...,
                },
                ...,
            },
            'dpTOverpTscaled_vs_qd0' : {
                '0.0eta0.2_10.0pT14.0' : {
                    ...
                },
                ...,
            }
            'dpTOverpTtimesavgOf1divpT_vs_qd0' : {
                '0.0eta0.2_10.0pT14.0' : {
                    ...
                },
                ...,
            }
            'dpTOverpTtimesmuOf1divpT_vs_qd0' : {
                '0.0eta0.2_10.0pT14.0' : {
                    ...
                },
                ...,
            }
        }
        """
        self.pT_corr_factor_dict = {
            "dpTOverpT_vs_qd0"                 : {},
            "dpTOverpTscaled_vs_qd0"           : {},
            "dpTOverpTtimesavgOf1divpT_vs_qd0" : {},
            "dpTOverpTtimesmuOf1divpT_vs_qd0"  : {},
        }
        for kb2d in self.KinBin2D_dict.values():
            bin_key = kb2d.get_bin_key(title_friendly=False)
            for kb2d_graph_type in kb2d.dpTOverpT_vs_qd0_fit_stats_dct.keys():
                # Collect correction factors based on all correction methods
                # from KB2Ds.
                self.pT_corr_factor_dict[kb2d_graph_type][bin_key] = {
                    "intercept"     : kb2d.dpTOverpT_vs_qd0_fit_stats_dct[kb2d_graph_type]["interc_and_err"][0],
                    "intercept_err" : kb2d.dpTOverpT_vs_qd0_fit_stats_dct[kb2d_graph_type]["interc_and_err"][1],
                    "slope"         : kb2d.dpTOverpT_vs_qd0_fit_stats_dct[kb2d_graph_type]["slope_and_err"][0],
                    "slope_err"     : kb2d.dpTOverpT_vs_qd0_fit_stats_dct[kb2d_graph_type]["slope_and_err"][1],
                    }

    def get_all_kb2d_plots(self, *kind):
        """Return all plots associated with the KinBin2Ds in this MuonCollection.
        
        Parameters
        ----------
        kind : str, variable args
            "graphs", "hists", "dpT/pT hists"
            If `kind` is empty, then all plots will be returned.
        """
        d = {"graphs" : "gr_dpTOverpT_vs_qd0",
            "dpT/pT hists"  : "h_dpTOverpT",
            "qd0 hists"  : "h_qd0",}
        all_plots = []
        if len(kind) > 0:
            # User requests specific plots. Grab them.
            for k in kind:
                attr = d[k]
                plts = [getattr(kb2d, attr) for kb2d in self.KinBin2D_dict.values()]
                all_plots.extend(plts)
        else:
            # Get all available plots.
            for attr in d.values():
                plts = [getattr(kb2d, attr) for kb2d in self.KinBin2D_dict.values()]
                all_plots.extend(plts)
        return all_plots

    def get_all_kb3d_plots(self, *kind):
        """Return all plots associated with the KinBin3Ds in this MuonCollection.
        
        Parameters
        ----------
        kind : str, variable args
            "graphs", "hists", "dpT/pT hists"
            If `kind` is empty, then all plots will be returned.
        """
        d = {"1/pT iterfit" : "frame_1OverpT",
             "dpT/pT iterfit" : "frame_dpTOverpT",
            "dpT/pT hists"  : "h_dpTOverpT",
            "qd0 hists"  : "h_qd0",}
        all_plots = []
        if len(kind) > 0:
            # User requests specific plots. Grab them.
            for k in kind:
                attr = d[k]
                plts = [getattr(kb3d, attr) for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()]
                all_plots.extend(plts)
        else:
            # Get all available plots.
            for attr in d.values():
                plts = [getattr(kb3d, attr) for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()]
                all_plots.extend(plts)
        return all_plots

    def get_all_plots(self):
        """Return all plots associated with this MuonCollection."""
        # plot_ls.extend([hist for hist in self.hist_inclusive_ls])
        # plot_ls.extend([kb.h_qd0 for kb in self.KinBin2D_dict.values()])
        # plot_ls.extend([kb.h_dpTOverpT for kb in self.KinBin2D_dict.values()])
        # plot_ls.extend([kb.gr_dpTOverpT_vs_qd0 for kb in self.KinBin2D_dict.values()])
        all_plots = [
            *self.multigraph_ls,
            *self.multigraph_1divpT_ls,
            *self.hist_inclusive_ls,
            *self.get_all_kb2d_plots(),
            *self.get_all_kb3d_plots(),
        ]
        return all_plots

    def make_multigraph(self, eta_min, eta_max, y_lim=None,
                              scale_by_1divpT=False,
                              scale_by_avgOf1divpT=False,
                              scale_by_muOf1divpT=False):
        """
        TODO: Update docstring.
        Return a TMultiGraph of all KinBin2D TGraphs in a given eta bin.

        #--- Deprecated section below! Now returns mg ---#
        Append mg to either:
        - self.multigraph_ls, or
        - self.multigraph_1divpT_ls.append(mg).
        #--- Deprecated section above! Now returns mg ---#

        Parameters
        ----------
        eta_min : float
            Min eta value of all KB2Ds to be plotted.
        eta_max : float
            Max eta value of all KB2Ds to be plotted.
        y_lim : 2-elem list
            [y_min, y_max] plot window.
        scale_by_1divpT : bool
            If True, make dpT/pT * 1/<pT> vs. q*d0 plot for all KB2Ds.
        """
        assert sum((scale_by_1divpT, scale_by_avgOf1divpT, scale_by_muOf1divpT)) <= 1
        mg = ROOT.TMultiGraph(f"{eta_min}eta{eta_max}","")
        if y_lim is None:
            y_min = -0.04
            y_max = 0.16
        else:
            y_min = y_lim[0]
            y_max = y_lim[1]
        # if scale_by_1divpT or scale_by_avgOf1divpT or scale_by_muOf1divpT:
        #     # Zoom in.
        #     mg.SetMinimum(-0.002)
        #     mg.SetMaximum(0.013)
        # else:
        #     mg.SetMinimum(-0.04)
        #     mg.SetMaximum(0.16)
        mg.SetMinimum(y_min)
        mg.SetMaximum(y_max)
        ct = 0
        for kb2d in self.KinBin2D_dict.values():
            # Only add the KinBin2Ds within this eta bin.
            if kb2d.eta_range != [eta_min, eta_max]:
                continue
            ct += 1
            if scale_by_1divpT:
                gr = kb2d.gr_dpTOverpTscaled_vs_qd0
                gr.fit_line = kb2d.fit_line_scaled
            elif scale_by_avgOf1divpT:
                gr = kb2d.gr_dpTOverpTtimesavgOf1divpT
                gr.fit_line = kb2d.fit_line_scaled_avgOf1divpT
            elif scale_by_muOf1divpT:
                gr = kb2d.gr_dpTOverpTtimesmuOf1divpT
                gr.fit_line = kb2d.fit_line_scaled_muOf1divpT
            else:
                gr = kb2d.gr_dpTOverpT_vs_qd0
                gr.fit_line = kb2d.fit_line
            # # gr.GetXaxis().SetLimits(-0.008, 0.008)  # Doesn't work maybe because gr 
            gr.eta_min = eta_min
            gr.eta_max = eta_max
            gr.pT_min = kb2d.pT_min
            gr.pT_max = kb2d.pT_max
            if ct == 1:
                x_label = gr.GetXaxis().GetTitle()
                y_label = gr.GetYaxis().GetTitle()
                title  = r"#splitline{AdHoc p_{T} corrections derived from %s}" % process_dct[self.prod_mode]["process"]
                title += r"{%.2f < #left|#eta#right| < %.2f}" % (eta_min, eta_max)
                all_titles = r"%s;%s;%s" % (title, x_label, y_label)
                mg.SetTitle(all_titles)
            color = color_dict_RooFit[ct]
            gr.SetMarkerColor(color)
            gr.SetLineColor(color)
            gr.fit_line.SetLineColor(color)
            gr.fit_line.text_color = color
            # The `0` allows for error bars outside of range to still plot.
            # The `p` plots the points.
            mg.Add(gr, "0p")
        return mg
        # if scale_by_1divpT:
        #     self.multigraph_1divpT_ls.append(mg)
        # elif fit_line_scaled_avgOf1divpT:
        #     self.multigraph_avgOf1divpT_ls.append(mg)
        # elif fit_line_scaled_muOf1divpT:
        #     self.multigraph_muOf1divpT_ls.append(mg)
        # else:
        #     self.multigraph_ls.append(mg)

        # def make_all_multigraphs(self, eta_ls,
        #                                scale_by_1divpT=False,
        #                                scale_by_avgOf1divpT=False,
        #                                scale_by_muOf1divpT=False):
        #     """
        #     Make a multigraph (dpT/pT vs. q*d0) for each eta bin in eta_ls.
        #     If scale_by_1divpT, then make dpT/pT * 1/<pT> vs. q*d0 plots.
        #     """
        #     for eta_min, eta_max in zip(eta_ls[:-1], eta_ls[1:]):
        #         self.make_multigraph(eta_min, eta_max,
        #                             scale_by_1divpT=scale_by_1divpT,
        #                             scale_by_avgOf1divpT=scale_by_avgOf1divpT,
        #                             scale_by_muOf1divpT=scale_by_muOf1divpT)

    def draw_mg_and_fits(self, mg, x_lim=None,
                        scale_by_1divpT=False,
                        scale_by_avgOf1divpT=False,
                        scale_by_muOf1divpT=False,
                        draw_leg=True):
        """Draw one TMultiGraph and all corresponding fit lines.
        
        TODO: Finish docstring.

        Parameters
        ----------
        x_lim : 2-elem list
            [x_min, x_max] plot window.
        """
        gr_ls = list(mg.GetListOfGraphs())
        fitline_ls = [gr.fit_line for gr in gr_ls]
        assert len(gr_ls) == len(fitline_ls)

        # Add box for all fit stats.
        n_graphs = len(gr_ls)
        # Give each graph text on TPaveStats this height:
        gr_height = 0.03
        max_height = n_graphs * gr_height
        y_max = 0.88
        y_min = y_max - max_height
        pave_x_min = 0.17
        # Make TPave (legend) wider.
        pave_x_max = 0.72 if scale_by_1divpT or scale_by_avgOf1divpT or scale_by_muOf1divpT else 0.58
        leg = ROOT.TPaveText(pave_x_min, y_min, pave_x_max, y_max, "NDC")  # NDC = normalized coord.
        leg.SetFillColor(0)
        leg.SetFillStyle(1001)  # Solid fill = 1001.
        leg.SetBorderSize(1) # Use 0 for no border.
        leg.SetTextAlign(11) # 11 is against left side, 22 is centered vert and horiz.
        leg.SetTextSize(0.015)

        # Don't show stats box on multigraph.
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptFit(0)
        # Draw all graphs at once using TMultiGraph.
        mg.Draw("a")
        ROOT.gPad.Modified()  # Necessary to change multigraph's x-axis.
        if x_lim is None:
            x_min = -0.012
            x_max = 0.012
        else:
            x_min = x_lim[0]
            x_max = x_lim[1]
        mg.GetXaxis().SetLimits(x_min, x_max)
        # Draw fit lines.
        # Turn stats box back on.
        # ROOT.gStyle.SetOptStat("iouRMe")
        # ROOT.gStyle.SetOptFit(1)
        for fit, gr in zip(fitline_ls, gr_ls):
            interc = fit.GetParameter(0)
            slope = fit.GetParameter(1)
            eqn  = r"%.0f < p_{T} < %.0f GeV: " % (gr.pT_min, gr.pT_max)
            if scale_by_1divpT:
                eqn += r"#Deltap_{T}/p_{T}#upoint (1/<p_{T}>) = %.3E/GeV + (%.3f/[cm#upointGeV]) #upoint qd_{0}" % (interc, slope)
            elif scale_by_avgOf1divpT:
                eqn += r"#Deltap_{T}/p_{T}#upoint (<1/p_{T}>) = %.3E/GeV + (%.3f/[cm#upointGeV]) #upoint qd_{0}" % (interc, slope)
            elif scale_by_muOf1divpT:
                eqn += r"#Deltap_{T}/p_{T}#upoint (#mu_{Gauss}(1/p_{T})) = %.3E/GeV + (%.3f/[cm#upointGeV]) #upoint qd_{0}" % (interc, slope)
            else:
                eqn += r"#Deltap_{T}/p_{T} = %.3E + (%.3f/cm) #upoint qd_{0}" % (interc, slope)
            txt = leg.AddText(eqn)
            txt.SetTextColor(fit.text_color)
            fit.Draw("same")
        if draw_leg:
            leg.Draw("same")
        return leg
        # if scale_by_1divpT:
        #     self.multigraph_1divpT_leg_ls.append(leg)
        # elif scale_by_avgOf1divpT:
        #     self.multigraph_avgOf1divpT_leg_ls.append(leg)
        # elif scale_by_muOf1divpT:
        #     self.multigraph_muOf1divpT_leg_ls.append(leg)
        # else:
        #     self.multigraph_leg_ls.append(leg)

    def draw_all_multigraphs(self, outpath_pdf, printer,
                                scale_by_1divpT=False,
                                scale_by_avgOf1divpT=False,
                                scale_by_muOf1divpT=False):
        """Draw a list of mgs and their fit lines to a canvas, 1 mg per page.
        
        Parameters
        ----------
        outpath_pdf : str
            Absolute file path to pdf of plots.
        printer : CanvasPrinter
            Should have a canvas already attached.
        """
        if scale_by_1divpT:
            multigraph_ls = self.multigraph_1divpT_ls#.copy()
        elif scale_by_avgOf1divpT:
            multigraph_ls = self.multigraph_avgOf1divpT_ls#.copy()
        elif scale_by_muOf1divpT:
            multigraph_ls = self.multigraph_muOf1divpT_ls#.copy()
        else:
            multigraph_ls = self.multigraph_ls#.copy()
        assert len(multigraph_ls) > 0
        for mg in multigraph_ls:
            self.draw_mg_and_fits(mg, scale_by_1divpT=scale_by_1divpT,
                                      scale_by_avgOf1divpT=scale_by_avgOf1divpT,
                                      scale_by_muOf1divpT=scale_by_muOf1divpT)
            printer.canv.Print(outpath_pdf)
            ROOT.gPad.Update()
        printer.make_plots_pretty(show_statsbox=True)
        # Now draw each multigraph's line, individually.
        for mg in multigraph_ls:
            gr_ls = list(mg.GetListOfGraphs())
            for gr in gr_ls:
                # color = ct + 1 if ct >= 5 else ct
                # gr.fit_line.SetLineColor(color)
                gr.Draw("ap")
                gr.fit_line.Draw("same")  # gr.fit_line accounts for scaling.
                printer.canv.Print(outpath_pdf)

    def write_m4muinfo_to_rootfile(self, outpath_rootfile, overwrite=False):#, write_geofit_vals=False, ):
        """Write m4mu and m4mu_corr info to root file for DSCB fits, etc."""
        # New file, TTree, and TH1Fs.
        check_overwrite(outpath_rootfile, overwrite)
        print(f"  Writing m4mu and m4mu_corr vals to root file:")
        print(f"  {outpath_rootfile}")
        outf = ROOT.TFile(outpath_rootfile, "recreate")
        newtree = ROOT.TTree("tree", "tree_m4mu_vals")

        # Set pointers for values to be stored in root file.
        ptr_m4mu = array('f', [0.])
        ptr_m4mu_corr = array('f', [0.])
        newtree.Branch("m4mu", ptr_m4mu, "m4mu/F")
        newtree.Branch("m4mu_corr", ptr_m4mu_corr, "m4mu_corr/F")

        # Loop over stored values inside MyMuonCollection.
        print(f"...Filling m4mu and m4mu_corr values...")
        assert len(self.m4mu_ls) == len(self.m4mu_corr_ls) != 0
        for ct,(m4mu,m4mu_corr) in enumerate(zip(self.m4mu_ls, self.m4mu_corr_ls)):
            m4mu_diff = m4mu_corr - m4mu
            rel_diff = m4mu_diff / m4mu

            if abs(rel_diff) > 0.05:
                print(f"event {ct}:")
                print(f"  m4mu={m4mu}, m4mu_corr={m4mu_corr}, rel_diff={rel_diff}")

            ptr_m4mu[0] = m4mu
            ptr_m4mu_corr[0] = m4mu_corr
            newtree.Fill()
        print(f"  Number of good events found: {len(self.m4mu_ls)}")

        # Save tree and close file.
        newtree.Write()
        outf.Close()

    def save_KB2Ds_separate_dcts(self, outdir, file_prefix="", overwrite=False, verbose=False):
        """
        Saves each KB2D as its own standalone dict.
        Useful when a MyMuonCollection is too large to hold in RAM.

        Parameters
        ----------
        outdir : str
            The dir path where the KB2D dicts will be stored.
        file_prefix : str
            The prefix name of the outfile.
            Will be appended with '__{KB2D_bin_key}.pkl'.
        """
        make_dirs(outdir, verbose=verbose)
        if len(file_prefix) > 0:
            if file_prefix[-1] + file_prefix[-2] not in "__":
                file_prefix = f"{file_prefix}__" 
        for kb2d in self.KinBin2D_dict.values():
            key = kb2d.get_bin_key(title_friendly=True)
            outpath = os.path.join(outdir, f"{file_prefix}{key}.pkl")
            check_overwrite(outpath, overwrite=overwrite)
            if verbose:
                print(f"[INFO] Saving KB2D: {key}")
            save_to_pkl(kb2d, outpath, overwrite=overwrite)

    def apply_pTcorr_to_all_muons(self, pT_corr_factor_dict, use_GeoFit_algo=False,
                                  force_zero_intercept=True, verbose=False):
        """Apply pT corrections to all muons stored in the KB2Ds.
        
        TODO: Update docstring.
        
        NOTE: `self.KinBin2D_dict` must be populated.

        Parameters
        ----------
        """
        corr_type = "GeoFit" if use_GeoFit_algo else "AdHoc"
        print(f"...Applying {corr_type} pT corr factors to MyMuonCollection.")
        for kb2d in self.KinBin2D_dict.values():
            if verbose:
                print(f"...Working on: eta={kb2d.eta_range}, pT={kb2d.pT_range}")
            mu_ls = kb2d.muon_ls
            assert isinstance(mu_ls, list)
            assert len(mu_ls) > 0
            for mu in kb2d.muon_ls:
                mu.pT_corr = correct_muon_pT(
                    mu.eta, mu.pT, mu.charge, mu.d0,
                    pT_corr_factor_dict, detection="auto",
                    force_zero_intercept=force_zero_intercept,
                    use_GeoFit_algo=use_GeoFit_algo,
                    print_all_muon_info=False,
                    verbose=verbose)