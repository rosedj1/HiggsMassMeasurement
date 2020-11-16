import time
import pickle
import ROOT as r
from array import array
# Local imports.
from Utils_Python.Selections import build_muons_from_HZZ4mu_event
from Utils_Python.Utils_Physics import calc_Hmass
from Utils_Python.Utils_Files import check_overwrite
from d0_Studies.d0_Utils.d0_cls import KinBin2D
from d0_Studies.d0_Utils.d0_fns import (correct_muon_pT, parse_etapT_key,
                                       get_binedges_from_keys, make_key_from_binedges)
from d0_Utils.d0_fns import find_bin_edges_of_value

class MyMuonCollection:
    """A class to handle a list of MyMuon objects."""
    
    def __init__(self, prod_mode_ls=[]):
        self.prod_mode_ls = prod_mode_ls
        self.muon_ls = []
        self.m4mu_ls = []
        self.m4mu_corr_ls = []

        self.KinBin2D_dict = {}
        self.pT_corr_factor_dict = None
        self.do_mu_pT_corr = False   # Will turn True when pT corr is performed.

        self.plot_ls = []            # Specific KinBin2D plots.
        self.hist_inclusive_ls = []  # All muons in same plot.
        self.kinbin2d_plot_ls = []   # dpTOverpT vs. qd0 plots.
        self.kinbin3d_plot_ls = []   # KinBin3D distributions, like qd0 and dpTOverpT.

    def get_prod_modes_str(self):
        """Return a title string of all the production modes of the muons."""
        return ', '.join(self.prod_mode_ls)

    def extract_muons_from_Higgs_file(self, infile_path, n_evts=-1, print_out_every=10000,
                                      do_mu_pT_corr=False, 
                                      pT_corr_factor_dict=None,
                                      verbose=False):
        """Fill self.muon_ls with all MyMuon muons from infile_path
        which pass Higgs selections. Also fills self.m4mu_ls per event.
        
        If do_mu_pT_corr, then use pT_corr_factor_dict to
        correct the per muon pT.

        Parameters
        ----------
        infile_path : str
            Absolute file path to root file.
        n_evts : int
            Number of events to scan over.
            Default of `-1` runs over all events.
        print_out_every : int
            Show which event number is being processed in batches
            of `print_out_every`.
        do_mu_pT_corr : bool
            If True, perform muon pT correction based on ad hoc d0 studies.
            If no pT_corr_factor_dict is given, will default to
            self.pT_corr_factor_dict.
        pT_corr_factor_dict : dict, optional
            Example:
            {'1.25eta1.5_75.0pT100.0': {
                'intercept': 0.0035652354181878324,
                'intercept_err': 0.0017553367492442711,
                'slope': 8.787601605129735,
                'slope_err': 1.6856807294196279},
            '1.25eta1.5_100.0pT1000.0': {
                'intercept': 0.002248465529246451,
                'intercept_err': 0.002921444195364061,
                'slope': 14.523473508569932,
                'slope_err': 3.4317499063149564}
            }
        verbose : bool
            Print juicy debug info.
        """
        if do_mu_pT_corr:
            self.do_mu_pT_corr = True

        print(f"...Opening root file:\n{infile_path}")
        f = r.TFile(infile_path, "read")
        t = f.Get("Ana/passedEvents")

        all_evts = t.GetEntries()
        if n_evts == -1:
            n_evts = all_evts
        print(f"...Running over {n_evts} events...")

        # Event loop.
        for evt_num in range(n_evts):
            if evt_num % print_out_every == 0:
                print(f"  Running over evt: {evt_num}")
            # Check if this event passes Higgs selections.
            # Build the 4 MyMuons from the event.
            mu_tup = build_muons_from_HZZ4mu_event(t, evt_num)
            if None in mu_tup: 
                continue

            if do_mu_pT_corr:
                # Correct each muon's pT according to the
                # ad hoc pT correction factors.
                # Save the muons with old and new pTs.
                if pT_corr_factor_dict is None:
                    pT_corr_factor_dict = self.pT_corr_factor_dict
                assert pT_corr_factor_dict is not None
                corr_mu_ls = []
                for mu in mu_tup:
                    # For all 4 muons in this event, correct the muon pT.
                    # Then evaluate the m4mu_corr for this event.
                    mu.pT_corr = correct_muon_pT(
                        mu.eta, mu.pT, mu.charge, mu.d0, 
                        pT_corr_factor_dict, detection="auto",
                        verbose=verbose
                    )
                    lorentzvec_mu_corr = r.Math.PtEtaPhiMVector(mu.pT_corr, mu.eta, mu.phi, mu.mass)
                    corr_mu_ls.append(lorentzvec_mu_corr)
                assert len(corr_mu_ls) == 4
                m4mu_corr = calc_Hmass(*corr_mu_ls)
                self.m4mu_corr_ls.append(m4mu_corr)

            # Add the good muons to the final muon list.
            self.muon_ls.extend(mu_tup)
            # Save the m4mu info.
            lorentzvector_mu_ls = [mu.get_LorentzVector("reco") for mu in mu_tup]
            m4mu = calc_Hmass(*lorentzvector_mu_ls)
            self.m4mu_ls.append(m4mu)
        # End evt loop.
        assert len(self.muon_ls) > 0, "[ERROR] No muons were found!"
        assert len(self.muon_ls) / len(self.m4mu_ls) == 4
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
        self.make_KinBin_hists(n_bins_dpTOverpT=100, x_lim_dpTOverpT=[-0.4, 0.4],
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
        for eta_min,eta_max in zip(eta_ls[:-1], eta_ls[1:]):
            for pT_min,pT_max in zip(pT_ls[:-1], pT_ls[1:]):
                key = self.make_bin_key(eta_min, eta_max, pT_min, pT_max, title_friendly=False)
                this_KinBin = KinBin2D(eta_range=[eta_min,eta_max], pT_range=[pT_min,pT_max])
                self.KinBin2D_dict[key] = this_KinBin
        print(f"[INFO] All KinBin2Ds have been made from given eta_ls and pT_ls.")

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
            kb2d.store_muon_info_in_KinBin3Ds(title_friendly=False)

    def make_KinBin_hists(self, n_bins_dpTOverpT=100, x_lim_dpTOverpT=[-0.4, 0.4],
                                n_bins_qd0=100, x_lim_qd0=[-0.01, 0.01]):
        """Have each KinBin2D fill and store its own self.h_dpTOverpT and
        self.h_qd0 hists.
        
        NOTE:
        - Plots get stored like, e.g.: self.h_dpTOverpT, self.h_qd0
        - Used to be called: make_empty_KinBin_hists()
        """
        for kb in self.KinBin2D_dict.values():
            kb.make_dpTOverpT_hist(n_bins=n_bins_dpTOverpT, x_lim=x_lim_dpTOverpT)
            kb.make_qd0_hist(n_bins=n_bins_qd0, x_lim=x_lim_qd0)
        print(f"Done making KinBin hists")

    def place_muons_into_KinBins(self, eta_ls=None, pT_ls=None, 
                                 pT_corr_factor_dict=None, verbose=False):
        """Put muon into correct KinBin2D, based on muon's (eta, pT) value.
        
        After muons have been assigned to KB2D, self.muon_ls is cleared.
        """
        if verbose: 
            print(f"...Putting {len(self.muon_ls)} muons into their respective KinBin2Ds...")
        # Store muons in KinBins.
        if (eta_ls is None) or (pT_ls is None):
            # Automatically decide into which KinBin2D to put muons
            # based on the pT_corr_factor_dict keys.
            assert pT_corr_factor_dict is not None
            for mu in self.muon_ls:
                binedge_tup = get_binedges_from_keys(mu.eta, mu.pT, pT_corr_factor_dict)
                key = make_key_from_binedges(binedge_tup)
                self.KinBin2D_dict[key].add_muon(mu)
        else:
            n_muons_skipped = 0
            for muon in self.muon_ls:
                # Decide which (eta, pT) bin this muon belongs to.
                eta_min, eta_max = find_bin_edges_of_value(abs(muon.eta), eta_ls)
                pT_min, pT_max = find_bin_edges_of_value(muon.pT, pT_ls)
                if any([x is None for x in (eta_min, eta_max, pT_min, pT_max)]):
                    if verbose:
                        msg = f"[WARNING] Muon found with strange (eta, pT) values: eta={[eta_min, eta_max]}, pT={[pT_min, pT_max]}"
                        print(msg)
                        print("Skipping this muon, since your bins cannot hold it!")
                    n_muons_skipped += 1
                    continue
                # Choose correct KinBin2D to put this muon.
                key = self.make_bin_key(eta_min, eta_max, pT_min, pT_max, title_friendly=False)
                self.KinBin2D_dict[key].add_muon(muon)
            # End muon loop.
            if verbose:
                n_mu = len(self.muon_ls)
                print(f"Skipped {n_muons_skipped}/{n_mu} muons ({n_muons_skipped/n_mu * 100.})%")
        print(
            f"[INFO] All muons assigned to KinBin2Ds.\n"
            f"       Overwriting MyMuonCollection.muon_ls to save space."
            )
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
        for kb2d in self.KinBin2D_dict.values():
            kb2d.do_itergausfit(bins_dpTOverpT=bins_dpTOverpT, bins_qd0=bins_qd0,
                       x_lim_dpTOverpT=x_lim_dpTOverpT, x_lim_qd0=x_lim_qd0,
                       fit_whole_range_first_iter=fit_whole_range_first_iter,
                       iters=iters, num_sigmas=num_sigmas,
                       switch_to_binned_fit=switch_to_binned_fit,
                       verbose=verbose, alarm_level=alarm_level,
                       use_mu_pT_corr=use_mu_pT_corr, only_draw_last=only_draw_last)

    def do_3D_iter_gaus_fits(self, bins_dpTOverpT=100, bins_qd0=100, 
                             x_lim_dpTOverpT=[-0.4,0.4], x_lim_qd0=[-0.01,0.01],
                             fit_whole_range_first_iter=True,
                             iters=1, num_sigmas=2, switch_to_binned_fit=2000,
                             alarm_level="warning",
                             verbose=False):
        """Store a dict of iterated Gaussian fit statistics to each KinBin3D.
        
        NOTE: For each KinBin3D, performs a fit on the: 
            - dpT/pT distribution
            - qd0 distribution 
        """
        for kb2d in self.KinBin2D_dict.values():
            kb2d.is_qd0_binned = True
            for kb3d in kb2d.KinBin3D_dict.values():
                kb3d.analyze_KinBin3D(bins_dpTOverpT, bins_qd0, x_lim_dpTOverpT, x_lim_qd0,
                                      fit_whole_range_first_iter=fit_whole_range_first_iter,
                                      iters=iters, num_sigmas=num_sigmas,
                                      switch_to_binned_fit=switch_to_binned_fit, 
                                      verbose=verbose, alarm_level=alarm_level)

    def fill_KinBin_hists(self):
        """Fill the hists of each KinBin2D with its own muon kinematic values."""
        print("...Filling KinBin2D histograms...")
        for kinbin in self.KinBin2D_dict.values():
            kinbin.fill_hists()

    def make_KinBin2D_graphs(self):
        """For each KinBin2D, use its stored muons to make a dpT/pT vs. qd0 graph."""
        print("...Making TGraphs for all KinBin2Ds...")
        for kinbin in self.KinBin2D_dict.values():
            kinbin.make_dpTOverpT_graph(color=4, do_fit=True)

    def make_inclusive_kinematic_plots(self):
        """Store self.hist_inclusive_ls as a list of histograms.
        
        Run over self.muon_ls and self.m4mu_ls to make
        and fill all inclusive plots.
        """
        print(f"Creating histograms of inclusive muon kinematics...")
        title_prefix = "inclusive "
        
        h_m4mu = r.TH1F("h_m4mu", r"%sm_{4#mu} from %s" % (title_prefix, self.get_prod_modes_str()), 140, 105, 140)
        h_m4mu.SetXTitle(r"m_{4#mu} (GeV)")
        h_m4mu.SetYTitle(r"Events / (%.1f GeV)" % h_m4mu.GetBinWidth(1))

        h_pT = r.TH1F("h_pT", title_prefix+"p_{T,#mu}^{reco}", 220, 0, 220)
        h_pT.SetXTitle(r"p_{T,#mu}^{reco} (GeV)")
        h_pT.SetYTitle(r"Events / (%.1f GeV)" % h_pT.GetBinWidth(1))

        h_pT_gen = r.TH1F("h_pT_gen", title_prefix+"p_{T,#mu}^{gen}", 220, 0, 220)
        h_pT_gen.SetXTitle(r"p_{T,#mu}^{gen} (GeV)")
        h_pT_gen.SetYTitle(r"Events / (%.1f GeV)" % h_pT_gen.GetBinWidth(1))

        h_eta = r.TH1F("h_eta", title_prefix+"#eta_{#mu}^{reco}", 100, -2.5, 2.5)
        h_eta.SetXTitle(r"#eta_{#mu}^{reco}")
        h_eta.SetYTitle(r"Events / (%.1f)" % h_eta.GetBinWidth(1))

        h_eta_gen = r.TH1F("h_eta_gen", title_prefix+"#eta_{#mu}^{gen}", 100, -2.5, 2.5)
        h_eta_gen.SetXTitle(r"#eta_{#mu}^{gen}")
        h_eta_gen.SetYTitle(r"Events / (%.1f)" % h_eta_gen.GetBinWidth(1))

        h_phi = r.TH1F("h_phi", title_prefix+"#phi_{#mu}^{reco}", 50, -4, 4)
        h_phi.SetXTitle(r"#phi_{#mu}^{reco}")
        h_phi.SetYTitle(r"Events / (%.1f)" % h_phi.GetBinWidth(1))

        h_phi_gen = r.TH1F("h_phi_gen", title_prefix+"#phi_{#mu}^{gen}", 50, -4, 4)
        h_phi_gen.SetXTitle(r"#phi_{#mu}^{gen}")
        h_phi_gen.SetYTitle(r"Events / (%.1f)" % h_phi_gen.GetBinWidth(1))

        h_qd0 = r.TH1F("h_qd0", title_prefix+"muon qd_{0}", 100, -0.01, 0.01)
        h_qd0.SetXTitle(r"qd_{0} (cm)")
        h_qd0.SetYTitle(r"Events / (%.4f cm)" % h_qd0.GetBinWidth(1))

        h_charge = r.TH1F("h_charge", title_prefix+"muon charge", 100, -10, 10)
        h_charge.SetXTitle(r"Charge of muon")
        h_charge.SetYTitle(r"Number of muons")

        h_dpTOverpT = r.TH1F("h_dpTOverpT", title_prefix+"(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", 100, -0.3, 0.3)
        h_dpTOverpT.SetXTitle(r"#Deltap_{T}^{reco}/p_{T}^{gen}")
        h_dpTOverpT.SetYTitle(r"Events / (%.4f)" % h_dpTOverpT.GetBinWidth(1))

        if self.do_mu_pT_corr:
            h_pT_corr = r.TH1F("h_pT_corr", title_prefix+"p_{T,#mu}^{reco, corr.}", 220, 0, 220)

            h_m4mu_corr = r.TH1F("h_m4mu_corr", title_prefix+"m_{4#mu}^{p_{T}, corr.} from %s" % self.get_prod_modes_str(), 140, 105, 140)
            h_m4mu_corr.SetXTitle(r"m_{4#mu}^{p_{T},corr.} (GeV)")
            h_m4mu_corr.SetYTitle(r"Events / (%.1f GeV)" % h_m4mu_corr.GetBinWidth(1))

            h_m4mu_diff = r.TH1F("h_m4mu_diff", title_prefix+"#Deltam_{4#mu} #equiv m_{4#mu}^{corr. p_{T}} - m_{4#mu}", 100, -10, 10)
            h_m4mu_diff.SetXTitle(r"#Deltam_{4#mu} (GeV)")
            h_m4mu_diff.SetYTitle(r"Events / (%.1f GeV)" % h_m4mu_diff.GetBinWidth(1))

            h_m4muvsm4mucorr = r.TH2F("h_m4muvsm4mucorr", title_prefix+"Correlation between m_{4#mu}^{corr. p_{T}} and m_{4#mu}", 
                                        100, 70, 170, 100, 70, 170)  # n_binX, X_Low,X_Hig, n_binY, Y_low, Y_high
            h_m4muvsm4mucorr.SetXTitle(r"m_{4#mu} (GeV)")
            h_m4muvsm4mucorr.SetYTitle(r"m_{4#mu}^{p_{T},corr.} (GeV)")

            h_rel_dpT_corr2gen = r.TH1F("h_rel_dpT_corr2gen", title_prefix+"(p_{T}^{corr} - p_{T}^{gen})/p_{T}^{gen}", 100, -0.3, 0.3)
            h_rel_dpT_corr2gen.SetXTitle(r"#Deltap_{T}^{corr}/p_{T}^{gen}")
            h_rel_dpT_corr2gen.SetYTitle(r"Events / (%.4f)" % h_rel_dpT_corr2gen.GetBinWidth(1))

            h_rel_dpT_corr2rec = r.TH1F("h_rel_dpT_corr2rec", title_prefix+"(p_{T}^{corr} - p_{T}^{reco})/p_{T}^{reco}", 100, -0.05, 0.05)
            h_rel_dpT_corr2rec.SetXTitle(r"#Deltap_{T}^{corr}/p_{T}^{reco}")
            h_rel_dpT_corr2rec.SetYTitle(r"Events / (%.4f)" % h_rel_dpT_corr2rec.GetBinWidth(1))

            hist_inclusive_ls = [
                h_m4mu, h_m4mu_corr, h_m4mu_diff, h_m4muvsm4mucorr,
                h_dpTOverpT, h_rel_dpT_corr2gen, h_rel_dpT_corr2rec,
                h_pT_gen, h_pT, h_eta_gen, h_eta, h_phi_gen, h_phi,
                h_qd0, h_charge
            ]
        else:
            hist_inclusive_ls = [
                h_m4mu, h_dpTOverpT,
                h_pT_gen, h_pT, h_eta_gen, h_eta, h_phi_gen, h_phi,
                h_qd0, h_charge
            ]

        for h in hist_inclusive_ls:
            h.Sumw2()
        
        print(f"Filling histograms...")
        # Loop over stored values inside MyMuonCollection.
        start = time.perf_counter()
        if self.do_mu_pT_corr:
            print(f"...Filling m4mu and m4mu_corr values...")
            for ct,(m4mu,m4mu_corr) in enumerate(zip(self.m4mu_ls, self.m4mu_corr_ls)):
                m4mu_diff = m4mu_corr - m4mu
                rel_diff = m4mu_diff / m4mu

                if abs(rel_diff) > 0.05:
                    msg = f"[INFO] Event {ct}: m4mu={m4mu}, m4mu_corr={m4mu_corr}, rel_diff={rel_diff}"
                    print(msg)

                h_m4mu.Fill(m4mu)
                h_m4mu_corr.Fill(m4mu_corr)
                h_m4mu_diff.Fill(m4mu_diff)
                h_m4muvsm4mucorr.Fill(m4mu, m4mu_corr)

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
            for m4mu in self.m4mu_ls:
                h_m4mu.Fill(m4mu)
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
        end = time.perf_counter()
        print(f"Filling histograms complete. Took {end - start:.3f} seconds.")
        self.hist_inclusive_ls = hist_inclusive_ls

    def save_to_pkl(self, obj, outpkl_path):
        """Save one obj to pickle."""
        with open(outpkl_path, 'wb') as output:
            pickle.dump(obj, output, protocol=2)
        print(f"[INFO] Pickle file written:\n{outpkl_path}\n")

    def overwrite_longlist_muon_info(self):
        """Overwrite each KinBin2D's and KinBin3D's long-listed attributes to save memory."""
        for kb2d in self.KinBin2D_dict.values():
            kb2d.overwrite_muon_info()
            for kb3d in kb2d.KinBin3D_dict.values():
                kb3d.overwrite_muon_info()
        
    def make_pT_corr_dict(self):
        """Store each KinBin2D's pT correction factors in a nested dict.

        NOTE: 
            - Each KinBin2D is associated with an (eta, pT) bin.
            - Each KinBin2D has its own set of pT corr factors (intecept, slope).
        
        Structure:
            {
                "0.0eta0.2_5.0pT7.0" : {
                    "intercept"     : 7.899725827947741e-05,
                    "intercept_err" : blah,
                    "slope"         : 0.6683781742967103,
                    "slope_err"     : blah,
                    },
                "0.0eta0.2_7.0pT10.0"  : {
                    "intercept"     : -2.526319693758452e-05,
                    "intercept_err" : blah,
                    "slope"         : 0.8709094429939193,
                    "slope_err"     : blah,
                    },
                    ...
            }
        """
        for kb2d in self.KinBin2D_dict.values():
            key = kb2d.get_bin_key(title_friendly=False)
            self.pT_corr_factor_dict[key] = {
                "intercept"     : kb2d.interc_and_err[0],
                "intercept_err" : kb2d.interc_and_err[1],
                "slope"         : kb2d.slope_and_err[0],
                "slope_err"     : kb2d.slope_and_err[1],
                }

    # DELETE BELOW
    # def make_plots_beforeafterpTcorr(self):
    #     """Per KinBin2D, store a TCanvas showing before/after pT corr."""
    #     assert self.do_mu_pT_corr, "Has pT correction been performed?"
    #     for kb2d in self.KinBin2D_dict.values():
    #         kb2d.make_beforeafterpTcorr_canvas()
    # DELETE ABOVE

    def collect_all_plots(self):
        """Stores all plots associated with this MuonCollection."""
        self.plot_ls.extend([hist for hist in self.hist_inclusive_ls])
        self.plot_ls.extend([kb.h_qd0 for kb in self.KinBin2D_dict.values()])
        self.plot_ls.extend([kb.h_dpTOverpT for kb in self.KinBin2D_dict.values()])
        self.plot_ls.extend([kb.gr_dpTOverpT_vs_qd0 for kb in self.KinBin2D_dict.values()])

    def collect_KinBin2D_plots(self):
        """Stores all plots associated with the KinBin2D objects in this MuonCollection."""
        self.kinbin2d_plot_ls.extend([kb2d.gr_dpTOverpT_vs_qd0 for kb2d in self.KinBin2D_dict.values()])

    def collect_KinBin3D_plots(self):
        """Stores all plots associated with the KinBin3D objects in this MuonCollection."""
        # self.kinbin3d_plot_ls.extend([kb3d.h_dpTOverpT for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
        self.kinbin3d_plot_ls.extend([kb3d.frame_dpTOverpT for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
        self.kinbin3d_plot_ls.extend([kb3d.h_qd0 for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
        # self.kinbin3d_plot_ls.extend([kb3d.frame_qd0 for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])

    def write_m4muinfo_to_rootfile(self, outpath_rootfile, overwrite=False):
        """Write m4mu and m4mu_corr info to root file for DSCB fits, etc."""
        # New file, TTree, and TH1Fs.
        check_overwrite(outpath_rootfile, overwrite)
        outf = r.TFile(outpath_rootfile, "recreate")
        newtree = r.TTree("tree", "tree_m4mu_vals")

        # Set pointers for values to be stored in root file.
        ptr_m4mu = array('f', [0.])
        ptr_m4mu_corr = array('f', [0.])
        newtree.Branch("m4mu", ptr_m4mu, "m4mu/F")
        newtree.Branch("m4mu_corr", ptr_m4mu_corr, "m4mu_corr/F")

        # Loop over stored values inside MyMuonCollection.
        print(f"...Filling m4mu and m4mu_corr values...")
        start = time.perf_counter()
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
        end = time.perf_counter()
        print(f"(Took {end - start:.3f} seconds.)")

        # Save tree and close file.
        newtree.Write()
        outf.Close()