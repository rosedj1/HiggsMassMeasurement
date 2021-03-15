import time
import numpy as np
# import matplotlib.pyplot as plt
import ROOT as r

from scipy.optimize import curve_fit
from array import array
# Package imports.
# Not all of these may be used here. Just saving time for now.
from Utils_Python.Utils_Files import make_dirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Plotting import (change_cmap_bkg_to_white, save_plots_to_outpath, make_1D_dist, get_stats_1Dhist, 
                                        get_stats_2Dhist, hist_y_label, make_2by2_subplots_for_ratioplots,
                                        make_stats_legend_for_1dhist, make_stats_legend_for_2dhist, 
                                        make_stats_legend_for_gaus_fit)
from Utils_Python.Utils_Physics import theta2pseudorap, pseudorap2theta, calc_dR, calc_dphi, perc_diff
from Utils_Python.Utils_StatsAndFits import (linear_func, gaussian_func, 
                                            fit_with_gaussian, fit_with_line, 
                                            iterative_fit_gaus, check_fit_convergence,
                                            prop_err_x_div_y, prop_err_on_dsigoversig,
                                            get_bestfit_vals_from_statsdict)
from Utils_Python.Utils_Collection_Helpers import weave_lists
from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit
from Utils_ROOT.ROOT_Plotting import make_TPave
from d0_Utils.d0_fns import (make_binning_array, centers_of_binning_array, get_subset_mask, correct_muon_pT,
                             make_kinem_subplot, combine_cut_list, calc_x_err_bins_from_bin_centers, 
                             find_equal_hist_regions_unbinned, find_bin_edges_of_value)
from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict
        
class KinBin2D:
    """A container to organize muons in an (eta, pT) bin.
    
    NOTE:
        This class consumes lots of memory because it stores
        the muons as ROOT.Math.LorentzVector objects and then
        DUPLICATES them after splitting a KinBin2D up into 
        multiple KinBin3Ds.
        It would be more efficient to then delete the KinBin2D.muon_ls
        since the information is now stored among all the KinBin3Ds.
        - I have attempted to overwrite self.muon_ls,
        but not sure if that clears the memory.
        """

    def __init__(self, eta_range, pT_range):
        """
        Parameters
        ----------
        eta_range = 2-elem list: [eta_min, eta_max]
        pT_range = 2-elem list: [pT_min, pT_max]
        """
        self.kinbin_dim = 2
        self.eta_range = eta_range
        self.pT_range = pT_range

        self.eta_min = eta_range[0]
        self.eta_max = eta_range[1]
        self.pT_min = pT_range[0]
        self.pT_max = pT_range[1]

        self.h_qd0 = None
        self.h_dpTOverpT = None
        self.gr_dpTOverpT_vs_qd0 = None

        self.muon_ls = []
        self.KinBin3D_dict = {}

        self.iters = None
        self.is_qd0_binned = None
        self.fit_func = None
        self.interc_and_err = (None, None)
        self.slope_and_err = (None, None)
        self.equalentry_qd0_bin_edges = None
        self.use_mu_pT_corr = None

        self.fit_stats_dict_dpTOverpT = None
        self.frame_dpTOverpT = None

    def get_bin_key(self, title_friendly=False):
        """Return a string which identifies this 2D bin.
        
        Optional parameter: title_friendly 
            Useful for making file names and title names. 
            If True, replaces '.' with 'p'.
        """
        key = f"{self.eta_min}eta{self.eta_max}_{self.pT_min}pT{self.pT_max}"
        return key if not title_friendly else key.replace(".", "p")

    def make_3Dbin_key(self, qd0_min, qd0_max, title_friendly=False):
        """Use this KinBin2D's (eta, pT) info and the given qd0 info to make a str key."""
        key_tmp = self.get_bin_key(title_friendly=title_friendly)
        key = f"{key_tmp}_{qd0_min}qd0{qd0_max}"
        if (title_friendly):
            return key.replace(".", "p").replace("-", "neg")
        else:
            return key
            
    def get_bin_edges(self):
        """Return a 4-tuple of the eta and pT bin edges."""
        return (self.eta_min, self.eta_max, self.pT_min, self.pT_max)

    def make_latex_2Dbin_cut_str(self):
        """Return a LaTeX raw string for ROOT that shows the (eta, pT) values for this bin.
        
        This should just be for display purposes, never for using as an identifying key.
        """
        return r"%.2f < #left|#eta#right| < %.2f, %.0f < p_{T} < %.0f GeV" % self.get_bin_edges()
    
    def make_latex_bin_cut_str(self):
        """Return a label with this KinBin's cuts in LaTeX form."""
        if self.kinbin_dim == 2:
            return self.make_latex_2Dbin_cut_str()
        elif self.kinbin_dim == 3:
            return self.make_latex_3Dbin_cut_str()
        else:
            msg = f"self.kinbin_dim value ({self.kinbin_dim}) not understood."
            raise ValueError(msg)

    def make_qd0_hist(self, n_bins=100, x_lim=[-0.01, 0.01]):
        """Make and fill a qd0 hist using the muon info in this KinBin.
        
        x_lim [x_min, x_max] determines the x-axis range for viewing.

        NOTE: Can work with either KinBin2D or KinBin3D objects.
        """
        key = self.get_bin_key(title_friendly=True)
        latex_root_bin_str = self.make_latex_bin_cut_str()
        x_label = r"qd_{0} (cm)"
        title = r"%s" % (latex_root_bin_str)
        h = r.TH1F(f"h_qd0_{key}", title, n_bins, x_lim[0], x_lim[1])
        h.Sumw2()
        h.StatOverflows(True)
        h.SetXTitle(x_label)
        h.SetYTitle(r"Events / (%.4f)" % h.GetBinWidth(1))
        # Fill the histogram.
        for val in [mu.charge * mu.d0 for mu in self.muon_ls]:
            h.Fill(val)
        self.h_qd0 = h

    def make_dpTOverpT_hist(self, n_bins=100, x_lim=[-0.12, 0.12]):
        """Make and fill a delta_pT/pT hist using the muon info in this KinBin.
        
        x_lim [x_min, x_max] determines the x-axis range for viewing.

        NOTE: Can work with either KinBin2D or KinBin3D objects.
        """
        key = self.get_bin_key(title_friendly=True)
        latex_root_bin_str = self.make_latex_bin_cut_str()
        x_label = r"(p_{T}^{REC} - p_{T}^{GEN})/p_{T}^{GEN}"
        title = r"%s" % (latex_root_bin_str)
        h = r.TH1F(f"h_dpTOverpT_{key}", title, n_bins, x_lim[0], x_lim[1])
        h.Sumw2()
        h.StatOverflows(True)
        h.SetXTitle(x_label)
        h.SetYTitle(r"Events / (%.4f)" % h.GetBinWidth(1))
        # Fill the histogram.
        for val in [mu.dpTOverpT for mu in self.muon_ls]:
            h.Fill(val)
        self.h_dpTOverpT = h

    def make_empty_hists(self, bins_qd0=100, bins_dpTOverpT=100):
        """Create empty histograms which will be filled with muon kinematic info."""
        # Prepare plot labels.
        latex_root_bin_str = self.make_latex_bin_cut_str()
        key = self.get_bin_key(title_friendly=False)
        make_title = lambda x : r"#bf{muon %s distribution} (%s)" % (x, latex_root_bin_str)
        title_qd0 = make_title(r"qd_{0}")
        title_dpTOverpT = make_title(r"#Deltap_{T}/p_{T}")
        # Make histograms.
        h_qd0 = r.TH1F(f"h_qd0_{key}", f"{title_qd0}", bins_qd0, -0.01, 0.01)
        h_qd0.SetXTitle(r"qd_{0} (cm)")
        h_qd0.SetYTitle(r"Events / (%.4f)" % h_qd0.GetBinWidth(1))
        self.h_qd0 = h_qd0
        # label_LaTeX_dict["qd0BS1"]["independent_label_ROOT"]
        h_dpTOverpT = r.TH1F(f"h_dpTOverpT_{key}", title_dpTOverpT, bins_dpTOverpT, -0.6,  0.6)
        h_dpTOverpT.SetXTitle(r"#Deltap_{T}/p_{T}")
        h_dpTOverpT.SetYTitle(r"Events / (%.4f)" % h_dpTOverpT.GetBinWidth(1))
        self.h_dpTOverpT = h_dpTOverpT

    def get_qd0_ls(self):
        """Return a list of q*d0 values for all muons in this KB2D."""
        return [mu.charge * mu.d0 for mu in self.muon_ls]

    def get_dpTOverpT_ls(self):
        """Return a list of q*d0 values for all muons in this KB2D."""
        return [(mu.pT - mu.gen_pT)/mu.gen_pT for mu in self.muon_ls]

    def add_muon(self, muon):
        """Add this MyMuon and its info to this KinBin2D."""
        self.muon_ls.extend([muon])
    
    def fill_hists(self):
        """Take the muons in this KinBin and store the
        qd0 and dpT/pT kinematics per muon in corresponding hists."""
        for mu in self.muon_ls:
            self.h_qd0.Fill(mu.charge * mu.d0)
            self.h_dpTOverpT.Fill(mu.dpTOverpT)

    def make_dpTOverpT_graph(self, color=4, do_fit=True):
        """Store a TGraph of dpT/pT vs. qd0 using the stored muons in this KinBin.
        Also draw and store a best-fit line on the graph.
        """
        latex_cut_str = self.make_latex_bin_cut_str()
        print(f"...Building TGraph for bin: {latex_cut_str}")
        # Prepare x and y vals.
        x_vals, y_vals, x_err_vals, y_err_vals = [], [], [], []
        if (self.is_qd0_binned):
            y_label = r"#mu_{Gauss}(#Deltap_{T}/p_{T})"
            for kb3d in self.KinBin3D_dict.values():
                x_vals.append(kb3d.qd0_avg)
                y_vals.append(kb3d.fit_stats_dict_dpTOverpT["mean_ls"][-1])
                y_err_vals.append(kb3d.fit_stats_dict_dpTOverpT["mean_err_ls"][-1])
        else:
            x_vals = self.get_qd0_ls()
            y_vals = self.get_dpTOverpT_ls()
            y_err_vals = np.zeros_like(x_vals)
        x_err_vals = np.zeros_like(x_vals)
        # Convert them to array.arrays afterward.
        x_arr = array('f', x_vals)
        y_arr = array('f', y_vals)
        x_err_arr = array('f', x_err_vals)
        y_err_arr = array('f', y_err_vals)
        # Make the graph.
        n_pts = len(x_arr)
        assert n_pts == len(y_arr) == len(x_err_arr) == len(y_err_arr) != 0
        # gr = r.TGraph(n_pts, x_arr, y_arr)
        gr = r.TGraphErrors(n_pts, x_arr, y_arr, x_err_arr, y_err_arr)
        # Make it pretty.
        graph_title = r"(%s)" % (latex_cut_str)
        gr.SetMarkerStyle(21) #25
        gr.SetMarkerSize(0.2)
        #         gr.Draw("AP")
        gr.SetLineColor(color)
        gr.SetLineWidth(1)
        gr.SetMarkerColor(color)
        gr.SetTitle(graph_title)
        gr.GetXaxis().SetTitle(r"qd_{0} [cm]")
        gr.GetYaxis().SetTitle(y_label)
        gr.GetXaxis().SetLimits(-0.005, 0.005)
        gr.GetYaxis().SetRangeUser(-0.04, 0.04)
        self.gr_dpTOverpT_vs_qd0 = gr
        #         gr.GetYaxis().SetTitleOffset(1.5)
        # Must do 2 Draw() calls to set TGraphs properly...
        # if (draw):
        #     gr.Draw("AP")
        
        if (do_fit):
            self.fit_line = self.fit_graph_with_line(self.gr_dpTOverpT_vs_qd0)

    def fit_graph_with_line(self, gr):
        """Return a linear fit function and after fitting it to graph."""
        if (self.is_qd0_binned):
            # FIXME
            x_min = -0.006# Make this automatic: min([gr.GetPointX(2)])
            x_max = 0.006
        else:
            qd0_ls = self.get_qd0_ls()
            x_min = min(qd0_ls)
            x_max = max(qd0_ls)
        y_min = -0.02
        y_max = 0.02
        # Make the fit function.
        key = self.get_bin_key(title_friendly=False)
        fit_func = r.TF1(f"f1_{key}", '[0]+[1]*x', x_min, x_max)
        fit_func.SetLineColor(2)
        fit_func.SetLineWidth(1)
        fit_func.SetLineStyle(1)
        # Fit it onto a histogram `h1`:
        gr.Fit(fit_func,'S')
        # The option 'S' saves the fit results into a pointer.
        # r.gStyle.SetOptFit(111)
        fit_func.Draw("same")
        self.interc_and_err = (fit_func.GetParameter(0), fit_func.GetParError(0))
        self.slope_and_err = (fit_func.GetParameter(1), fit_func.GetParError(1))
        self.chi2 = fit_func.GetChisquare()
        self.NDF = fit_func.GetNDF()
        return fit_func

    def make_empty_equalentry_KinBin3Ds(self, regions, algo=("normal", -1), verbose=False, title_friendly=False):
        """Use the muons in this KinBin2D (eta, pT) to create a dict of equal-entry KinBin3Ds (eta, pT, qd0).

        NOTE: The qd0 bin edges are automatically determined in such a way that 
            there is an equal number of muons in each qd0 bin.
        
        regions : int
            The number of equal-entry regions to split the qd0 axis into.
            Each equal-entry region corresponds to a KinBin3D (eta_range, pT_range, qd0_range)

        Structure of KinBin3D dict:
            eta_ls = [0.2, 0.4]
            pT_ls  = [7.0, 10.0]  # NOTE: these lists do not need to be the same length!
            regions = 3  
            
            => Will make 3 KinBin3Ds each with an (eta, pT, qd0) bin: 
                {
                    "0.2eta0.4_7.0pT10.0_-0.026519qd0-0.0014461" : KinBin3D(0.2to0.4, 7.0to10.0, -0.026519to-0.0014461),
                    "0.2eta0.4_7.0pT10.0_-0.0014461qd0-0.0014522" : KinBin3D(0.2to0.4, 7.0to10.0, -0.0014461to-0.0014522),
                    "0.2eta0.4_7.0pT10.0_-0.0014522qd00.025578" : KinBin3D(0.2to0.4, 7.0to10.0, -0.026519to0.025578),
                }
        """
        self.equalentry_qd0_bin_edges, regions = find_equal_hist_regions_unbinned(self.get_qd0_ls(), regions, algo=algo, verbose=verbose)
        for qd0_min,qd0_max in zip(self.equalentry_qd0_bin_edges[:-1], self.equalentry_qd0_bin_edges[1:]):
            kb3d = KinBin3D(eta_range=self.eta_range, pT_range=self.pT_range, qd0_range=[qd0_min,qd0_max])
            key = self.make_3Dbin_key(qd0_min, qd0_max, title_friendly=title_friendly)
            self.KinBin3D_dict[key] = kb3d

    def store_muon_info_in_KinBin3Ds(self, title_friendly=False, verbose=False):
        """Save KinBin2D muon info in correct KinBin3Ds, based on (eta, pT, qd0) values of each muon."""
        for muon in self.muon_ls:
            # Decide which (eta, pT, qd0) bin this muon belongs to.
            qd0 = muon.charge * muon.d0
            # qd0_min, qd0_max = find_bin_edges_of_value(qd0, qd0_bin_edges_tmp)
            qd0_min, qd0_max = find_bin_edges_of_value(qd0, self.equalentry_qd0_bin_edges, verbose=verbose)
            if any([x is None for x in (qd0_min, qd0_max)]):
                msg = f"[WARNING] This muon has a strange qd0 value: qd0={[qd0_min, qd0_max]}"
                print(msg)
                print("Skipping this muon, since your bins cannot hold it!")
                continue
            # Put muon into correct KinBin3D.
            key = self.make_3Dbin_key(qd0_min, qd0_max, title_friendly=title_friendly)
            # self.KinBin3D_dict[key].store_muon_info(muon)
            self.KinBin3D_dict[key].add_muon(muon)

    def overwrite_muon_info(self, delete_all=True):
        """Overwrite long-listed attributes to save memory.
        
        NOTE: Not even sure if this actually clears the memory. Probably?
        - Look into Garbage Collection.
        """
        if delete_all:
            self.muon_ls = "overwritten"

    def do_itergausfit(self, data=None, bins_dpTOverpT=100, bins_qd0=100,
                       x_lim_dpTOverpT=[-0.4,0.4], x_lim_qd0=[-0.01,0.01],
                       fit_whole_range_first_iter=True,
                       iters=1, num_sigmas=2, marker_color=None, line_color=None,
                       switch_to_binned_fit=2000, verbose=False, alarm_level="warning",
                       use_mu_pT_corr=False, only_draw_last=False):
        """Perform an unbinned iterated Gaussian fit on the dpT/pT
        distribution. Store the statistics and best-fit vals of this KinBin.
        
        FIXME: Update Parameter descriptions.

        NOTE:
        - Should also work for child classes, like KinBin3D.
        - Requires: len(self.muon_ls) > 0.

        Parameters
        ----------
        data : list or array-like
            The data to be plotted and on which Gaussian fits will be done.
        bins_dpTOverpT : int
            The number of bins along the x-axis of the dpTOverpT distribution.
        bins_qd0 : int, DEPRECATED
            The number of bins along the x-axis of the qd0 distribution.
        x_lim_dpTOverpT : 2-elem list, optional
            The min and max x-axis values to show on the dpTOverpT plot: [x_min, x_max]
        x_lim_qd0 : 2-elem list, optional
            The min and max x-axis values to show on the qd0 plot: [x_min, x_max]
        fit_whole_range_first_iter : bool, optional
            If True, the entire x-axis range will be fit over for first iteration.
            Otherwise, first fit range is: mean(data) +- num_sigmas * rms(data)
        iters : int
            Number of Gaussian fit iterations to perform.
        num_sigmas : float
            FIXME: Finish this.
        marker_color : ROOT color or int
            Force the marker colors to be this.
            Default is kBlue+2 for before correction and kRed+2 for after correction.
        line_color : ROOT color or int
            Force the line color to be this.
            Default is blue for before correction and red for after correction.
            Useful when only showing the last iterated fit.
        switch_to_binned_fit : int, optional
            The max number of entries in an array on which an UNBINNED fit should be performed.
            If n_entries in array > switch_to_binned_fit, then a binned fit will be done.
        verbose : bool
            If True, print juicy debug info.
        use_mu_pT_corr : bool
            If True, use the corrected muon pT (from d0 studies).
            Make a distribution of (pT_rec_corr - pT_gen)/pT_gen.
        """
        self.n_entries = len(self.muon_ls)
        try:
            assert self.n_entries > 0
        except AssertionError:
            # Move along people, there's nothing to see here.
            print("There are no muons found in this KB2D. .")
            return
        if data is not None:
            try:
                assert len(data) > 0
            except AssertionError:
                print("There are no data to be plotted.")
                return
        # Make sure that the MyMuon has the requested attribute.
        msg = f"MyMuon doesn't have attribute `{attr}`"
        assert getattr(self.muon_ls[0], attr), msg

        self.iters = iters
        if verbose:
            start = time.perf_counter()
            print(
                f"...Analyzing KinBin{self.kinbin_dim}D: {self.get_bin_key()}...\n"
                f"...which has {self.n_entries} muons.\n"
                f"...Performing {iters} Gaussian fit iterations using RooFit...\n\n"
                )
        # self.qd0_avg = self.calc_avg_qd0()

        # Create self.h_qd0 and self.h_dpTOverpT.
        # make_empty_hists(self, bins_qd0=100, bins_dpTOverpT=100)
        # self.make_dpTOverpT_hist(n_bins=bins_dpTOverpT, x_lim=x_lim_dpTOverpT)
        # self.make_qd0_hist(n_bins=bins_qd0, x_lim=x_lim_qd0)
        marker_color_before_corr = r.kBlue+2 if marker_color is None else marker_color
        marker_color_after_corr = r.kRed+2 if marker_color is None else marker_color
        line_color_before_corr = 4 if line_color is None else line_color
        line_color_after_corr = 2 if line_color is None else line_color

        if use_mu_pT_corr:
            self.use_mu_pT_corr = use_mu_pT_corr
            # Perform iter Gauss fit on corrected reco pT vals.
            dpTOverpT_corr_ls = [(mu.pT_corr - mu.gen_pT)/mu.gen_pT for mu in self.muon_ls]

            self.fit_stats_dict_dpTOverpT_corr, self.frame_dpTOverpT_corr = RooFit_iterative_gaus_fit(
                                dpTOverpT_corr_ls, 
                                binned_fit=False, switch_to_binned_fit=switch_to_binned_fit, 
                                iters=iters, num_sigmas=num_sigmas,
                                n_bins=bins_dpTOverpT, x_lim=x_lim_dpTOverpT,
                                fit_whole_range_first_iter=fit_whole_range_first_iter,
                                xframe=None,
                                x_label=r"(p_{T}^{REC,corr.} - p_{T}^{GEN})/p_{T}^{GEN}", 
                                title="%s" % self.make_latex_bin_cut_str(), 
                                units="",
                                marker_color=marker_color_pT_corr, force_last_line_color=line_color_after_corr,
                                only_draw_last=only_draw_last, 
                                verbose=verbose)

        # Perform iter Gauss fit on reco pT vals.
        # dpTOverpT_ls = [mu.dpTOverpT for mu in self.muon_ls]
        # data = [getattr(mu, attr) for mu in self.muon_ls]

        self.fit_stats_dict_dpTOverpT, self.frame_dpTOverpT = RooFit_iterative_gaus_fit(
                                data, 
                                binned_fit=False, switch_to_binned_fit=switch_to_binned_fit, 
                                iters=iters, num_sigmas=num_sigmas,
                                n_bins=bins_dpTOverpT, x_lim=x_lim_dpTOverpT,
                                fit_whole_range_first_iter=fit_whole_range_first_iter,
                                xframe=None,
                                x_label=r"(p_{T}^{REC} - p_{T}^{GEN})/p_{T}^{GEN}", 
                                title="%s" % self.make_latex_bin_cut_str(), 
                                units="",
                                marker_color=marker_color_before_corr, force_last_line_color=line_color_before_corr,
                                only_draw_last=only_draw_last, 
                                verbose=verbose)

        # Print info on fit stats convergence.
        for stat in ["mean_ls", "mean_err_ls", "std_ls", "std_err_ls"]:
            if use_mu_pT_corr:
                check_fit_convergence(self.fit_stats_dict_dpTOverpT_corr[stat],
                                  max_perc_diff=5,
                                  compare_to_last=3,
                                  alarm_level=alarm_level)
            check_fit_convergence(self.fit_stats_dict_dpTOverpT[stat],
                                  max_perc_diff=5,
                                  compare_to_last=3,
                                  alarm_level=alarm_level)

        if verbose:
            end = time.perf_counter()
            print(
                f"Completed {iters} Iterated Gaussian fits "
                f"on KinBin{self.kinbin_dim}D: {self.get_bin_key()}\n"
                f"...took {(end - start):.6f} seconds.\n"
                )
            print("Now saving best-fit vals for this KinBin2D.")
        self.store_bestfit_vals()
        
    def store_bestfit_vals(self):
        """Store best-fit values from iterated Gauss fits of this KinBin2D."""
        bf_mean, bf_mean_err, bf_std, bf_std_err = get_bestfit_vals_from_statsdict(self.fit_stats_dict_dpTOverpT)
        self.bestfit_mean_beforecorr = bf_mean
        self.bestfit_meanerr_beforecorr = bf_mean_err
        self.bestfit_std_beforecorr = bf_std
        self.bestfit_stderr_beforecorr = bf_std_err

        if self.use_mu_pT_corr:
            bf_mean_corr, bf_mean_err_corr, bf_std_corr, bf_std_err_corr = get_bestfit_vals_from_statsdict(self.fit_stats_dict_dpTOverpT_corr)
            self.bestfit_mean_aftercorr = bf_mean_corr
            self.bestfit_meanerr_aftercorr = bf_mean_err_corr
            self.bestfit_std_aftercorr = bf_std_corr
            self.bestfit_stderr_aftercorr = bf_std_err_corr

            # Calculate percent improvement in sigma and its error.
            self.sigma_perc_improve = (bf_std_corr - bf_std) / bf_std * 100.
            self.sigma_perc_improve_err = prop_err_on_dsigoversig(
                                bf_std,
                                bf_std_corr,
                                bf_std_err,
                                bf_std_err_corr
                        ) * 100. # sig1, sig2, sig_err1, sig_err2

    def make_beforeafterpTcorr_frames(self):
        """Store the paves and frames for before/after pT correction."""
        frame_uncorr = self.frame_dpTOverpT.Clone()
        frame_corr   = self.frame_dpTOverpT_corr.Clone()

        # Move the hist stats boxes off of one another.
        statsbox_uncorr = frame_uncorr.findObject("stats")
        statsbox_uncorr.SetX1NDC(0.70)
        statsbox_uncorr.SetX2NDC(0.90)
        statsbox_uncorr.SetY1NDC(0.75)  # y min.
        statsbox_uncorr.SetY2NDC(0.90)
        statsbox_uncorr.SetTextColor(r.kBlue+2)

        statsbox_corr = frame_corr.findObject("stats")
        statsbox_corr.SetX1NDC(0.70)
        statsbox_corr.SetX2NDC(0.90)
        statsbox_corr.SetY1NDC(0.60)  # y min.
        statsbox_corr.SetY2NDC(0.75)
        statsbox_corr.SetTextColor(r.kRed+2)

        # # Move the fit stats boxes off of one another.
        fit_stats_uncorr = frame_uncorr.findObject("TPave")
        fit_stats_uncorr.SetX1NDC(0.13)
        fit_stats_uncorr.SetX2NDC(0.38)
        fit_stats_uncorr.SetY1NDC(0.77)  # y min.
        fit_stats_uncorr.SetY2NDC(0.88)
        fit_stats_uncorr.SetBorderSize(1)
        fit_stats_uncorr.SetTextColor(r.kBlue)
        line2change = fit_stats_uncorr.GetLineWith("Fit")
        oldtext = line2change.GetTitle().rstrip(":")
        line2change.SetTitle(r"%s before p_{T} corrections:" % oldtext)

        fit_stats_corr = frame_corr.findObject("TPave")
        fit_stats_corr.SetX1NDC(0.13)
        fit_stats_corr.SetX2NDC(0.38)
        fit_stats_corr.SetY1NDC(0.66)  # y min.
        fit_stats_corr.SetY2NDC(0.77)
        fit_stats_corr.SetBorderSize(1)
        fit_stats_corr.SetTextColor(r.kRed)
        # Change some of the text in this pave.
        line2change = fit_stats_corr.GetLineWith("Fit")
        oldtext = line2change.GetTitle().rstrip(":")
        line2change.SetTitle(r"%s after p_{T} corrections:" % oldtext)
        line2change = fit_stats_corr.GetLineWith("mu")
        newtext = line2change.GetTitle().replace("mu", r"mu_{corr.}")
        line2change.SetTitle(newtext)
        line2change = fit_stats_corr.GetLineWith("sigma")
        newtext = line2change.GetTitle().replace("sigma", r"sigma_{corr.}")
        line2change.SetTitle(newtext)

        # Make a box showing the improvement.
        pave_sig_improve = make_TPave(w=0.25, h=0.06, topright_corner_pos=(0.38, 0.66))
        txt_improve  = r"#frac{#left|#sigma_{corr.} - #sigma#right|}{#sigma} = "
        txt_improve += r"%.4g #pm %.4g" % (abs(self.sigma_perc_improve), self.sigma_perc_improve_err) + "%"
        pave_sig_improve.AddText(txt_improve)

        # Save everything in this KinBin2D.
        self.frame_dpTOverpT_clone = frame_uncorr
        self.frame_dpTOverpT_corr_clone = frame_corr
        self.statsbox_uncorr = statsbox_uncorr
        self.statsbox_corr = statsbox_corr
        self.fit_stats_uncorr = fit_stats_uncorr
        self.fit_stats_corr = fit_stats_corr
        self.pave_sig_improve = pave_sig_improve
        # Put two frames onto each other.
        # frame_uncorr.addObject(frame_corr)
        # frame_uncorr.Draw("same")
        # cname = self.get_bin_key(title_friendly=False)
        # c = r.TCanvas(cname, cname)
        # c.cd()
        # c.Draw()

        # c.Update()
        # self.frame_beforeafterpTcorr = frame_uncorr

    def draw_beforeafterpTcorr(self):
        """Draw a plot showing before/after pT corr for this KinBin2D."""
        self.frame_dpTOverpT_clone.Draw()
        self.frame_dpTOverpT_corr_clone.Draw("same")
        self.statsbox_uncorr.Draw("same")
        self.statsbox_corr.Draw("same")
        self.fit_stats_uncorr.Draw("same")
        self.fit_stats_corr.Draw("same")
        self.pave_sig_improve.Draw("same")

    def check_muon_belongs(self, mymuon):
        """Return True if mymuon has abs(eta) and pT that correspond to this KB2D."""
        return True if (self.eta_min < abs(mymuon.eta)) and (abs(mymuon.eta) < self.eta_max) else False

class KinBin3D(KinBin2D):
    """A container to organize muons in an (eta, pT, qd0) bin."""

    def __init__(self, eta_range=None, pT_range=None, qd0_range=None, 
        n_entries=0, kinem=None, fit_stats_dict=None, fit_type=None, 
        pT_stats_ls=None, qd0_stats_ls=None, cut_str=None):
        """
        fit_stats_dict : dict
            {
                "mean_ls"     : [mean_fit1, mean_fit2, ...]
                "mean_err_ls" : [mean_err_fit1, mean_err_fit2, ...]
                "std_ls"      : [sigma_fit1, sigma_fit2, ...]
                "std_err_ls"  : [sigma_err_fit1, sigma_err_fit2, ...]
            }
        fit_type : str
            Either "binned" or "unbinned"
        """
        self.kinbin_dim = 3
        self.eta_range = eta_range
        self.eta_min = eta_range[0]
        self.eta_max = eta_range[1]
        self.pT_range = pT_range
        self.pT_min = pT_range[0]
        self.pT_max = pT_range[1]
        self.qd0_range = qd0_range
        self.qd0_min = qd0_range[0]
        self.qd0_max = qd0_range[1]

        self.n_entries = n_entries
        self.qd0_avg = None

        self.muon_ls = []

        self.kinem = kinem
        self.fit_stats_dict = fit_stats_dict  # An older attribute from an older code.
        self.fit_stats_dict_dpTOverpT = None
        self.fit_stats_dict_qd0 = None
        self.fit_type = fit_type
        self.pT_stats_ls = pT_stats_ls
        self.qd0_stats_ls = qd0_stats_ls
        self.cut_str = cut_str

        self.h_dpTOverpT = None
    
    def get_bin_key(self, title_friendly=False):
        """Return a string which identifies this 3D bin.
        
        Optional parameter: title_friendly 
            Useful for making file names and title names. 
            If True, replaces '.' with 'p'.
        """
        key = (
            f"{self.eta_min}eta{self.eta_max}_"
            f"{self.pT_min:.1f}pT{self.pT_max:.1f}_"
            f"{self.qd0_min:.6f}qd0{self.qd0_max:.6f}"
            )
        return key if not title_friendly else key.replace(".", "p")

    def add_muon(self, muon):
        """Add this muon as a ROOT.Math.LorentzVector into this KinBin3D."""
        self.muon_ls.extend([muon])
        
    def make_latex_3Dbin_cut_str(self):
        """Return a LaTeX raw string for ROOT that shows the (eta, pT, qd0) values for this bin."""
        cut =  r"%.2f < #left|#eta#right| < %.2f, " % (self.eta_min, self.eta_max)
        cut += r"%.0f < p_{T} < %.0f GeV, " % (self.pT_min, self.pT_max)
        cut += r"%.3E < qd_{0} < %.3E cm" % (self.qd0_min, self.qd0_max)
        return cut
    
    def analyze_KinBin3D(self, bins_dpTOverpT, bins_qd0,
                         x_lim_dpTOverpT, x_lim_qd0,
                         binned_fit=False, fit_whole_range_first_iter=True,
                         iters=1, num_sigmas=2,
                         switch_to_binned_fit=2000, verbose=False, alarm_level="warning",
                         use_data_in_xlim=False):
        """Analyze and store dpT/pT and qd0 muon info within this KinBin3D.
        Also perform iterated Gaussian fits on the dpT/pT dist.
        
        Parameters
        ----------
        bins_dpTOverpT : int
            The number of bins along the x-axis of the dpTOverpT distribution.
        bins_qd0 : int
            The number of bins along the x-axis of the qd0 distribution.
        x_lim_dpTOverpT : 2-elem list, optional
            The min and max x-axis values to show on the dpTOverpT plot: [x_min, x_max]
        x_lim_qd0 : 2-elem list, optional
            The min and max x-axis values to show on the qd0 plot: [x_min, x_max]
        fit_whole_range_first_iter : bool, optional
            If True, the entire x-axis range will be fit over for first iteration.
            Otherwise, first fit range is: mean(data) +- num_sigmas * rms(data)
        iters : int
            Number of Gaussian fit iterations to perform.
        switch_to_binned_fit : int, optional
            The max number of entries in an array on which an UNBINNED fit should be performed.
            If n_entries in array > switch_to_binned_fit, then a binned fit will be done.
        verbose : bool
            If True, print juicy debug info.
        """
        self.n_entries = len(self.muon_ls)
        if verbose:
            print(
                f"[INFO] Analyzing KinBin3D:"
                f"[INFO]   {self.get_bin_key()} with {self.n_entries} entries\n"
                f"[INFO] Performing {iters} Gaussian fit iterations using RooFit"
                )
        self.qd0_avg = self.calc_avg_qd0()
        self.make_dpTOverpT_hist(n_bins=bins_dpTOverpT, x_lim=x_lim_dpTOverpT)
        self.make_qd0_hist(n_bins=bins_qd0, x_lim=x_lim_qd0)
        
        self.fit_stats_dict_dpTOverpT, self.frame_dpTOverpT = RooFit_iterative_gaus_fit(
                                self.get_dpTOverpT_ls(), 
                                binned_fit=binned_fit, switch_to_binned_fit=switch_to_binned_fit, 
                                iters=iters, num_sigmas=num_sigmas, 
                                n_bins=bins_dpTOverpT, x_lim=x_lim_dpTOverpT,
                                fit_whole_range_first_iter=fit_whole_range_first_iter,
                                xframe=None, x_label=r"(p_{T}^{REC} - p_{T}^{GEN})/p_{T}^{GEN}", 
                                title=r"%s" % self.make_latex_bin_cut_str(), 
                                units="", marker_color=1,
                                force_last_line_color=None, only_draw_last=False,
                                verbose=verbose, view_plot=False,
                                use_data_in_xlim=use_data_in_xlim)
        for stat in ["mean_ls", "mean_err_ls", "std_ls", "std_err_ls"]:
            check_fit_convergence(self.fit_stats_dict_dpTOverpT[stat],
                                  max_perc_diff=5,
                                  compare_to_last=3,
                                  alarm_level=alarm_level)

    def calc_avg_qd0(self):
        """Return the unweighted mean of the qd0 values of all muons in this KinBin3D."""
        return np.mean(self.get_qd0_ls())

    def correct_muon_ls_pT(self, pT_corr_factor_dct,
                           eta_binedge_ls, pT_binedge_ls,
                           verbose=False):
        """Return a list of the MyMuon objects in this KinBin3D
        with their pT corrected using the ad hoc d0 method."""
        corr_mu_ls = []
        for mu in self.muon_ls:
            mu.pT_corr = correct_muon_pT(mu.eta, mu.pT, mu.charge, mu.d0, 
                        pT_corr_factor_dct, 
                        eta_binedge_ls, pT_binedge_ls,
                        verbose=verbose)
            corr_mu_ls.append(mu)
        return corr_mu_ls

class OrganizerKB2D:
    """Organizes KinBin2D (eta, pT), for careful phase space analysis."""

    def __init__(self):
        """Create organizer to manage multiple KinBin2Ds."""
        pass

    def read_KB2D_from_pkl(self, inpath_pkl):
        """Return pickled dict of KinBin2Ds from inpath_pkl. Also store it."""
        self.kb2d_dct = open_pkl(inpath_pkl)
        # return self.kb2d_dct

class KinBin3DOrganizer():
    """
    ~~~ New and improved KinBinOrganizer ~~~
    Organizes KinBin3D objects by similar properties. Returns lists of KinBin objs.

    NOTE: This entire class should focus on the kinbin_obj being a dictionary.
        This provides MUCH better organization and much less 
        automatic figuring out how everything fits together.
    """
    def __init__(self, kinbin_obj):
        """Identify whether kinbin_obj is a dict or list and process it accordingly."""
        if isinstance(kinbin_obj, list):
            self.kinbin_ls = kinbin_obj
            self.identify_eta_pT_bins()
            self.make_KinBin_dict()
            
        elif isinstance(kinbin_obj, dict):
            self.kinbin_dict = kinbin_obj
            # FIXME: unpack dict?
            
        else:
            err = (
                f"  kinbin_obj type ({type(kinbin_obj)}) not recognized.\n"
                f"  Should be either type `list` or `dict`."
            )
            raise TypeError(err)
        
    def identify_eta_pT_bins(self):
        """
        Given the KinBin_obj that was given to this KinBin3DOrganizer, 
        automatically identify the unique eta ranges and pT ranges of all KinBin3D objs.

        Gives this KinBin3DOrganizer attributes like:
            self.eta_range_ls = [ [0.0,0.2], [0.2,0.4], ..., [2.3,2.4] ]
            self.pT_range_ls  = [ [5.0,7.0], [7.0,10.0], ..., [200.0, 1000.0] ]
            
            self.eta_ls     = [0.0, 0.2, 0.4, ..., 2.3, 2.4]
            self.pT_ls      = [5.0, 7.0, 10.0, ..., 200.0, 1000.0]

            self.eta_min_ls = [0.0, 0.2, 0.4, ..., 2.3]
            self.pT_min_ls  = [5.0, 7.0, 10.0, ..., 200.0]
        """
        # Identify available eta and pT values in kinbin_ls.
        eta_range_unique = set([tuple(kb.eta_range) for kb in self.kinbin_ls])  # Nested lists have to become tuples first.
        pT_range_unique = set([tuple(kb.pT_range) for kb in self.kinbin_ls])
        sorted_eta_range_unique = sorted(list(eta_range_unique))
        sorted_pT_range_unique = sorted(list(pT_range_unique))
        self.eta_range_ls = [list(tup) for tup in sorted_eta_range_unique]
        self.pT_range_ls = [list(tup) for tup in sorted_pT_range_unique]
        
        e_min = []
        e_max = []
        p_min = []
        p_max = []
        for kb in self.kinbin_ls:
            e_min.append(kb.eta_range[0])
            e_max.append(kb.eta_range[1])
            p_min.append(kb.pT_range[0])
            p_max.append(kb.pT_range[1])
        eta_min_unique = set(e_min)
        eta_max_unique = set(e_max)
        pT_min_unique = set(p_min)
        pT_max_unique = set(p_max)
        
        self.eta_ls = sorted(list(eta_min_unique | eta_max_unique))
        self.pT_ls = sorted(list(pT_min_unique | pT_max_unique))
        
        self.eta_min_ls = self.eta_ls[:-1]
        self.pT_min_ls  = self.pT_ls[:-1]
        
        print(f"Number of KinBins: {len(self.kinbin_ls)}")
        print(f"           eta_ls: {self.eta_ls}")
        print(f"            pT_ls: {self.pT_ls}")
        
    def make_KinBin_dict(self):
        """
        Creates a dictionary of KinBin lists. 
        
        Takes in HUGE list of a KinBin3D objects, a list of eta bin edges, and a 
        list of pT bin edges, and sorts KinBin3D objs based on which ones 
        share the same ranges. 
        eta, pT are keys, KinBin3D lists are the values.
        """
        # Make the dict.
        kinbin_dict = {}
        for eta_min in self.eta_min_ls:
            kinbin_dict[eta_min] = {}

            for pT_min in self.pT_min_ls:
                kinbin_dict[eta_min][pT_min] = [kb for kb in self.kinbin_ls if (kb.eta_range[0] == eta_min) and (kb.pT_range[0] == pT_min)]
        
        self.kinbin_dict = kinbin_dict
    
    def find_similar_KinBins(self, eta_range, pT_range):
        """
        Takes in list of a KinBin3D objects and, depending on what eta and pT range you give it,
        will return a subset of the original list which falls within those eta and pT ranges.
        """
        kb_same_eta = [kb for kb in self.kinbin_ls if kb.eta_range == eta_range]
        kb_same_eta_same_pT = [kb for kb in kb_same_eta if kb.pT_range == pT_range]
        
        return kb_same_eta_same_pT
    
    def find_all_KinBin_ls_const_range(self, const_range=[-1,-1], const_bin="pT"):

        ls_2D = []
        if const_bin in "eta":
            for pT_range in self.pT_range_ls:
                ls_2D.append( self.find_similar_KinBins(const_range, pT_range) )
                
        elif const_bin in "pT":
            for eta_range in self.eta_range_ls:
                ls_2D.append( self.find_similar_KinBins(eta_range, const_range) )
                
        return ls_2D
    
    def get_plotting_vals_from_KinBin_ls(self, kinbin_ls):
        
        # Make sure each KinBin3D is in the same pT and eta range.
        # Works for a SINGLE GRAPH LINE.
        # Must first pass in a list of KinBin3D:
        assert len(set([kb.pT_range for kb in kinbin_ls][0])) == 2
        assert len(set([kb.eta_range for kb in kinbin_ls][0])) == 2

        qd0_avg_ls = [kb.qd0_stats_ls[1] for kb in kinbin_ls]
        dpToverpT_bestfitmean_ls = [kb.fit_stats_dict["mean_ls"][-1] for kb in kinbin_ls]
        dpToverpT_bestfitmean_err_ls = [kb.fit_stats_dict["mean_err_ls"][-1] for kb in kinbin_ls]

        return qd0_avg_ls, dpToverpT_bestfitmean_ls, dpToverpT_bestfitmean_err_ls

class GraphLineKinBin3D:
    """
    One of the lines drawn on a graph. Contains all the info that went into building this line. 

    NOTE: This entire class should focus on the kinbin_obj being a dictionary
        This provides MUCH better organization and much less 
        automatic figuring out how everything fits together.
    """
    import numpy as np
    # Eventually just give this a list of KinBin3D obj.
    def __init__(self, x_vals, y_vals, x_err_vals=None, y_err_vals=None, eta_range=None, pT_range=None):
        self.x_vals = x_vals
        self.y_vals = y_vals
        self.x_err_vals = x_err_vals
        self.y_err_vals = y_err_vals
        self.eta_range = eta_range
        self.pT_range = pT_range
        
    def draw_graph(self, x_label="", y_label="", title="", kbin_ls=None, scale_by_1divpT=False, const_bin="pT", 
                   ax=None, x_lim=None, y_lim=None, count=1, verbose=False, 
                   fit_line=True, x_fit_range=None, legend_str_yequals="", legend_str_xvar=""):
        """
        Draws data points (values of: kinem_x, kinem_y) to an axes object. 
        In particular, used for making dpT/pT vs. q*d0 plots, but could probably be generalized.
        
        Parameters
        ----------
        x_label : str
            The x-axis label. If no x_label is given, then an automatic one 
            is generated based on kinem_x.
        y_label : str
            The y-axis label. If no y_label is given, then an automatic one 
            is generated based on kinem_y.
        title : str
        
        kbin_ls : list
            Should contain KinBin3D objects all within the same eta, pT bin. 

        scale_by_1divpT : bool
            If True, scale point k on this line by avg(pT_k),
            where pT_k is the average pT of cube k. 
        const_bin : 


        binning_type : str
            Must be either 'eta' or 'pT'. Used for proper labeling of title and legend.
        kbin_example : KinematicBin object
            This KinematicBin contains all the cut information necessary for proper
            legend and axes labeling.
        ax : axes object
            The axes on which to draw the graph. 
            If an axes is not provided, a default one is made.
        count : int
            A key to a dictionary of colors. 
            Values of the dict are color strings, like: 'black', 'red', etc. 
        """
        # FIXME: Each kb in kbin_ls has the same eta and pT range.
        # So is this necessary? 
        if const_bin in "pT":
            # Should be the case that each KinBin has same pT range.
            tup_ls = [tuple(kb.pT_range) for kb in kbin_ls]
        elif const_bin in "eta":
            # Similar argument for eta range.
            tup_ls = [tuple(kb.eta_range) for kb in kbin_ls]
        assert len(set(tup_ls)) == 1
        
        kbin_example = kbin_ls[0]
        eta_range = kbin_example.eta_range
        pT_range = kbin_example.pT_range

        if ax is None:
            f, ax = plt.subplots(figsize=(12.8, 9.6))
            
        # al=1  # alpha=0 is transparent
        # elw=1  # error bar line width
        # ms=1.5  # marker size
        # fontsize_legend = 6
        
        ecolor=color_dict[count]
        mec=color_dict[count]  # Marker edge color.
        mfc=color_dict[count]  # Marker face color.
        # cs=1.5  # cap size
        # mew=0.7  # marker edge width
        markerstyle = "None"  # For ax.errobar() it needs to be "None".

        if len(x_label) == 0:
            # Need x_label_base for legend.
            x_label_base = r"avg $q(\mu^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{BS} }$ [cm]"
            x_label = x_label_base + "\n" + r"in $(\left| \eta \right|, p_{T}, q*d_{0})$ cube"

        if len(y_label) == 0:
            y_label = r"iter. Gaus fit $\mu( \Delta p_{T} \ / p_{T}^{\mathrm{REC}})$"
        
        eta_text = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$" % (eta_range[0], eta_range[1])
        pT_text = r"$%.1f < p_{T}^{\mathrm{REC}} < %.1f$ [GeV]" % (pT_range[0], pT_range[1])
        
        if (const_bin in "pT"):
            label_text = eta_text
            title = pT_text
        elif (const_bin in "eta"):
            label_text = pT_text
            title = eta_text
        else:
            pass
        
        if (scale_by_1divpT):
            avg_pT_for3Dcubes = np.asarray([kb.pT_stats_ls[1] for kb in kbin_ls])#   float(list(kbin_example.pT_stats_ls[0])[1])
            avg_pT_err_for3Dcubes = np.asarray([kb.pT_stats_ls[2] for kb in kbin_ls])#  float(list(kbin_example.pT_stats_ls[0])[2])
            
            pT_scale_text = r" $* \frac{1}{\mathrm{avg}(p_{T}^{\mathrm{REC}})}$"
            y_label += pT_scale_text
            
            # Propagate errors on Gaus(mu) / avg(pT)
            gaus_mu_vals = np.asarray(self.y_vals)
            gaus_mu_err_vals = np.asarray(self.y_err_vals)
            
            y_vals, y_err_vals = prop_err_x_div_y(gaus_mu_vals, avg_pT_for3Dcubes, gaus_mu_err_vals, avg_pT_err_for3Dcubes)
            if (verbose):
                print(
                    f"About to propagate error...\n"
                    f"  gaus_mu_vals: {gaus_mu_vals}\n"
                    f"  avg_pT_for3Dcube: {avg_pT_for3Dcubes}\n"
                    f"  gaus_mu_err_vals: {gaus_mu_err_vals}\n"
                    f"  avg_pT_err_for3Dcube: {avg_pT_err_for3Dcubes}\n"
                    f"  y_vals: {y_vals}\n"
                    f"  y_err_vals: {y_err_vals}\n"
                )
        else:
            y_vals = self.y_vals
            y_err_vals = self.y_err_vals
            
        x_vals = self.x_vals
        
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        if x_lim is not None:
            ax.set_xlim(x_lim)
        if y_lim is not None:
            ax.set_ylim(y_lim)
        
        if (fit_line):
            ax.errorbar(x_vals, y_vals, xerr=self.x_err_vals, yerr=y_err_vals, fmt=markerstyle,
                    color=color_dict[count], 
                        # elinewidth=elw, ms=ms, 
                        # mec=mec, 
                        # capsize=cs, mew=mew,
                        mfc=mfc, ecolor=ecolor)
            ax = self.do_linear_fit(x_vals, y_vals, x_fit_range=x_fit_range, ax=ax, count=count, 
                                    leg_label_text=label_text, legend_str_yequals=r"$\Delta p_T/p_T$", legend_str_xvar=x_label_base, 
                                    scale_by_1divpT=scale_by_1divpT)
        else:
            ax.errorbar(x_vals, y_vals, xerr=self.x_err_vals, yerr=y_err_vals, fmt=markerstyle,
                        color=color_dict[count], 
                            # elinewidth=elw, ms=ms, 
                            # mec=mec, 
                            # capsize=cs, mew=mew,
                            mfc=mfc, ecolor=ecolor, label=label_text)  # This has label text here.
            ax.legend(loc="upper left", framealpha=al)#, fontsize=fontsize_legend)
        
        # Don't show d0 cuts and the cuts of whatever binning type (like "eta") is being used.
    #         tmp_dict = kbin_example.cut_dict.copy()
    #         for key in list(kbin_example.cut_dict.keys()):
    #             if (binning_type in key) or ("d0" in key):
    #                 del tmp_dict[key]
                        
    #         sorted_cut_ls = [value for (key, value) in sorted(tmp_dict.items())]
    #         cut_str = combine_cut_list(sorted_cut_ls)
    #         textbox_text = "Selections:\n" + cut_str  # Don't show the d0 cut text. Luckily it is the first by alphabetical sorting. 
            
    #         if count == 1:
    #             ax.text(0.05, 0.85, textbox_text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)

    def do_linear_fit(self, x_vals, y_vals, x_fit_range, ax=None, count=1, leg_label_text="", legend_str_yequals="", legend_str_xvar="", scale_by_1divpT=False):
        """
        NOTE: This method finds the best-fit line and saves the parameters:
            self.best_fit_params = [best_fit_intercept, best_fit_slope]
        """
        if ax is None:
            f, ax = plt.subplots()

        popt, popt_err, pcov = fit_with_line(x_vals, y_vals)
        best_fit_intercept = popt[0]
        best_fit_slope = popt[1]
        
        self.best_fit_params = [best_fit_intercept, best_fit_slope]

        if x_fit_range is None:
            # Get it automatically from the x_vals.
            x_fit_range = [min(x_vals), max(x_vals)]        
        x_fit_range = np.linspace(x_fit_range[0], x_fit_range[1], 50)
        
        y_fit_vals = linear_func(x_fit_range, best_fit_intercept, best_fit_slope)
        if (scale_by_1divpT):
            label_text = r", %s * 1/avg$(p_{T})$ = %.3E/[GeV] + %.3E/([GeV*cm]) * (%s) " % (legend_str_yequals, best_fit_intercept, best_fit_slope, legend_str_xvar)
        else:
            label_text = r", %s = %.3E + %.3E/[cm] * (%s) " % (legend_str_yequals, best_fit_intercept, best_fit_slope, legend_str_xvar)
            
        final_leg_text = leg_label_text + label_text

        ax.plot(x_fit_range, y_fit_vals, color=color_dict[count], marker="", label=final_leg_text, alpha=0.7)
        
        ax.legend(loc="upper left", framealpha=1)#, fontsize=text_size_legend)
        return ax

class GraphLine():
    """
    OLD VERSION!

    One of the lines drawn on a graph. Contains all the info that went into building this line. 
    """
    def __init__(self, x_vals, y_vals, y_err_vals=np.zeros(0)):
        self.x_vals = x_vals
        self.y_vals = y_vals
    #         self.x_err_vals = x_err_vals
        self.y_err_vals = y_err_vals
        
    def draw_graph(self, kinem_x, kinem_y, x_label="", y_label="", binning_type="", kbin_example=None, ax=None, count=1):
        """
        Draws data points (values of: kinem_x, kinem_y) to an axes object. 
        In particular, used for making dpT/pT vs. d0q plots, but could probably be generalized.
        
        Parameters
        ----------
        kinem_x : str
            The full name of the kinematic variable plotted on the x-axis.
            Should be a key in the label_LaTeX_dict.
        kinem_y : str
            The full name of the kinematic variable plotted on the y-axis.
            Should be a key in the label_LaTeX_dict.
        x_label : str
            The x-axis label. If no x_label is given, then an automatic one 
            is generated based on kinem_x.
        y_label : str
            The y-axis label. If no y_label is given, then an automatic one 
            is generated based on kinem_y.
        binning_type : str
            Must be either 'eta' or 'pT'. Used for proper labeling of title and legend.
        kbin_example : KinematicBin object
            This KinematicBin contains all the cut information necessary for proper
            legend and axes labeling.
        ax : axes object
            The axes on which to draw the graph. 
            If an axes is not provided, a default one is made.
        count : int
            A key to a dictionary of colors. 
            Values of the dict are color strings, like: 'black', 'red', etc. 
        """
        if binning_type not in ["eta", "pT"]:
            raise ValueError("[ERROR] Wrong `binning_type` specified. Must be either 'pT' or 'eta'. Stopping now.")
            
        if ax is None:
            f, ax = plt.subplots(figsize=(12.8, 9.6))
            
    #         #--- Check that things make sense: ---#
    #         # Example: If binning in eta, make sure each HistInfo object has identical pT_cuts as every other.
    #         wrong_eta_binning = (binning_type in 'eta') and len(set([(hist.pT_range[0], hist.pT_range[1]) for hist in entire_HistInfo_list])) != 1
    #         # Do same thing for binning in pT.
    #         wrong_pT_binning = (binning_type in 'pT') and len(set([(hist.eta_range[0], hist.eta_range[1]) for hist in entire_HistInfo_list])) != 1
    #         if (wrong_eta_binning or wrong_pT_binning):
    #             err_msg = f"\n\nBinning type ({binning_type}) specified, "
    #             err_msg += f"but not all graphs share same {binning_type}_range. Stopping now."
    #             raise RuntimeError(err_msg)

        al=1  # alpha=0 is transparent
        elw=1  # error bar line width
        ms=4  # marker size
        ecolor=color_dict[count]
        mec=color_dict[count]  # Marker edge color.
        mfc=color_dict[count]  # Marker face color.
        cs=1  # cap size
        mew=0.7  # marker edge width

        if len(x_label) == 0:
            x_label = label_LaTeX_dict[kinem_x]["independent_label"]
            unit_x = label_LaTeX_dict[kinem_x]["units"]
            if len(unit_x) > 0:
                x_label += " [{}]".format(unit_x)
        if len(y_label) == 0:
            y_label  = label_LaTeX_dict[kinem_y]["independent_label"]
            y_label += " (iterated Gaus fit mean)"
        title = label_LaTeX_dict[binning_type + "1"]["independent_label"] + " Binning"
    #         if binning_type == "eta":
    #             title = r"$\left| $" + title + r"$\right| $"
        ax.set_title(title)
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        
        # The "x-errors" are calculated automatically to be 1/2 the distance to the next data point. 
        low_x_err, high_x_err = calc_x_err_bins_from_bin_centers(self.x_vals)
        
        label_text = kbin_example.cut_dict[binning_type]
        if (kbin_example.verbose):
            print("Drawing graph, binning in {}:".format(binning_type))
            print(kbin_example.cut_dict[binning_type] + "\n")
        ax.errorbar(self.x_vals, self.y_vals, xerr=[low_x_err, high_x_err], yerr=self.y_err_vals, fmt='s', label=label_text,
                #color=color_dict[count], 
                    elinewidth=elw, ms=ms, mec=mec, capsize=cs, mew=mew, mfc=mfc, ecolor=ecolor)
        ax.legend(loc="lower right", framealpha=al)#, fontsize=text_size_legend)
        
        # Don't show d0 cuts and the cuts of whatever binning type (like "eta") is being used.
        tmp_dict = kbin_example.cut_dict.copy()
        for key in list(kbin_example.cut_dict.keys()):
            if (binning_type in key) or ("d0" in key):
                del tmp_dict[key]
                    
        sorted_cut_ls = [value for (key, value) in sorted(tmp_dict.items())]
        cut_str = combine_cut_list(sorted_cut_ls)
        textbox_text = "Selections:\n" + cut_str  # Don't show the d0 cut text. Luckily it is the first by alphabetical sorting. 
        
        if count == 1:
            ax.text(0.05, 0.85, textbox_text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        
    def do_linear_fit(self, ax=None):
        if ax is None:
            f, ax = plt.subplots(figsize=(12.8,9.6)) 
        
        # Do fit. 
        # Draw fit on axes.
        # return optimized parameters
        self.popt_linear = None

class KinematicBinVaex():
    # FIXME: I attempted to make it work with vaex. It almost works!
    def __init__(self, 
                vaex_df, 
    #                  inpath_dataframe, 
                 n_evts, 
                 massZ_cut_ls, eta_cut_ls, pT_cut_ls, qd0_cut_ls, d0_type="BS", dR_cut=0.02, 
                 use_ptotal_instead=False, verbose=False):
        """
        Pass in a DataFrame (DF) and specify the eta and pT cuts to create a subset of DF.
        
        Parameters
        ----------
        vaex_df : pandas.DataFrame
            ROOT file converted into a DataFrame. Columns are branches. Rows are events.
        n_evts : int
            Number of events to search over - not guaranteed to find this many events! 
            Use '-1' to loop over all events in the df. 
        eta_cut_ls : list or array-like of floats
            A list of [eta_min, eta_max]. Example: [0.9, 1.8]
        pT_cut_ls : list or array-like of floats
            A list of [pT_min, pT_max]. Example: [5, 20]
        qd0_cut_ls : list or array-like of floats
            A list of values of d0*charge to cut on: [qd0_min, qd0_max]. E.g. [-0.01, 0.01]
        d0_type : str
            Which d0 to cut on: "BS" or "PV"
        use_ptotal_instead : bool
            Cut on total momentum instead of pT. 
            In most places, pT could stand for either p or pT (depending on 'use_ptotal_instead').
                Just be careful because not all places are adapted for p!
        dR_cut : float
            A cut to save events in which muon1 and muon2 both have dR < dR_cut.
        verbose : bool
            If True, get debug info and see where you are while the code runs.
            
        NOTE:
            The methods further down are more developed than the methods closer to __init__(). 
            Therefore, clean up and consolidate the earlier methods. 
        """
        if n_evts == -1:
            n_evts = len(vaex_df)
        df = vaex_df[:n_evts]    # Original DF. 
        print("length of df:",len(df))
        self.n_evts_asked_for = n_evts
        self.n_evts_found = -999
        self.kinem_vals_after_selection = {}
        self.stats_dict = {}
        self.verbose = verbose
        
        self.cuts = ""
        self.sorted_cut_ls = []
        self.cut_dict = {}
        self.use_ptotal_instead = use_ptotal_instead
        self.p_str = "p" if (self.use_ptotal_instead) else "pT"
        self.p_str_latex = r"$p^{\mathrm{REC}}$" if (self.use_ptotal_instead) else r"$p_{T}^{\mathrm{REC}}$"

        self.massZ_min = massZ_cut_ls[0]
        self.massZ_max = massZ_cut_ls[1]        
        self.eta_min   = eta_cut_ls[0]
        self.eta_max   = eta_cut_ls[1]
        self.pT_min    = pT_cut_ls[0] 
        self.pT_max    = pT_cut_ls[1]         
        self.qd0_min   = qd0_cut_ls[0]
        self.qd0_max   = qd0_cut_ls[1]
        self.d0_type   = d0_type
        self.dR_cut    = dR_cut   
                
        self.apply_initial_cuts(df, verbose)
        
    def apply_initial_cuts(self, df, verbose):
        """
        Creates a subset of the original DataFrame in which initial cuts are applied.
        Cuts:
            pT
            eta
            q*d0
            massZ
            dR
        """
    # VDF is good up to here at least.
        # Cuts:
        
        
        print("Show1:", df.column_names)
        # Create masks.
        print("Getting mask massZ...") 
        self.mask_massZ = mask_massZ = self.get_mask_massZ(df)
        print("Getting mask eta...")
        self.mask_eta1, self.mask_eta2 = mask_eta1, mask_eta2 = self.get_mask_eta(df)
        print("Getting mask pT...")
        self.mask_pT1,  self.mask_pT2  = mask_pT1,  mask_pT2  = self.get_mask_pT(df)
        print("Getting mask qd0...")
        self.mask_qd01, self.mask_qd02 = mask_qd01, mask_qd02 = self.get_mask_qd0(df)    
        print("Getting mask dR...")
        self.mask_dR1,  self.mask_dR2  = mask_dR1,  mask_dR2  = self.get_mask_dR(df) 
        
        # Combine masks.
        print("Individual masks...")
        self.mask_kinembin_lep1 = mask_massZ & mask_dR1 & mask_eta1 & mask_pT1 & mask_qd01
        self.mask_kinembin_lep2 = mask_massZ & mask_dR2 & mask_eta2 & mask_pT2 & mask_qd02

        # Keep all events in which either muon1 passed all selections or muon2 passed all. 
        print("Combining masks...")
        self.all_masks = self.mask_kinembin_lep1 | self.mask_kinembin_lep2
        print("Applying all masks...")
        # Apply masks and update DataFrame.
    #         self.binned_df = df[self.all_masks]
        df_cuts_applied = df.select(self.all_masks)
        print("vdf attached to object")
        self.n_evts_found = len(df_cuts_applied)
        print("empty fucker")
        # The cut_dict has been filled. Now convert it to an ordered list (alphabetically).
        self.sorted_cut_ls = sorted_cut_ls = [value for (key, value) in sorted(self.cut_dict.items())]
        self.cuts = combine_cut_list(sorted_cut_ls)
        
        if (self.verbose): 
            perc = self.n_evts_found / float(self.n_evts_asked_for) * 100.
            print("[INFO] Events found: {} ({:.3f}% of total events scanned)".format(self.n_evts_found, perc))
            print(r"using cuts: {}".format(self.cuts) + "\n")
      
    def get_mask_qd0(self, df):
        d0_type = self.d0_type
        if d0_type == "PV":
            mask_qd01 = (self.qd0_min < df['qd0PV1']) & (df['qd0PV1'] < self.qd0_max)
            mask_qd02 = (self.qd0_min < df['qd0PV2']) & (df['qd0PV2'] < self.qd0_max)
        elif d0_type == "BS":
            mask_qd01 = (self.qd0_min < df['qd0BS1']) & (df['qd0BS1'] < self.qd0_max)
            mask_qd02 = (self.qd0_min < df['qd0BS2']) & (df['qd0BS2'] < self.qd0_max)
            
        cuts_qd0 = r"$%.3f < q(\mu)*d_{0}^{\mathrm{%s}} < %.3f$" % (self.qd0_min, self.d0_type, self.qd0_max)
        key = "qd0{}".format(self.d0_type)
        self.cut_dict[key] = cuts_qd0
        print("mask_qd01:", mask_qd01.sum())
        print("mask_qd02:", mask_qd02.sum())
        return mask_qd01, mask_qd02
 
    def get_mask_pT(self, df):
        if (self.use_ptotal_instead):
            mask_pT1 = (self.pT_min < df['p1']) & (df['p1'] < self.pT_max) 
            mask_pT2 = (self.pT_min < df['p2']) & (df['p2'] < self.pT_max)
        else:
            mask_pT1 = (self.pT_min < df['pT1']) & (df['pT1'] < self.pT_max) 
            mask_pT2 = (self.pT_min < df['pT2']) & (df['pT2'] < self.pT_max)
        
        cuts_p = r"$%d <$ %s $< %d$ GeV" % (self.pT_min, self.p_str_latex, self.pT_max)  # The string brings in its own '$'.
        self.cut_dict[self.p_str] = cuts_p
        print("mask_pT1:", mask_pT1.sum())
        print("mask_pT2:", mask_pT2.sum())
        return mask_pT1, mask_pT2

    def get_mask_eta(self, df):
        mask_eta1 = (self.eta_min < abs(df['eta1'])) & (abs(df['eta1']) < self.eta_max)
        mask_eta2 = (self.eta_min < abs(df['eta2'])) & (abs(df['eta2']) < self.eta_max)   
        
        cuts_eta = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$" % (self.eta_min, self.eta_max)
        self.cut_dict["eta"] = cuts_eta
        print("mask_eta1:", mask_eta1.sum())
        print("mask_eta2:", mask_eta2.sum())
        return mask_eta1, mask_eta2

    def get_mask_dR(self, df):
        print("df['delta_R1']:", df['delta_R1'][:30])
        print("df['delta_R2']:", df['delta_R2'][:30])
        print("self.dR_cut:", self.dR_cut)
        mask_dR1 = (df['delta_R1'] < self.dR_cut)
        mask_dR2 = (df['delta_R2'] < self.dR_cut)
        
        cuts_dR = r"$\Delta R < %.3f$" % (self.dR_cut)
        self.cut_dict["delta_R"] = cuts_dR
        print("mask_dR1:", mask_dR1.sum())
        print("mask_dR2:", mask_dR2.sum())
        return mask_dR1, mask_dR2
    
    def get_mask_massZ(self, df):
        mask_massZ = (self.massZ_min < df['massZ']) & (df['massZ'] < self.massZ_max)
        
        cuts_massZ = r"$%.1f < m_{\mu\mu} < %.1f$ GeV" % (self.massZ_min, self.massZ_max)
        self.cut_dict["massZ"] = cuts_massZ
        print("mask_massZ:", mask_massZ.sum())
        return mask_massZ    
    
    def apply_mask_get_data(self, df, kinem, lep_selection_type="", weave=False):
        """
        Return the kinematic values of the certain leptons by grouping them in different ways. 
            
            lep_selection_type = "1" -- get kinem values in which muon1 passes all selection criteria 
                                       (muon 2 may or may not pass selections).
            lep_selection_type = "2" -- get kinem values in which muon2 passes all selection criteria.
            lep_selection_type = "both" -- BOTH muons must pass selections to get data from event.
            lep_selection_type = "either" -- Either muon1, or muon2, or both must pass selections to get data from event.
            lep_selection_type = "independent" -- Get kinematic values for muon1 and muon2, with no restriction on which event they came from.
                Note: If kinem ends in "1" or "2", then the kinematic values of the other lepton are automatically grabbed.
                
        E.g. Apply a boolean mask for event selection and retrieve all "delta_R1" values.

        Parameters
        ----------
        df : vaex DataFrame
            DF with cuts already applied (massZ, pT, etc.).
        kinem : str
            A complete branch name in the DataFrame or root file. 
                E.g. "pT1", "genLep_pt2", "massZ"
        lep_selection_type : int
            The lepton's mask you want to apply.
        weave : bool
            Weave lep1 kinematic values and lep2 kinematic values together, so that slicing doesn't just give lep1. 
            Only relevant for "independent" selection.
        
        Returns
        -------
        kinem_vals : array 
            Kinematic values with chosen mask applied. 
            All values satisfy selection criteria for this kinematic bin.
        """
        mask1 = self.mask_kinembin_lep1
        mask2 = self.mask_kinembin_lep2
        
        if lep_selection_type == "1":
            # Only select events in which muon1 passes selections.
            mask = mask1
        elif lep_selection_type == "2":
            # Only select events in which muon2 passes selections.
            mask = mask2
        elif lep_selection_type == "both":
            # Only select events in which BOTH muon1 and muon2 pass selections.
            mask = mask1 & mask2
        elif lep_selection_type == "either":
            # Only select events in which BOTH muon1 and muon2 pass selections.
            mask = mask1 | mask2
        elif lep_selection_type == "independent":
            # Go through all muons in all events, without regard for other muon. 
            if kinem[-1] in ["1", "2"]:
                # Lep1 (or lep2) kinematic detected. Go find the other lepton's kinematic values.
                # FIXME: the variable massZ_vtxChi2 will be wrongly caught by this 'if' statement!
                kinem1 = kinem[:-1] + "1"
                kinem2 = kinem[:-1] + "2"
                kinem_vals1 = df[kinem1][mask1].values
                kinem_vals2 = df[kinem2][mask2].values
            else:
                # The kinematic doesn't depend on lep1 or lep2, like: massZ, GENmass2l, etc.
                kinem_vals1 = df[kinem][mask1].values
                kinem_vals2 = df[kinem][mask2].values
            
            if (weave):
                # Weave values together so that when slicing (like [:5]), you don't just grab kinem_vals1. 
                kinem_vals = np.array( weave_lists(kinem_vals1, kinem_vals2) )
            else: 
                kinem_vals = np.append(kinem_vals1, kinem_vals2)
        
            return kinem_vals
        
        else: 
            raise ValueError("[ERROR] `lep_selection_type` was not specified properly. Stopping now.")
        
        # A selection, other than "independent" was chosen.
        kinem_vals = df[kinem][mask].values
        
        return kinem_vals
    
    def make_2D_plot(self, 
                     df, 
                     x_kinem, y_kinem, 
                     x_bin_limits=[0, 1, 0.1], y_bin_limits=[0, 1, 0.1],
                     lep_selection_type="",
                     run_over_only_n_evts=-1, 
                     title="",
                     exclusive=True,
                     save_plot=False, save_as_png=False, outpath="",
                     ax=None):
        """
        Make a 2D plot. Two examples:
            (1) dphi vs. deta  
            (2) dphi vs. dtheta
        User can specify the binning along either axis. 
        This method plots only the muons which pass the selection 
        (as opposed to taking any event in which at least 1 muon pass kinematic bin criteria).
        
        Parameters
        ----------
        x_kinem : str
            The PARTIAL name of the kinematical variable to be plotted along x-axis. 
            Only works for kinematics which end with '1' or '2'.
            - Example: x_kinem="delta_theta" (for which there are two branches: "delta_theta1", "delta_theta2")
        y_kinem : str
            The PARTIAL name of the kinematical variable to be plotted along y-axis. 
            Only works for kinematics which end with '1' or '2'.
            - Example: y_kinem="delta_eta" (for which there are two branches: "delta_eta1", "delta_eta2")
        x_bin_limits : list or array-like of floats
            The bin limits on the horizontal axis. [bin_min_left_edge, bin_max_right_edge, bin_width]
            - Example: [-2.5, 2.5, 0.1]
        y_bin_limits : list or array-like of floats
            The bin limits on the vertical axis. [bin_min_left_edge, bin_max_right_edge, bin_width]
            - Example: [-2.5, 2.5, 0.1]
        lep_selection_type : str
            What kind of selection to perform on the leptons. Choices:
            #UPDATE
        run_over_only_n_evts : int
            Number of events to plot. Use '-1' to use all events in this kinembin.
        title : str
            Alternate title to put on plot. Overrides the default one made in this method.
        exclusive : bool
            Means "only put muons which passed all selections in this plot".
            FIXME: It must be set to True for now...
        save_plot : bool
            If True, save the plot as a pdf and possibly a png.
        save_as_png : bool
            If True, save the plot as a png.
        outpath : str
            Path to save plot.
        """           
        x_kinem1 = x_kinem + "1"
        x_kinem2 = x_kinem + "2"
        y_kinem1 = y_kinem + "1"
        y_kinem2 = y_kinem + "2"
        
        x_vals = self.apply_mask_get_data(x_kinem1, lep_selection_type=lep_selection_type, weave=True)
        y_vals = self.apply_mask_get_data(y_kinem2, lep_selection_type=lep_selection_type, weave=True)
        if run_over_only_n_evts != -1:
            x_vals = x_vals[:run_over_only_n_evts]
            y_vals = y_vals[:run_over_only_n_evts]

        # A special case to make comparison of (delta_phi vs. delta_theta) easy with (delta_phi vs. delta_eta).
        if (x_kinem1[:-1] == "delta_theta") and (y_kinem1[:-1] == "delta_phi"):
            x_vals *= -1

        #--- Make plots ---#
        if (ax is None):
            f, ax = plt.subplots(figsize=(12.8, 9.6))
        
        x_2D_bins, x_2D_bin_width = make_binning_array(x_bin_limits)
        y_2D_bins, y_2D_bin_width = make_binning_array(y_bin_limits) 
        
        # Plot 1: dphi vs. deta
        if lep_selection_type not in ["1","2"]:
            x_label = label_LaTeX_dict[x_kinem1]["independent_label"]
            y_label = label_LaTeX_dict[y_kinem1]["independent_label"]
        else:
            x_label_1 = label_LaTeX_dict[x_kinem1]["label"]
            x_label_2 = label_LaTeX_dict[x_kinem2]["label"]
            y_label_1 = label_LaTeX_dict[y_kinem1]["label"]
            y_label_2 = label_LaTeX_dict[y_kinem2]["label"]

            x_unit = label_LaTeX_dict[x_kinem2]["units"]
            y_unit = label_LaTeX_dict[y_kinem2]["units"]

            def prep_2D_label(label_1, label_2, unit, bin_width):
                label = "{},   {}".format(label_1, label_2)
                label += "\n" + "(bin width: {:.2E})".format(bin_width)
                if len(unit) > 0:
                    label =label.rstrip(")")
                    label += " {})".format(unit)  
                return label

            x_label = prep_2D_label(x_label_1, x_label_2, x_unit, x_2D_bin_width)
            y_label = prep_2D_label(y_label_1, y_label_2, y_unit, y_2D_bin_width)

        ax.set_xlabel(x_label)#, fontsize=label_size)
        ax.set_ylabel(y_label)#, fontsize=label_size)
        
        if len(title) > 0:
            title += "\n"
        cuts = "Selection type = {}:\n".format(lep_selection_type) 
        cuts += r"{}".format(self.cuts)
    #         ax.set_title(title + cuts)#, fontsize=label_size)
    
        # Stats: 
    #         stat_text_x = 0.1
        stat_text_x = 0.2
        stat_text_y = 0.83
        
        stats_ls = get_stats_2Dhist(x_vals, y_vals)
        leg_label = cuts + "\n" + make_stats_legend_for_2dhist(stats_ls)
        ax.text(stat_text_x, stat_text_y, leg_label, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        
        newcmp = change_cmap_bkg_to_white('rainbow')
        bin_vals, x_bin_edges, y_bin_edges, im = ax.hist2d(x_vals, y_vals, bins=[x_2D_bins, y_2D_bins], cmap=newcmp)
        plt.colorbar(im, ax=ax)

    def plot_kinem_genrec_comparison(self,
                              kinem_gen, kinem_rec, 
                              lep_selection_type="independent", 
                              x_limits=[-1, 1],
                              bin_limits=[-0.5,0.5,0.5], 
                              run_over_only_n_evts=-1,
                              ax=None, ax_ratio=None, log_scale=False
                              ):
            """
            FIXME: Need to implement under/overflow bins!

            Plots differences in kinematics. Includes a ratio plot at the bottom.
            
            Parameters
            ----------
            kinem_gen : str
                The generator-level kinematical variable from the column of DF.
                E.g. "genLep_pt1" or "genLep_eta2", etc.
            kinem_rec : str
                The reconsructed-level kinematical variable from the column of DF.
                E.g. "pT1" or "eta2", etc.
            x_limits : 2-element list
                The x_min and x_max to show along x_axis, for viewing purposes: [x_min, x_max]
            bin_limits : 3-element list
                [first_bin_left_edge, last_bin_right_edge, bin_width]
            ax : axes object
                An external axes object to pass in, on which the main plot will be plotted.
                If None, a default one will get created.
            ax_ratio : axes object
                An external axes object to pass in, on which the ratio plot will be plotted.
                If None, a default one will get created.
            log_scale : bool
                Set y-axis to log scale.
            """

            df = self.binned_df

            x_bin_arr, x_bin_width = make_binning_array(bin_limits)
            x_bin_centers_arr = shift_binning_array(x_bin_arr)

            #--- Get data ---#
            # Make sure you're plotting gen and rec of same lepton.
            if kinem_gen[-1] != kinem_rec[-1]:
                print("[WARNING] It seems you are plotting lepton 1 kinematics vs. lepton 2's!")
            if kinem_gen[-1] != lep_selection_type:
                err_msg = "[ERROR] You want to plot lep{} kinematics but you specified lep{} selection type.".format(kinem_gen[-1], lep_selection_type)
                raise ValueError(err_msg)

            # Chooses either lep1 or lep2. 
            data_rec = self.apply_mask_get_data(kinem_rec, lep_selection_type)
            data_gen = self.apply_mask_get_data(kinem_gen, lep_selection_type)

            if run_over_only_n_evts != -1:
                data_rec = data_rec[:run_over_only_n_evts]
                data_gen = data_gen[:run_over_only_n_evts]
                
            # Gen and Reco stats:
            stats_ls_gen = get_stats_1Dhist(data_gen)
            stats_ls_rec = get_stats_1Dhist(data_rec)

            #----------------#
            #--- Plot It. ---#
            #----------------#
            if (ax is None) or (ax_ratio is None):
                fig = plt.figure(figsize=(10,8))
                ax = fig.add_axes([0.17,0.33,0.825,0.54])  # [low_left_corner_x, low_left_corner_y, width, height]
                ax_ratio = fig.add_axes([0.17,0.12,0.825,0.20])

            leg_label_gen = make_stats_legend_for_1dhist(stats_ls_gen)
            leg_label_rec = make_stats_legend_for_1dhist(stats_ls_rec)

            leg_label_gen = leg_label_gen.replace("Entries", "GEN Entries")
            leg_label_rec = leg_label_rec.replace("Entries", "REC Entries")

            hist_bin_vals_gen, bin_edges_gen, _ = ax.hist(data_gen, bins=x_bin_arr, histtype='step', color='g', lw=2, label=leg_label_gen)
            hist_bin_vals_rec, bin_edges_rec, _ = ax.hist(data_rec, bins=x_bin_arr, histtype='step', color='b', label=leg_label_rec)

            hist_bin_vals_gen_modified = hist_bin_vals_gen.copy()
            hist_bin_vals_gen_modified[hist_bin_vals_gen == 0] = 0.0000000001
            ratio_vals = (hist_bin_vals_rec - hist_bin_vals_gen) / hist_bin_vals_gen_modified

            ax_ratio.errorbar(x_bin_centers_arr, ratio_vals, xerr=x_bin_width/2, color='black', ms=0, capsize=0, mew=0, ecolor='black', drawstyle='steps-mid', alpha=1)

            # Pretty up the plot. 
            ax_ratio.grid(which='major',axis='x')
            ax_ratio.grid(which='major',axis='y', ls='-')
            ax_ratio.grid(which='minor',axis='y')

            # Hide first tick label on y-axis, since it overlaps with ratio plot's tick label.
            a=ax.get_yticks().tolist()
            a[0]=''
            ax.set_yticklabels(a)

            # Hide main plot's x tick labels.
            plt.setp(ax.get_xticklabels(), visible=False)

            # Only show a few of the tick labels on ratio plot.
            n_tick_labels = 5
            ax_ratio.yaxis.set_major_locator(plt.MaxNLocator(n_tick_labels))
            ax_ratio.axhline(c='r', lw=2, ls='-')

            textsize_legend = 10
            textsize_axislabels = 12
            textsize_title = 12

            unit_gen = label_LaTeX_dict[kinem_gen]["units"]
            unit_rec = label_LaTeX_dict[kinem_rec]["units"]
            
            y_label = hist_y_label(x_bin_width, unit_gen)
            
            # Remember that ax shouldn't have an x-label; it's covered by the ax_ratio. 
            ax.set_ylabel(y_label, fontsize=textsize_axislabels)
    #             ax.set_title(r"Selection: {}".format(self.cuts), fontsize=textsize_title)
            
            textbox_text = "Selection type = {}:\n".format(lep_selection_type) + self.cuts
            ax.text(0.025, 0.83, textbox_text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        
            plt.minorticks_on()

            x_min = x_limits[0]
            x_max = x_limits[1]
            ax.set_xlim([x_min, x_max])
            ax_ratio.set_xlim([x_min, x_max])
            ax_ratio.set_ylim([-0.12, 0.12])

            x_label_gen = label_LaTeX_dict[kinem_gen]["label"]
            x_label_rec = label_LaTeX_dict[kinem_rec]["label"]
            
            x_label = r"{}, {}".format(x_label_gen, x_label_rec)
            if len(unit_gen) > 0:
                x_label += " [{}]".format(unit_gen)
            ax_ratio.set_xlabel(x_label, fontsize=textsize_axislabels)
            ax_ratio.set_ylabel(r'$\frac{ (\mathrm{REC} - \mathrm{GEN}) }{ \mathrm{GEN} }$', fontsize=textsize_axislabels*1.5)

            y_max = max(max(hist_bin_vals_gen), max(hist_bin_vals_rec))
            ax.set_ylim([0, y_max*1.4])

            ax.legend(loc='upper right', framealpha=0.9, fontsize=textsize_legend)#, horizontalalignment='right')
            
            return ax, ax_ratio

    def plot_1D_kinematics(self, kinem="", lep_selection_type="independent", x_limits=[0,0], bin_limits=[0,0,0], run_over_only_n_evts=-1, ax=None, x_label="", y_label="", title="", y_max=-1, log_scale=False, iter_gaus=(False, 0)):
        """
        Make a histogram of a kinematical variable in the DF.  
        FIXME: 
            - Currently only works with kinem variables that end with `1` or `2`. 
            [ ] Eventually plot other quantities like massZ, etc.
        
        Parameters
        ----------
        kinem : str
            The full name of the kinematical branch/column ("delta_eta1", "delta_R2", "massZ", etc.)
        lep_selection_type : int
            Indicates the kind of selection to apply on leptons:
                lep_selection_type = "1"    -- Plot kinem values in which only lepton1 passed selections.
                lep_selection_type = "2"    -- Plot kinem values in which only lepton2 passed selections.
                lep_selection_type = "both" -- Plot kinem values in which lepton1 AND lepton2 passed selections.
                lep_selection_type = "either" -- Plot kinem values in which lepton1 OR lepton2 passed selections.
        x_limits : 2-element list
            The x_min and x_max to show along x_axis, for viewing purposes: [x_min, x_max]
        bin_limits : 3-element list
            [first_bin_left_edge, last_bin_right_edge, bin_width]   
        run_over_only_n_evts : int
            Number of events to put into histogram. 
        ax : axes object
            An external axes object to pass in, on which the main plot will be plotted.
            If None, a default one will get created.
        x_label : str
            Override the default x-label that would be produced.
        y_label : str
            Override the default y-label that would be produced.
        title : str
            Override the default title that would be produced.
        y_max : float
            The max y-value on vertical axis for viewing purposes.
            If set to -1, then automatic detection is used. 
        log_scale : bool
            Use log-scale on vertical axis.
        iter_gaus : 2-tuple
            If True, perform an iterative gaussian fit on the core N times.
            Syntax: (switch, N)
        """
        #--- Consistency checks ---#
        if ("BS" in kinem) or ("PV" in kinem):
            if self.d0_type not in kinem:
                err_msg = "[ERROR] The kinematic '{}' was specified but d0_type is '{}'.\nStopping now".format(kinem, self.d0_type)
                raise ValueError(err_msg)
                
        df = self.binned_df
        
        #--- Get data ---#
        data = self.apply_mask_get_data(kinem=kinem, lep_selection_type=lep_selection_type, weave=True)  # kinem must be a full name
            
        if run_over_only_n_evts != -1:
            data = data[:run_over_only_n_evts]
            
        if ax is None:
            # Axes doesn't exist yet. Make it.
            fig, ax = plt.subplots(figsize=(12.8,9.6))
            
        if bin_limits == [0,0,0]:
            # No bin limits specified, so use default binning for this kinematical variable.
            bin_limits = label_LaTeX_dict[kinem]["default_bin_limits"]
             
        if x_limits == [0,0]:
            # No x-limits specified, so use default x-limits for this kinematical variable.
            x_limits = label_LaTeX_dict[kinem]["default_x_limits"]
            
        unit = label_LaTeX_dict[kinem]["units"]
            
        x_bins, binw = make_binning_array(bin_limits)

        if len(x_label) == 0:
                # Both muons are being chosen. Change the labels.
            key = "label" if lep_selection_type in ["1","2"] else "independent_label"
            x_label = label_LaTeX_dict[kinem][key]
            if len(unit) > 0:
                x_label += " [{}]".format(unit)
            
        if len(y_label) == 0:
            # User didn't specify label, so make it.
            y_label = hist_y_label(binw, unit)
            
    #         if len(title) == 0:
    #             title = "Selections:\n" + self.cuts
        
        ax, bin_vals, bin_edges, stats = make_1D_dist(ax, data, x_limits, x_bins, 
                                                      x_label, y_label, title, y_max=-1, log_scale=False)
        
        # Nested dictionaries.
        # Initializing this particular kinematic bin's dictionary of stats.
        self.stats_dict[kinem] = {}
        self.stats_dict[kinem]['hist_stats'] = stats
        self.stats_dict[kinem]['bin_vals'] = bin_vals
        self.stats_dict[kinem]['bin_edges'] = bin_edges

        textbox_text = r"Selection type = {}:".format(lep_selection_type) + "\n"
        if lep_selection_type in ["1","2"]:
            textbox_text = textbox_text.replace("= ", r"$\mu$")
        textbox_text += self.cuts
        ax.text(0.03, 0.87, textbox_text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        
        if iter_gaus[0]:
            # Do iterative fitting procedure.
            # Produce a dictionary of stats for these fits.
            best_guess_coeff = max(bin_vals)  # The Gaus coeff is usually around the max height of Gaus curve.
            best_guess_mean = stats[1]
            best_guess_stdev = stats[3]
            stats_dict, ax = iterative_fit_gaus(iter_gaus[1], bin_edges, bin_vals, 
                                                param_guess=[best_guess_coeff, best_guess_mean, best_guess_stdev],
                                                param_bounds=([0,-1000,-100], [999999999,1000,100]),
                                                ax=ax, draw_on_axes=True, verbose=self.verbose)
            # Use plotted kinem as the key for this dict of stats. 
            self.stats_dict[kinem]['fit_stats'] = stats_dict         
            
    def count_in_pT_eta_region_exclusive(self):
        """
        Return the number of muons in the DF subset which pass the eta and pT cuts. 
        'Exclusive' because only those muons which pass such cuts are counted 
        (as opposed to counting both muons from event just because one pass the cuts).
        """
        pass

class KinematicBinOLD():
    """FIXME: Created in the spring of 2020, back when I was plotting with plt.
    
    My first class built around pandas and vaex.
    """

    def __init__(self, 
                df_original, 
    #                  inpath_dataframe, 
                 n_evts, 
                 massZ_cut_ls, eta_cut_ls, pT_cut_ls, d0q_cut_ls, d0_type="BS", dR_cut=0.02, 
                 use_ptotal_instead=False, verbose=False):
        """
        Pass in a DataFrame (DF) and specify the eta and pT cuts to create a subset of DF.
        
        Parameters
        ----------
        df_original : pandas.DataFrame
            ROOT file converted into a DataFrame. Columns are branches. Rows are events.
        n_evts : int
            Number of events to search over - not guaranteed to find this many events! 
            Use '-1' to loop over all events in the df. 
        eta_cut_ls : list or array-like of floats
            A list of [eta_min, eta_max]. Example: [0.9, 1.8]
        pT_cut_ls : list or array-like of floats
            A list of [pT_min, pT_max]. Example: [5, 20]
        d0q_cut_ls : list or array-like of floats
            A list of values of d0*charge to cut on: [d0q_min, d0q_max]. E.g. [-0.01, 0.01]
        d0_type : str
            Which d0 to cut on: "BS" or "PV"
        use_ptotal_instead : bool
            Cut on total momentum instead of pT. 
            In most places, pT could stand for either p or pT (depending on 'use_ptotal_instead').
                Just be careful because not all places are adapted for p!
        dR_cut : float
            A cut to save events in which muon1 and muon2 both have dR < dR_cut.
        verbose : bool
            If True, get debug info and see where you are while the code runs.
            
        NOTE:
            The methods further down are more developed than the methods closer to __init__(). 
            Therefore, clean up and consolidate the earlier methods. 
        """
        if n_evts == -1:
            n_evts = len(df_original)
        df = df_original.loc[:n_evts].copy()    # Original DF. 
        
    #         if n_evts > 0:
    #             df = pd.read_hdf(inpath_dataframe, stop=n_evts)
    #         elif n_evts == -1:
    #             # Read in entire DF.
    #             df = pd.read_hdf(inpath_dataframe)
    #         else: 
    #             raise ValueError("[ERROR] self.n_evts={} specified incorrectly.".format(self.n_evts))
            # Apparently this way is even slower than just copying the DF...
    
        self.n_evts_asked_for = n_evts
        self.n_evts_found = -999
        self.kinem_vals_after_selection = {}
        self.stats_dict = {}
        self.binned_df = None  # Will become the collection of events in which mu1 or mu2 passes all cuts.
        self.verbose = verbose
        
        self.cuts = ""
        self.sorted_cut_ls = []
        self.cut_dict = {}
        self.use_ptotal_instead = use_ptotal_instead
        self.p_str = "p" if (self.use_ptotal_instead) else "pT"
        self.p_str_latex = r"$p^{\mathrm{REC}}$" if (self.use_ptotal_instead) else r"$p_{T}^{\mathrm{REC}}$"

        self.massZ_min = massZ_cut_ls[0]
        self.massZ_max = massZ_cut_ls[1]        
        self.eta_min   = eta_cut_ls[0]
        self.eta_max   = eta_cut_ls[1]
        self.pT_min    = pT_cut_ls[0] 
        self.pT_max    = pT_cut_ls[1]         
        self.d0q_min   = d0q_cut_ls[0]
        self.d0q_max   = d0q_cut_ls[1]
        self.d0_type   = d0_type
        self.dR_cut    = dR_cut   
                
        self.apply_initial_cuts(df, verbose)
        
    def apply_initial_cuts(self, df, verbose):
        """
        Creates a subset of the original DataFrame in which initial cuts are applied.
        Cuts:
            pT
            eta
            d0
            massZ
            dR
        """
        # Cuts:
        eta_min   = self.eta_min  
        eta_max   = self.eta_max  
        pT_min    = self.pT_min   
        pT_max    = self.pT_max   
        d0q_min    = self.d0q_min
        d0q_max    = self.d0q_max
        dR_cut    = self.dR_cut   
        massZ_min = self.massZ_min
        massZ_max = self.massZ_max

        #--------------------------------------#
        #--- Save and manipulate variables. ---#
        #--------------------------------------#
        # GEN info. 
        eta1_gen_ser = df['genLep_eta1']
        eta2_gen_ser = df['genLep_eta2']
        phi1_gen_ser = df['genLep_phi1']
        phi2_gen_ser = df['genLep_phi2']
        pT1_gen_ser  = df['genLep_pt1']
        pT2_gen_ser  = df['genLep_pt2']

        # RECO info.
        eta1_rec_ser = df['eta1']
        eta2_rec_ser = df['eta2']
        phi1_rec_ser = df['phi1']
        phi2_rec_ser = df['phi2']
        pT1_rec_ser  = df['pT1']
        pT2_rec_ser  = df['pT2']
        
        # Store other variables.
        df.loc[:,'genLep_theta1'] = theta1_gen_ser = pseudorap2theta(eta1_gen_ser)
        df.loc[:,'genLep_theta2'] = theta2_gen_ser = pseudorap2theta(eta2_gen_ser)
        df.loc[:,'theta1'] = theta1_rec_ser = pseudorap2theta(eta1_rec_ser)
        df.loc[:,'theta2'] = theta2_rec_ser = pseudorap2theta(eta2_rec_ser)
        
        df.loc[:,'genLep_p1'] = pT1_gen_ser / np.sin( theta1_gen_ser )
        df.loc[:,'genLep_p2'] = pT2_gen_ser / np.sin( theta2_gen_ser )
        df.loc[:,'p1'] = pT1_rec_ser / np.sin( theta1_rec_ser )  # Total momentum
        df.loc[:,'p2'] = pT2_rec_ser / np.sin( theta2_rec_ser )  # Total momentum
                        
        df.loc[:,'delta_eta1'] = deta1_ser = eta1_rec_ser - eta1_gen_ser
        df.loc[:,'delta_eta2'] = deta2_ser = eta2_rec_ser - eta2_gen_ser
        df.loc[:,'delta_theta1'] = dtheta1_ser = theta1_rec_ser - theta1_gen_ser
        df.loc[:,'delta_theta2'] = dtheta2_ser = theta2_rec_ser - theta2_gen_ser

        # Remember that delta_phi requires special treatment:
        # -pi < delta_phi < pi
        df.loc[:,'delta_phi1'] = dphi1_ser = calc_dphi(phi1_rec_ser, phi1_gen_ser)
        df.loc[:,'delta_phi2'] = dphi2_ser = calc_dphi(phi2_rec_ser, phi2_gen_ser)

        df.loc[:,'delta_R1'] = dR1_ser = calc_dR(deta1_ser, dphi1_ser)
        df.loc[:,'delta_R2'] = dR2_ser = calc_dR(deta2_ser, dphi2_ser)
        df.loc[:,'delta_Rtheta1'] = dRtheta1_ser = calc_dR(dtheta1_ser, dphi1_ser)
        df.loc[:,'delta_Rtheta2'] = dRtheta2_ser = calc_dR(dtheta2_ser, dphi2_ser)
        
        df.loc[:,'delta_pT1'] = dpT1 = pT1_rec_ser - pT1_gen_ser
        df.loc[:,'delta_pT2'] = dpT2 = pT2_rec_ser - pT2_gen_ser
        
        df.loc[:,'delta_pToverpT1'] = dpT1 / pT1_gen_ser
        df.loc[:,'delta_pToverpT2'] = dpT2 / pT2_gen_ser
        
        df.loc[:,'delta_pToverRecpT1'] = dpT1 / pT1_rec_ser
        df.loc[:,'delta_pToverRecpT2'] = dpT2 / pT2_rec_ser     
        
        df.loc[:,'d0BSq1'] = df['d0BS1'] * df['Id1'].replace(13,-1).replace(-13,1)
        df.loc[:,'d0BSq2'] = df['d0BS2'] * df['Id2'].replace(13,-1).replace(-13,1)
        df.loc[:,'d0PVq1'] = df['d0PV1'] * df['Id1'].replace(13,-1).replace(-13,1)
        df.loc[:,'d0PVq2'] = df['d0PV2'] * df['Id2'].replace(13,-1).replace(-13,1)
        
        # Create masks.
        self.mask_massZ = mask_massZ = self.get_mask_massZ(df)
        
        self.mask_eta1, self.mask_eta2 = mask_eta1, mask_eta2 = self.get_mask_eta(df)
        self.mask_pT1,  self.mask_pT2  = mask_pT1,  mask_pT2  = self.get_mask_pT(df)
        self.mask_d0q1, self.mask_d0q2 = mask_d0q1, mask_d0q2 = self.get_mask_d0q(df)    
        self.mask_dR1,  self.mask_dR2  = mask_dR1,  mask_dR2  = self.get_mask_dR(df) 
        
        # Combine masks.
        self.mask_kinembin_lep1 = mask_massZ & mask_dR1 & mask_eta1 & mask_pT1 & mask_d0q1
        self.mask_kinembin_lep2 = mask_massZ & mask_dR2 & mask_eta2 & mask_pT2 & mask_d0q2

        # Keep all events in which either muon1 passed all selections or muon2 passed all. 
        self.all_masks = self.mask_kinembin_lep1 | self.mask_kinembin_lep2
        
        # Apply masks and update DataFrame.
    #        self.binned_df = df[self.all_masks].copy()
        self.binned_df = df[self.all_masks]
    #        del df
        self.n_evts_found = len(self.binned_df)
        
        # The cut_dict has been filled. Now convert it to an ordered list (alphabetically).
        self.sorted_cut_ls = sorted_cut_ls = [value for (key, value) in sorted(self.cut_dict.items())]
        self.cuts = combine_cut_list(sorted_cut_ls)
        
        if (self.verbose): 
            perc = self.n_evts_found / float(self.n_evts_asked_for) * 100.
            print("[INFO] Events found: {} ({:.3f}% of total events scanned)".format(self.n_evts_found, perc))
            print(r"using cuts: {}".format(self.cuts) + "\n")
      
    def get_mask_d0q(self, df):
        d0_type = self.d0_type
        if d0_type == "PV":
            mask_d0q1 = (self.d0q_min < df['d0PVq1']) & (df['d0PVq1'] < self.d0q_max)
            mask_d0q2 = (self.d0q_min < df['d0PVq2']) & (df['d0PVq2'] < self.d0q_max)
        elif d0_type == "BS":
            mask_d0q1 = (self.d0q_min < df['d0BSq1']) & (df['d0BSq1'] < self.d0q_max)
            mask_d0q2 = (self.d0q_min < df['d0BSq2']) & (df['d0BSq2'] < self.d0q_max)
            
        cuts_d0q = r"$%.3f < d_{0}^{\mathrm{%s}}*q(\mu) < %.3f$" % (self.d0q_min, self.d0_type, self.d0q_max)
        key = "d0{}q".format(self.d0_type)
        self.cut_dict[key] = cuts_d0q
        return mask_d0q1, mask_d0q2
        
    def get_mask_pT(self, df):
        if (self.use_ptotal_instead):
            mask_pT1 = (self.pT_min < df['p1']) & (df['p1'] < self.pT_max) 
            mask_pT2 = (self.pT_min < df['p2']) & (df['p2'] < self.pT_max)
        else:
            mask_pT1 = (self.pT_min < df['pT1']) & (df['pT1'] < self.pT_max) 
            mask_pT2 = (self.pT_min < df['pT2']) & (df['pT2'] < self.pT_max)
        
        cuts_p = r"$%d <$ %s $< %d$ GeV" % (self.pT_min, self.p_str_latex, self.pT_max)  # The string brings in its own '$'.
        self.cut_dict[self.p_str] = cuts_p
        return mask_pT1, mask_pT2

    def get_mask_eta(self, df):
        mask_eta1 = (self.eta_min < abs(df['eta1'])) & (abs(df['eta1']) < self.eta_max)
        mask_eta2 = (self.eta_min < abs(df['eta2'])) & (abs(df['eta2']) < self.eta_max)   
        
        cuts_eta = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$" % (self.eta_min, self.eta_max)
        self.cut_dict["eta"] = cuts_eta
        return mask_eta1, mask_eta2

    def get_mask_dR(self, df):
        mask_dR1 = (df['delta_R1'] < self.dR_cut)
        mask_dR2 = (df['delta_R2'] < self.dR_cut)
        
        cuts_dR = r"$\Delta R < %.3f$" % (self.dR_cut)
        self.cut_dict["delta_R"] = cuts_dR
        return mask_dR1, mask_dR2
    
    def get_mask_massZ(self, df):
        mask_massZ = (self.massZ_min < df['massZ']) & (df['massZ'] < self.massZ_max)
        
        cuts_massZ = r"$%.1f < m_{\mu\mu} < %.1f$ GeV" % (self.massZ_min, self.massZ_max)
        self.cut_dict["massZ"] = cuts_massZ
        return mask_massZ    
    
    def apply_mask_get_data(self, kinem, lep_selection_type="", weave=False):
        """
        Return the kinematic lepton data which pass selections in this kinematic bin.
        E.g. Apply a boolean mask for event selection and retrieve all "delta_R1" values.

        Parameters
        ----------
        kinem : str
            A complete branch name in the DataFrame or root file. 
                E.g. "pT1", "genLep_pt2", "massZ"
        lep_selection_type : int
            The lepton's mask you want to apply.
            
            lep_selection_type = "1" -- get kinem values in which muon1 passes all selection criteria 
                                       (muon 2 may or may not pass selections).
            lep_selection_type = "2" -- get kinem values in which muon2 passes all selection criteria.
            lep_selection_type = "both" -- BOTH muons must pass selections to get data from event.
            lep_selection_type = "either" -- Either muon1, or muon2, or both must pass selections to get data from event.
            lep_selection_type = "independent" -- Get kinematic values for muon1 and muon2, with no restriction on which event they came from.
                Note: If kinem ends in "1" or "2", then the kinematic values of the other lepton are automatically grabbed.
        weave : bool
            Weave lep1 kinematic values and lep2 kinematic values together, so that slicing doesn't just give lep1. 
            Only relevant for "independent" selection.
        
        Returns
        -------
        kinem_vals : array 
            Kinematic values with chosen mask applied. 
            All values satisfy selection criteria for this kinematic bin.
        """
        df = self.binned_df
        mask1 = self.mask_kinembin_lep1
        mask2 = self.mask_kinembin_lep2
        
        if lep_selection_type == "1":
            # Only select events in which muon1 passes selections.
            mask = mask1
        elif lep_selection_type == "2":
            # Only select events in which muon2 passes selections.
            mask = mask2
        elif lep_selection_type == "both":
            # Only select events in which BOTH muon1 and muon2 pass selections.
            mask = mask1 & mask2
        elif lep_selection_type == "either":
            # Only select events in which BOTH muon1 and muon2 pass selections.
            mask = mask1 | mask2
        elif lep_selection_type == "independent":
            # Go through all muons in all events, without regard for other muon. 
            if kinem[-1] in ["1", "2"]:
                # Lep1 (or lep2) kinematic detected. Go find the other lepton's kinematic values.
                # FIXME: the variable massZ_vtxChi2 will be wrongly caught by this 'if' statement!
                kinem1 = kinem[:-1] + "1"
                kinem2 = kinem[:-1] + "2"
                kinem_vals1 = df[kinem1][mask1].values
                kinem_vals2 = df[kinem2][mask2].values
            else:
                # The kinematic doesn't depend on lep1 or lep2, like: massZ, GENmass2l, etc.
                kinem_vals1 = df[kinem][mask1].values
                kinem_vals2 = df[kinem][mask2].values
            
            if (weave):
                # Weave values together so that when slicing (like [:5]), you don't just grab kinem_vals1. 
                kinem_vals = np.array( weave_lists(kinem_vals1, kinem_vals2) )
            else: 
                kinem_vals = np.append(kinem_vals1, kinem_vals2)
        
            return kinem_vals
        
        else: 
            raise ValueError("[ERROR] `lep_selection_type` was not specified properly. Stopping now.")
        
        # A selection, other than "independent" was chosen.
        kinem_vals = df[kinem][mask].values
        
        return kinem_vals
    
    def make_2D_plot(self, 
                     x_kinem, y_kinem, 
                     x_bin_limits=[0, 1, 0.1], y_bin_limits=[0, 1, 0.1],
                     lep_selection_type="",
                     run_over_only_n_evts=-1, 
                     title="",
                     exclusive=True,
                     save_plot=False, save_as_png=False, outpath="",
                     ax=None):
        """
        Make a 2D plot. Two examples:
            (1) dphi vs. deta  
            (2) dphi vs. dtheta
        User can specify the binning along either axis. 
        This method plots only the muons which pass the selection 
        (as opposed to taking any event in which at least 1 muon pass kinematic bin criteria).
        
        Parameters
        ----------
        x_kinem : str
            The PARTIAL name of the kinematical variable to be plotted along x-axis. 
            Only works for kinematics which end with '1' or '2'.
            - Example: x_kinem="delta_theta" (for which there are two branches: "delta_theta1", "delta_theta2")
        y_kinem : str
            The PARTIAL name of the kinematical variable to be plotted along y-axis. 
            Only works for kinematics which end with '1' or '2'.
            - Example: y_kinem="delta_eta" (for which there are two branches: "delta_eta1", "delta_eta2")
        x_bin_limits : list or array-like of floats
            The bin limits on the horizontal axis. [bin_min_left_edge, bin_max_right_edge, bin_width]
            - Example: [-2.5, 2.5, 0.1]
        y_bin_limits : list or array-like of floats
            The bin limits on the vertical axis. [bin_min_left_edge, bin_max_right_edge, bin_width]
            - Example: [-2.5, 2.5, 0.1]
        lep_selection_type : str
            What kind of selection to perform on the leptons. Choices:
            #UPDATE
        run_over_only_n_evts : int
            Number of events to plot. Use '-1' to use all events in this kinembin.
        title : str
            Alternate title to put on plot. Overrides the default one made in this method.
        exclusive : bool
            Means "only put muons which passed all selections in this plot".
            FIXME: It must be set to True for now...
        save_plot : bool
            If True, save the plot as a pdf and possibly a png.
        save_as_png : bool
            If True, save the plot as a png.
        outpath : str
            Path to save plot.
        """           
        x_kinem1 = x_kinem + "1"
        x_kinem2 = x_kinem + "2"
        y_kinem1 = y_kinem + "1"
        y_kinem2 = y_kinem + "2"
        
        x_vals = self.apply_mask_get_data(x_kinem1, lep_selection_type=lep_selection_type, weave=True)
        y_vals = self.apply_mask_get_data(y_kinem2, lep_selection_type=lep_selection_type, weave=True)
        if run_over_only_n_evts != -1:
            x_vals = x_vals[:run_over_only_n_evts]
            y_vals = y_vals[:run_over_only_n_evts]

        # A special case to make comparison of (delta_phi vs. delta_theta) easy with (delta_phi vs. delta_eta).
        if (x_kinem1[:-1] == "delta_theta") and (y_kinem1[:-1] == "delta_phi"):
            x_vals *= -1

        #--- Make plots ---#
        if (ax is None):
            f, ax = plt.subplots(figsize=(12.8, 9.6))
        
        x_2D_bins, x_2D_bin_width = make_binning_array(x_bin_limits)
        y_2D_bins, y_2D_bin_width = make_binning_array(y_bin_limits) 
        
        # Plot 1: dphi vs. deta
        if lep_selection_type not in ["1","2"]:
            x_label = label_LaTeX_dict[x_kinem1]["independent_label"]
            y_label = label_LaTeX_dict[y_kinem1]["independent_label"]
        else:
            x_label_1 = label_LaTeX_dict[x_kinem1]["label"]
            x_label_2 = label_LaTeX_dict[x_kinem2]["label"]
            y_label_1 = label_LaTeX_dict[y_kinem1]["label"]
            y_label_2 = label_LaTeX_dict[y_kinem2]["label"]

            x_unit = label_LaTeX_dict[x_kinem2]["units"]
            y_unit = label_LaTeX_dict[y_kinem2]["units"]

            def prep_2D_label(label_1, label_2, unit, bin_width):
                label = "{},   {}".format(label_1, label_2)
                label += "\n" + "(bin width: {:.2E})".format(bin_width)
                if len(unit) > 0:
                    label =label.rstrip(")")
                    label += " {})".format(unit)  
                return label

            x_label = prep_2D_label(x_label_1, x_label_2, x_unit, x_2D_bin_width)
            y_label = prep_2D_label(y_label_1, y_label_2, y_unit, y_2D_bin_width)

        ax.set_xlabel(x_label)#, fontsize=label_size)
        ax.set_ylabel(y_label)#, fontsize=label_size)
        
        if len(title) > 0:
            title += "\n"
        cuts = "Selection type = {}:\n".format(lep_selection_type) 
        cuts += r"{}".format(self.cuts)
    #         ax.set_title(title + cuts)#, fontsize=label_size)
    
        # Stats: 
    #         stat_text_x = 0.1
        stat_text_x = 0.2
        stat_text_y = 0.83
        
        stats_ls = get_stats_2Dhist(x_vals, y_vals)
        leg_label = cuts + "\n" + make_stats_legend_for_2dhist(stats_ls)
        ax.text(stat_text_x, stat_text_y, leg_label, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        
        newcmp = change_cmap_bkg_to_white('rainbow')
        bin_vals, x_bin_edges, y_bin_edges, im = ax.hist2d(x_vals, y_vals, bins=[x_2D_bins, y_2D_bins], cmap=newcmp)
        plt.colorbar(im, ax=ax)

    def plot_kinem_genrec_comparison(self,
                              kinem_gen, kinem_rec, 
                              lep_selection_type="independent", 
                              x_limits=[-1, 1],
                              bin_limits=[-0.5,0.5,0.5], 
                              run_over_only_n_evts=-1,
                              ax=None, ax_ratio=None, log_scale=False
                              ):
            """
            FIXME: Need to implement under/overflow bins!

            Plots differences in kinematics. Includes a ratio plot at the bottom.
            
            Parameters
            ----------
            kinem_gen : str
                The generator-level kinematical variable from the column of DF.
                E.g. "genLep_pt1" or "genLep_eta2", etc.
            kinem_rec : str
                The reconsructed-level kinematical variable from the column of DF.
                E.g. "pT1" or "eta2", etc.
            x_limits : 2-element list
                The x_min and x_max to show along x_axis, for viewing purposes: [x_min, x_max]
            bin_limits : 3-element list
                [first_bin_left_edge, last_bin_right_edge, bin_width]
            ax : axes object
                An external axes object to pass in, on which the main plot will be plotted.
                If None, a default one will get created.
            ax_ratio : axes object
                An external axes object to pass in, on which the ratio plot will be plotted.
                If None, a default one will get created.
            log_scale : bool
                Set y-axis to log scale.
            """

            df = self.binned_df

            x_bin_arr, x_bin_width = make_binning_array(bin_limits)
            x_bin_centers_arr = centers_of_binning_array(x_bin_arr)

            #--- Get data ---#
            # Make sure you're plotting gen and rec of same lepton.
            if kinem_gen[-1] != kinem_rec[-1]:
                print("[WARNING] It seems you are plotting lepton 1 kinematics vs. lepton 2's!")
            if kinem_gen[-1] != lep_selection_type:
                err_msg = "[ERROR] You want to plot lep{} kinematics but you specified lep{} selection type.".format(kinem_gen[-1], lep_selection_type)
                raise ValueError(err_msg)

            # Chooses either lep1 or lep2. 
            data_rec = self.apply_mask_get_data(kinem_rec, lep_selection_type)
            data_gen = self.apply_mask_get_data(kinem_gen, lep_selection_type)

            if run_over_only_n_evts != -1:
                data_rec = data_rec[:run_over_only_n_evts]
                data_gen = data_gen[:run_over_only_n_evts]
                
            # Gen and Reco stats:
            stats_ls_gen = get_stats_1Dhist(data_gen)
            stats_ls_rec = get_stats_1Dhist(data_rec)

            #----------------#
            #--- Plot It. ---#
            #----------------#
            if (ax is None) or (ax_ratio is None):
                fig = plt.figure(figsize=(10,8))
                ax = fig.add_axes([0.17,0.33,0.825,0.54])  # [low_left_corner_x, low_left_corner_y, width, height]
                ax_ratio = fig.add_axes([0.17,0.12,0.825,0.20])

            leg_label_gen = make_stats_legend_for_1dhist(stats_ls_gen)
            leg_label_rec = make_stats_legend_for_1dhist(stats_ls_rec)

            leg_label_gen = leg_label_gen.replace("Entries", "GEN Entries")
            leg_label_rec = leg_label_rec.replace("Entries", "REC Entries")

            hist_bin_vals_gen, bin_edges_gen, _ = ax.hist(data_gen, bins=x_bin_arr, histtype='step', color='g', lw=2, label=leg_label_gen)
            hist_bin_vals_rec, bin_edges_rec, _ = ax.hist(data_rec, bins=x_bin_arr, histtype='step', color='b', label=leg_label_rec)

            hist_bin_vals_gen_modified = hist_bin_vals_gen.copy()
            hist_bin_vals_gen_modified[hist_bin_vals_gen == 0] = 0.0000000001
            ratio_vals = (hist_bin_vals_rec - hist_bin_vals_gen) / hist_bin_vals_gen_modified

            ax_ratio.errorbar(x_bin_centers_arr, ratio_vals, xerr=x_bin_width/2, color='black', ms=0, capsize=0, mew=0, ecolor='black', drawstyle='steps-mid', alpha=1)

            # Pretty up the plot. 
            ax_ratio.grid(which='major',axis='x')
            ax_ratio.grid(which='major',axis='y', ls='-')
            ax_ratio.grid(which='minor',axis='y')

            # Hide first tick label on y-axis, since it overlaps with ratio plot's tick label.
            a=ax.get_yticks().tolist()
            a[0]=''
            ax.set_yticklabels(a)

            # Hide main plot's x tick labels.
            plt.setp(ax.get_xticklabels(), visible=False)

            # Only show a few of the tick labels on ratio plot.
            n_tick_labels = 5
            ax_ratio.yaxis.set_major_locator(plt.MaxNLocator(n_tick_labels))
            ax_ratio.axhline(c='r', lw=2, ls='-')

            textsize_legend = 10
            textsize_axislabels = 12
            textsize_title = 12

            unit_gen = label_LaTeX_dict[kinem_gen]["units"]
            unit_rec = label_LaTeX_dict[kinem_rec]["units"]
            
            y_label = hist_y_label(x_bin_width, unit_gen)
            
            # Remember that ax shouldn't have an x-label; it's covered by the ax_ratio. 
            ax.set_ylabel(y_label, fontsize=textsize_axislabels)
    #             ax.set_title(r"Selection: {}".format(self.cuts), fontsize=textsize_title)
            
            textbox_text = "Selection type = {}:\n".format(lep_selection_type) + self.cuts
            ax.text(0.025, 0.83, textbox_text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        
            plt.minorticks_on()

            x_min = x_limits[0]
            x_max = x_limits[1]
            ax.set_xlim([x_min, x_max])
            ax_ratio.set_xlim([x_min, x_max])
            ax_ratio.set_ylim([-0.12, 0.12])

            x_label_gen = label_LaTeX_dict[kinem_gen]["label"]
            x_label_rec = label_LaTeX_dict[kinem_rec]["label"]
            
            x_label = r"{}, {}".format(x_label_gen, x_label_rec)
            if len(unit_gen) > 0:
                x_label += " [{}]".format(unit_gen)
            ax_ratio.set_xlabel(x_label, fontsize=textsize_axislabels)
            ax_ratio.set_ylabel(r'$\frac{ (\mathrm{REC} - \mathrm{GEN}) }{ \mathrm{GEN} }$', fontsize=textsize_axislabels*1.5)

            y_max = max(max(hist_bin_vals_gen), max(hist_bin_vals_rec))
            ax.set_ylim([0, y_max*1.4])

            ax.legend(loc='upper right', framealpha=0.9, fontsize=textsize_legend)#, horizontalalignment='right')
            
            return ax, ax_ratio

    def plot_1D_kinematics(self, kinem="", lep_selection_type="independent", x_limits=[0,0], bin_limits=[0,0,0], run_over_only_n_evts=-1, ax=None, x_label="", y_label="", title="", y_max=-1, log_scale=False, iter_gaus=(False, 0)):
        """
        Make a histogram of a kinematical variable in the DF.  
        FIXME: 
            - Currently only works with kinem variables that end with `1` or `2`. 
            [ ] Eventually plot other quantities like massZ, etc.
        
        Parameters
        ----------
        kinem : str
            The full name of the kinematical branch/column ("delta_eta1", "delta_R2", "massZ", etc.)
        lep_selection_type : int
            Indicates the kind of selection to apply on leptons:
                lep_selection_type = "1"    -- Plot kinem values in which only lepton1 passed selections.
                lep_selection_type = "2"    -- Plot kinem values in which only lepton2 passed selections.
                lep_selection_type = "both" -- Plot kinem values in which lepton1 AND lepton2 passed selections.
                lep_selection_type = "either" -- Plot kinem values in which lepton1 OR lepton2 passed selections.
        x_limits : 2-element list
            The x_min and x_max to show along x_axis, for viewing purposes: [x_min, x_max]
        bin_limits : 3-element list
            [first_bin_left_edge, last_bin_right_edge, bin_width]   
        run_over_only_n_evts : int
            Number of events to put into histogram. 
        ax : axes object
            An external axes object to pass in, on which the main plot will be plotted.
            If None, a default one will get created.
        x_label : str
            Override the default x-label that would be produced.
        y_label : str
            Override the default y-label that would be produced.
        title : str
            Override the default title that would be produced.
        y_max : float
            The max y-value on vertical axis for viewing purposes.
            If set to -1, then automatic detection is used. 
        log_scale : bool
            Use log-scale on vertical axis.
        iter_gaus : 2-tuple
            If True, perform an iterative gaussian fit on the core N times.
            Syntax: (switch, N)
        """
        #--- Consistency checks ---#
        if ("BS" in kinem) or ("PV" in kinem):
            if self.d0_type not in kinem:
                err_msg = "[ERROR] The kinematic '{}' was specified but d0_type is '{}'.\nStopping now".format(kinem, self.d0_type)
                raise ValueError(err_msg)
                
        df = self.binned_df
        
        #--- Get data ---#
        data = self.apply_mask_get_data(kinem=kinem, lep_selection_type=lep_selection_type, weave=True)  # kinem must be a full name
            
        if run_over_only_n_evts != -1:
            data = data[:run_over_only_n_evts]
            
        if ax is None:
            # Axes doesn't exist yet. Make it.
            fig, ax = plt.subplots(figsize=(12.8,9.6))
            
        if bin_limits == [0,0,0]:
            # No bin limits specified, so use default binning for this kinematical variable.
            bin_limits = label_LaTeX_dict[kinem]["default_bin_limits"]
             
        if x_limits == [0,0]:
            # No x-limits specified, so use default x-limits for this kinematical variable.
            x_limits = label_LaTeX_dict[kinem]["default_x_limits"]
            
        unit = label_LaTeX_dict[kinem]["units"]
            
        x_bins, binw = make_binning_array(bin_limits)

        if len(x_label) == 0:
                # Both muons are being chosen. Change the labels.
            key = "label" if lep_selection_type in ["1","2"] else "independent_label"
            x_label = label_LaTeX_dict[kinem][key]
            if len(unit) > 0:
                x_label += " [{}]".format(unit)
            
        if len(y_label) == 0:
            # User didn't specify label, so make it.
            y_label = hist_y_label(binw, unit)
            
    #         if len(title) == 0:
    #             title = "Selections:\n" + self.cuts
        
        ax, bin_vals, bin_edges, stats = make_1D_dist(ax, data, x_limits, x_bins, 
                                                      x_label, y_label, title, y_max=-1, log_scale=False)
        
        # Nested dictionaries.
        # Initializing this particular kinematic bin's dictionary of stats
        # for this plotted kinematic. 
        self.stats_dict[kinem] = {}
        self.stats_dict[kinem]['hist_stats'] = stats
        self.stats_dict[kinem]['bin_vals'] = bin_vals
        self.stats_dict[kinem]['bin_edges'] = bin_edges

        textbox_text = r"Selection type = {}:".format(lep_selection_type) + "\n"
        if lep_selection_type in ["1","2"]:
            textbox_text = textbox_text.replace("= ", r"$\mu$")
        textbox_text += self.cuts
        ax.text(0.03, 0.87, textbox_text, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes)
        
        if (iter_gaus[0]):
            # Do iterative fitting procedure.
            # Produce a dictionary of stats for these fits.
            best_guess_coeff = max(bin_vals)  # The Gaus coeff is usually around the max height of Gaus curve.
            best_guess_mean = stats[1]
            best_guess_stdev = stats[3]
            stats_dict, ax = iterative_fit_gaus(iter_gaus[1], bin_edges, bin_vals, 
                                                param_guess=[best_guess_coeff, best_guess_mean, best_guess_stdev],
                                                param_bounds=([0,-1000,-100], [999999999,1000,100]),
                                                ax=ax, draw_on_axes=True, verbose=self.verbose)
            # Use plotted kinem as the key for this dict of stats. 
            self.stats_dict[kinem]['fit_stats'] = stats_dict         
            
    def count_in_pT_eta_region_exclusive(self):
        """
        Return the number of muons in the DF subset which pass the eta and pT cuts. 
        'Exclusive' because only those muons which pass such cuts are counted 
        (as opposed to counting both muons from event just because one pass the cuts).
        """
        pass

class KinBinOrganizerOLD():
    """
    FIXME: Created in the spring of 2020, back when I was plotting with plt.
    Creates and organizes KinematicBin objects which all have exactly the same cuts, except for d0 cuts. 
    """
    def __init__(self, 
                 df, 
                 n_evts_scan=1000, 
                 massZ_cut_ls=[60,120],
                 eta_cut_ls=[0.0, 0.3], 
                 pT_cut_ls=[20, 60], 
                 d0_type='BS', 
                 dR_cut=0.05, 
                 use_ptotal_instead=False, 
                 verbose=True):
        """
        These cuts will apply to all KinematicBin objects within this KinBinOrganizer container.
        """
        self.df = df 
        self.n_evts_scan = n_evts_scan
        self.massZ_cut_ls = massZ_cut_ls
        self.eta_cut_ls = eta_cut_ls
        self.pT_cut_ls = pT_cut_ls
        self.d0_type = d0_type
        self.dR_cut = dR_cut
        self.use_ptotal_instead = use_ptotal_instead
        self.verbose = verbose
        
    def make_kbin_ls_over_d0_range(self, d0_bin_limits):
        """
        Make a list of KinematicBin objects, all with identical cuts EXCEPT d0q cuts. 
        """
        d0_bin_arr, d0_bin_width = make_binning_array(d0_bin_limits)
        
        # if d0_bin_width < 0.0005:
        #     err_msg = f"WARNING: d0_bin_width ({d0_bin_width}) is too small (d0_bin_width < 0.0005).\nStopping now."
        #     raise ValueError(err_msg)    
            
        self.d0_bin_arr = d0_bin_arr
        self.d0_bin_arr_shifted = shift_binning_array(d0_bin_arr)
        
        kbin_ls = []
        for elem in range(len(d0_bin_arr)-1):
            # Make a kbin for each d0 bin.
            d0_this = d0_bin_arr[elem]
            d0_next = d0_bin_arr[elem+1]

            kb = KinematicBin(df_original=self.df, 
                              n_evts=self.n_evts_scan, 
                              massZ_cut_ls=self.massZ_cut_ls,
                              eta_cut_ls=self.eta_cut_ls, 
                              pT_cut_ls=self.pT_cut_ls, 
                              d0q_cut_ls=[d0_this, d0_next],
                              d0_type=self.d0_type,
                              dR_cut=self.dR_cut,
                              use_ptotal_instead=self.use_ptotal_instead, 
                              verbose=self.verbose)
            kbin_ls.append(kb)
            
        self.kbin_ls = kbin_ls

    def plot_dpToverpT_for_kbin_ls(self, kinem="delta_pToverpT1", lep_selection_type='independent', 
                                   x_limits=[-0.3, 0.3], bin_limits=[-0.3, 0.3, 0.004], 
                                   run_over_only_n_evts=-1, 
                                   ax=None, y_max=-1, log_scale=False, 
                                   iter_gaus=(False, 3),
                                   make_pdf=False, pdf_obj=None ):
        """
        The main purpose of making these plots is to fill the stats_dict for each kbin.
        
        Parameters
        ----------
        make_pdf : 
        
        """
        for kb in self.kbin_ls:
            kb.plot_1D_kinematics(kinem=kinem, lep_selection_type=lep_selection_type, 
                                  x_limits=x_limits, bin_limits=bin_limits, run_over_only_n_evts=run_over_only_n_evts, 
                                  ax=ax, x_label="", y_label="", title="", y_max=y_max, log_scale=log_scale, 
                                  iter_gaus=iter_gaus)
            if (make_pdf):
                pdf_obj.savefig()
            plt.close()
        # All kbins have filled their stats_dict.
            
    #     def get_dpToverpT_iter_gaus_fit_means(self, kinem):
    def get_iter_gaus_fit_stats(self, kinem):
        """"""
        # Get graph values.
        self.hist_mean_ls     = []
        self.hist_mean_err_ls = []
        self.fit_mean_ls      = []
        self.fit_mean_err_ls  = []
        
        for kb in self.kbin_ls:
    #             self.hist_mean_ls.append(kb.stats_dict['delta_pToverpT1']['hist_stats'][1])
    #             self.hist_mean_err_ls.append(kb.stats_dict['delta_pToverpT1']['hist_stats'][2])
    #             self.fit_mean_ls.append(kb.stats_dict['delta_pToverpT1']['fit_stats']['mean_ls'][-1])
    #             self.fit_mean_err_ls.append(kb.stats_dict['delta_pToverpT1']['fit_stats']['mean_err_ls'][-1])
            self.hist_mean_ls.append(kb.stats_dict[kinem]['hist_stats'][1])
            self.hist_mean_err_ls.append(kb.stats_dict[kinem]['hist_stats'][2])
            self.fit_mean_ls.append(kb.stats_dict[kinem]['fit_stats']['mean_ls'][-1])
            self.fit_mean_err_ls.append(kb.stats_dict[kinem]['fit_stats']['mean_err_ls'][-1])
