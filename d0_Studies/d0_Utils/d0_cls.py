import numpy as np
import matplotlib.pyplot as plt

# Not all of these may be used here. Just saving time for now.
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import (change_cmap_bkg_to_white, save_plots_to_outpath, make_1D_dist, get_stats_1Dhist, 
                                    get_stats_2Dhist, hist_y_label, make_2by2_subplots_for_ratioplots,
                                    make_stats_legend_for_1dhist, make_stats_legend_for_2dhist, 
                                    make_stats_legend_for_gaus_fit)
from PyUtils.Utils_Physics import theta2pseudorap, pseudorap2theta, calc_dR, calc_dphi
from PyUtils.Utils_StatsAndFits import gaussian, fit_with_gaussian, iterative_fit_gaus
from PyUtils.Utils_Collection_Helpers import weave_lists
from d0_Utils.d0_fns import (make_binning_array, centers_of_binning_array, get_subset_mask, 
                             make_kinem_subplot, combine_cut_list, calc_x_err_bins)
from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict

class HistInfo:
    """
    Right now, just the most basic object possible.
    """
    def __init__(self):
        pass
    
class HistProxy:
    """
    This class needs serious salvaging. Has not been touched in a long time. 
    May have useful methods. 
    """
    def __init__(self,
#                  infile,
                 tree,
                 year,
                 charge,
                 bspv,
                 binning_style,  # 'pT' or 'eta',
                 n_bins,
                 x_axis_bounds_list, # 2-element list
                 d0_bin_bounds_list,
                 eta_bin_bounds_list,
                 pT_cuts_internal_list,
                 massZ_cut,
                 test_in_tier2):
#                 make_big_pdf):

#         self.infile    = infile
        self.tree      = tree
        self.year      = year
        self.charge    = charge
        self.bspv      = bspv
        self.bin_style = binning_style
        self.x_bounds = x_axis_bounds_list
        self.d0_bin_bounds = d0_bin_bounds_list  # A 2-element list.
        self.eta_bin_bounds = eta_bin_bounds_list
        self.pT_cuts    = pT_cuts_internal_list  # Was pT_low_str and pT_high_str
        self.massZ_cut = massZ_cut
        self.n_bins    = n_bins
        self.test_in_tier2 = test_in_tier2
        self.make_big_pdf = make_big_pdf
        
        self.ID        = charge_dict[charge]
        
        self.draw_hist()
    
    def draw_hist(self):
        
#         infile    = self.infile
        year      = self.year
        tree      = self.tree
        charge    = self.charge
        bspv      = self.bspv
        bin_style = self.bin_style
        pT_cuts    = self.pT_cuts  # Was pT_low_str and pT_high_str
        massZ_cut = self.massZ_cut
        n_bins    = self.n_bins
        ID        = self.ID
        
        x_min = self.x_bounds[0]
        x_max = self.x_bounds[1]
        this_d0_str = self.d0_bin_bounds[0]
        next_d0_str = self.d0_bin_bounds[1]
        this_eta_str = str(self.eta_bin_bounds[0])
        next_eta_str = str(self.eta_bin_bounds[1])
        pT_low_str = str(self.pT_cuts[0])
        pT_high_str = str(self.pT_cuts[1])
               
        
        c = ROOT.TCanvas()
        c.cd()

        if bin_style in 'eta':
            cuts_per_event  = "(Id1 == %s && %s < d0%s1 && d0%s1 < %s && %s < eta1 && eta1 < %s && %s < pT1 && pT1 < %s) || " % (ID, this_d0_str, bspv, bspv, next_d0_str, this_eta_str, next_eta_str, pT_low_str, pT_high_str)
            cuts_per_event += "(Id2 == %s && %s < d0%s2 && d0%s2 < %s && %s < eta2 && eta2 < %s && %s < pT2 && pT2 < %s)"     % (ID, this_d0_str, bspv, bspv, next_d0_str, this_eta_str, next_eta_str, pT_low_str, pT_high_str)
            hist_name = "deltapT_mu%s_%s__%s_pT_%s__%s_eta_%s__%s_d0_%s" % (charge, bspv, pT_low_str, pT_high_str, this_eta_str, next_eta_str, this_d0_str, next_d0_str)
        # WARNING: 'this_pT_str' is NOT the same as 'pT_low_str'
        # 'this_pT_str' is for binning over pT, whereas 'pT_low_str' is for implicit pT cuts during eta binning.

        #----- NEEDS DEBUGGING AND TESTING -----#
        elif bin_style in 'pT':
            cuts_per_event  = "(Id1 == %s && %s < d0%s1 && d0%s1 < %s && %s < pT1 && pT1 < %s) || " % (ID, this_d0_str, bspv, bspv, next_d0_str, this_pT_str, next_pT_str)
            cuts_per_event += "(Id2 == %s && %s < d0%s2 && d0%s2 < %s && %s < pT2 && pT2 < %s)"     % (ID, this_d0_str, bspv, bspv, next_d0_str, this_pT_str, next_pT_str)
            hist_name = "deltapT_mu%s_%s__%s_pT_%s__%s_d0_%s" % (charge, bspv, this_pT_str, next_pT_str, this_d0_str, next_d0_str) 
#        elif bin_style in 'eta_with_pT_cut'
        #---------------------------------------#

        if len(massZ_cut) != 0:
            cuts_per_event = "("+cuts_per_event+") && %s" % (massZ_cut)

        hist_name = make_str_title_friendly(hist_name)  
        hist_name = hist_name.replace('+','pos') # Must do this because fn() above doesn't work for some reason?
        hist_name = hist_name.replace('-','neg')

        h = ROOT.TH1F(hist_name, cuts_per_event, n_bins, x_min, x_max)
        
        # Fill up the histo with either lep1's deltapT info, if it has the ID of interest.
        # Otherwise use the other lep's deltapT info.
        tricky_root_expr  = "10000*(pT1-genLep_pt1)/genLep_pt1/genLep_pt1*(Id1==%s) + " % ID
        tricky_root_expr += "10000*(pT2-genLep_pt2)/genLep_pt2/genLep_pt2*(Id2==%s)" % ID
        #--- Confirmed that this works by using a Google Spreadsheet! ---#
        
        # @@@@@ WARNING: this is just for 2016 file!!! @@@@@
        tree.Draw("%s >> %s"%(tricky_root_expr, hist_name), cuts_per_event, "")
        
        latex_name  = "10^{4} #times (p_{T}^{RECO}-p_{T}^{GEN})/(p_{T}^{GEN})^{2}"
        bw = bin_width/10000.
        h.GetXaxis().SetTitle(latex_name)
        h.GetYaxis().SetTitle("Events / [%.4f GeV^{-1}]" % bw)
        h.GetXaxis().SetTitleOffset(1.3)
        h.GetYaxis().SetTitleOffset(1.3)
        h.SetTitle(cuts_per_event)
        # h1.SetAxisRange(0.0, 0.1, "X")   
        # h1.SetLabelSize(0.03, "Y")                        
        # h1.SetLineColor(1)
        h.Draw("hist 9 same")
        # h1.Draw("e1 hist 9 same")
        c.Draw()
        
        mean = h.GetMean()
        mean_err = h.GetMeanError()
        stdev = h.GetStdDev()
        stdev_err = h.GetStdDevError()
        
        self.c = c
        self.h = h
        self.tree = tree
        self.hist_name = hist_name
        self.cuts_per_event = cuts_per_event
        self.h_mean = mean
        self.h_mean_err = mean_err
        self.h_stdev = stdev
        self.h_stdev_err = stdev_err
        
        # Save plots.
        global outpath_plots_deltapT_dist
        global d0_min 
        global d0_max
        global make_plots_deltapT_vs_d0
        
        if (self.test_in_tier2):
            fullpath = os.path.join(outpath_plots_deltapT_dist, hist_name)
            c.SaveAs(fullpath + '.png')
            c.SaveAs(fullpath + '.pdf')  
            
        #----- I will probably have to make_big_pdf outside this class.
#        elif (self.make_big_pdf):
#            
#            if bin_style in 'eta':
#                pdf_title = ("deltapT_dist_MC%s_mu%s__%s_eta_%s__%s_d0%s_%s_increm%s" % (year, charge, this_eta_str, next_eta_str, str(d0_min), str(bspv), str(d0_max), str(d0_bin_width)))
#                pdf_title = ("deltapT_dist_MC%s_mu%s_%s__%s_pT_%s__%s_eta_%s__%s_d0%s_%s_increm%s" % (year, charge, bspv, pT_low_str, pT_high_str, this_eta_str, next_eta_str, str(d0_min), str(bspv), str(d0_max), str(d0_bin_width)))
#            elif bin_style in 'pT':
#                if len(this_pT_str) < 4:  # Turn: 5p0 --> 05p0, for plot-ordering purposes.
#                    this_pT_str = '0'+this_pT_str  
#                pdf_title = ("deltapT_dist_MC%s_mu%s_%s__%s_pT_%s__%s_d0%s_%s_increm%s" % (year, charge, bspv, this_pT_str, next_pT_str, str(d0_min), str(bspv), str(d0_max), str(d0_bin_width)))
#            
#            pdf_title = make_str_title_friendly(pdf_title)
#            fullpath = os.path.join(outpath_plots_deltapT_dist, pdf_title)
#            
#            if this_d0_str in str(d0_min):
#                c.Print(fullpath + '.pdf[')
#            c.Print(fullpath + '.pdf')
#            if next_d0_str in str(d0_max):
#                c.Print(fullpath + '.pdf]')

        # Extract mean, RMS, and store for later.
        #--- FIXME! Implement a fit with CBxBW + exp
        #--- FIXME! Implement a fit with Voigtian
        
            
        if (verbose):
#             print("deltapT_mean_list:", deltapT_mean_list)
#             print("deltapT_mean_err_list:", deltapT_mean_err_list)
            if bin_style in 'eta':
                print("Completed %s<pT<%s, %s<eta<%s, %s<d0%s<%s \n" % (pT_low_str, pT_high_str, this_eta_str, next_eta_str, this_d0_str, bspv, next_d0_str))
                print("Here's the important stored info:")
                print("self.c", self.c)
                print("self.h", self.h)
                print("self.tree", self.tree)
                print("self.hist_name", self.hist_name)
                print("self.cuts_per_event", self.cuts_per_event)
                print("self.h_mean", self.h_mean)
                print("self.h_mean_err", self.h_mean_err)
                print("self.h_stdev", self.h_stdev)
                print("self.h_stdev_err", self.h_stdev_err)
            if bin_style in 'pT':
                print("Completed %s<pT<%s, %s<d0%s<%s \n" % (this_pT_str, next_pT_str, this_d0_str, bspv, next_d0_str))
                
        if (make_plots_deltapT_vs_d0):
            return mean, mean_err
        
class KinematicBin():

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
    
    
#--------------------------------#

class KinBinOrganizer():
    """
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
        
        if d0_bin_width < 0.0005:
            err_msg = f"WARNING: d0_bin_width ({d0_bin_width}) is too small (d0_bin_width < 0.0005).\nStopping now."
            raise ValueError(err_msg)    
            
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
            
class GraphLine():
    """
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
        low_x_err, high_x_err = calc_x_err_bins(self.x_vals)
        
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

class KinBin3D():
    def __init__(self, eta_range, pT_range, qd0_range, n_entries, kinem, fit_stats_dict, pT_stats_ls, qd0_stats_ls):
        self.eta_range = eta_range
        self.pT_range = pT_range
        self.qd0_range = qd0_range
        self.n_entries = n_entries
        self.kinem = kinem
        self.fit_stats_dict = fit_stats_dict
        self.pT_stats_ls = pT_stats_ls,
        self.qd0_stats_ls = qd0_stats_ls

    # def get_best_fit_results():
    #     pass
        
    # def 