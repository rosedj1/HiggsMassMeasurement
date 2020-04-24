import numpy as np
import matplotlib.pyplot as plt

# Not all of these may be used here. Just saving time for now.
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Plotting import (change_cmap_bkg_to_white, save_plots_to_outpath, make_1D_dist, get_stats_1Dhist, 
                                    get_stats_2Dhist, hist_y_label, make_2by2_subplots_for_ratioplots,
                                    add_underoverflow_entries, make_stats_legend_for_1dhist, make_stats_legend_for_2dhist, 
                                    make_stats_legend_for_gaus_fit)
from PyUtils.Utils_Physics import theta2pseudorap, pseudorap2theta, calc_dR, calc_dphi
from PyUtils.Utils_StatsAndFits import gaussian, fit_with_gaussian, iterative_fit_gaus
from PyUtils.Utils_Collection_Helpers import weave_lists
from d0_Utils.d0_fns import (make_binning_array, shift_binning_array, get_subset_mask, make_kinem_subplot)
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
        
        
        
        
class KinemBinnedEtaPt():

    def __init__(self, df_original, n_evts, eta_cut_ls, pT_cut_ls, use_ptotal_instead=False, dR_cut=0.02, verbose=False):
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
        
        self.n_evts_asked_for = n_evts
        self.n_evts_found = -999
        self.kinem_vals_after_selection = {}
        
        self.use_ptotal_instead = use_ptotal_instead
        self.p_str = "p" if (self.use_ptotal_instead) else "pT"
        self.p_str_latex = r"$p^{\mathrm{REC}}$" if (self.use_ptotal_instead) else r"$p_{T}^{\mathrm{REC}}$"
        
        self.eta_min   = eta_cut_ls[0]
        self.eta_max   = eta_cut_ls[1]
        self.pT_min    = pT_cut_ls[0] 
        self.pT_max    = pT_cut_ls[1]         
        self.dR_cut    = dR_cut    
        self.massZ_min = 60
        self.massZ_max = 120    
        
        self.apply_pT_eta_cuts(df, verbose)
        self.save_exclusive_masks()
        
    def apply_pT_eta_cuts(self, df, verbose=False):
        """
        Creates a subset of the original DataFrame.
        The subset contains only the events in which either muon1 passes 
        all eta and pT cuts, or muon2 does, or both muons do. 
        """
        # Cuts:
        eta_min   = self.eta_min  
        eta_max   = self.eta_max  
        pT_min    = self.pT_min   
        pT_max    = self.pT_max   
        dR_cut    = self.dR_cut   
        massZ_min = self.massZ_min
        massZ_max = self.massZ_max

        #----------------#
        #--- Analysis ---#
        #----------------#
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
        
        df.loc[:,'d0BSxq1'] = df['d0BS1'] * df['Id1'].replace(13,-1).replace(-13,1)
        df.loc[:,'d0BSxq2'] = df['d0BS2'] * df['Id2'].replace(13,-1).replace(-13,1)
        df.loc[:,'d0PVxq1'] = df['d0PV1'] * df['Id1'].replace(13,-1).replace(-13,1)
        df.loc[:,'d0PVxq2'] = df['d0PV2'] * df['Id2'].replace(13,-1).replace(-13,1)
        
        # Create masks.
        mask_dR = (dR1_ser < self.dR_cut) & (dR2_ser < self.dR_cut)

        mask_massZ = (massZ_min < df['massZ']) & (df['massZ'] < massZ_max)
        
        if (self.use_ptotal_instead):
            mask_pT1 = (pT_min < df['p1']) & (df['p1'] < pT_max) 
            mask_pT2 = (pT_min < df['p2']) & (df['p2'] < pT_max)
        else:
            mask_pT1 = (pT_min < df['pT1']) & (df['pT1'] < pT_max) 
            mask_pT2 = (pT_min < df['pT2']) & (df['pT2'] < pT_max)
        
        mask_eta1 = (eta_min < abs(df['eta1'])) & (abs(df['eta1']) < eta_max)
        mask_eta2 = (eta_min < abs(df['eta2'])) & (abs(df['eta2']) < eta_max)        

        # Combine masks.
        mask_kinembin_lep1 = mask_eta1 & mask_pT1
        mask_kinembin_lep2 = mask_eta2 & mask_pT2

        all_masks = mask_dR & mask_massZ & (mask_kinembin_lep1 | mask_kinembin_lep2)
        
        # Apply masks and update DataFrame.
        self.binned_df = df[all_masks]
        
        self.cuts =  r"$%.1f < m_{\mu\mu} < %.1f$ GeV" % (self.massZ_min, self.massZ_max)
        self.cuts += r",   $%.2f < \left| \eta \right| < %.2f$" % (self.eta_min, self.eta_max)
        self.cuts += r",   $%d <$ %s $< %d$ GeV" % (self.pT_min, self.p_str_latex, self.pT_max)  # The string brings in its own '$'.
        self.cuts += r",   $\Delta R < %.3f$" % (self.dR_cut)
        self.n_evts_found = len(self.binned_df)
        
        if (verbose): 
            perc = self.n_evts_found / float(self.n_evts_asked_for) * 100.
            print(r"Events found: {} ({:.2f}% of total events), using cuts: {}".format(self.n_evts_found, perc, self.cuts))
            
    def save_exclusive_masks(self):
        """
        Save two boolean masks:
            mask 1 -> shows which events muon1 passes
            mask 2 -> shows which events muon2 passes
        Keeping the masks separate like this is what is "exlusive" about the masks. 
        """
        df = self.binned_df
    
        # Create Masks
        if (self.use_ptotal_instead):
            mask_pT1 = (self.pT_min < df['p1']) & (df['p1'] < self.pT_max) 
            mask_pT2 = (self.pT_min < df['p2']) & (df['p2'] < self.pT_max)
        else:
            mask_pT1 = (self.pT_min < df['pT1']) & (df['pT1'] < self.pT_max) 
            mask_pT2 = (self.pT_min < df['pT2']) & (df['pT2'] < self.pT_max)
        
        mask_eta1 = (self.eta_min < abs(df['eta1'])) & (abs(df['eta1']) < self.eta_max)
        mask_eta2 = (self.eta_min < abs(df['eta2'])) & (abs(df['eta2']) < self.eta_max)        

        # Combine masks and save them.
        self.mask_kinembin_lep1 = mask_eta1 & mask_pT1
        self.mask_kinembin_lep2 = mask_eta2 & mask_pT2
        
    def make_2D_plot(self, 
                     x_kinem, y_kinem, 
                     x_bin_limits=[0, 1, 0.1], y_bin_limits=[0, 1, 0.1],
                     run_over_only_n_evts=-1, 
                     title="",
                     exclusive=True,
                     save_plot=False, save_as_png=False, verbose=False, outpath="",
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
            The kinematical variable to be plotted along x-axis. 
            Only works for kinematics which end with '1' or '2'.
            - Example: x_kinem="delta_theta" (for which there are two branches: "delta_theta1", "delta_theta2")
        y_kinem : str
            The kinematical variable to be plotted along y-axis. 
            Only works for kinematics which end with '1' or '2'.
            - Example: y_kinem="delta_eta" (for which there are two branches: "delta_eta1", "delta_eta2")
        x_bin_limits : list or array-like of floats
            The bin limits on the horizontal axis. [bin_min_left_edge, bin_max_right_edge, bin_width]
            - Example: [-2.5, 2.5, 0.1]
        y_bin_limits : list or array-like of floats
            The bin limits on the vertical axis. [bin_min_left_edge, bin_max_right_edge, bin_width]
            - Example: [-2.5, 2.5, 0.1]
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
        verbose : bool
            Get debug and code progress info.
        outpath : str
            Path to save plot.
            
        FIXME:
        - Maybe generalize this into a method to make any 2D plot.
        """           
        if run_over_only_n_evts == -1:
            run_over_only_n_evts = self.n_evts_found       

        n = run_over_only_n_evts
        
        x_kinem1 = x_kinem + "1"
        x_kinem2 = x_kinem + "2"
        y_kinem1 = y_kinem + "1"
        y_kinem2 = y_kinem + "2"
        
        x_vals = self.apply_mask_get_data(x_kinem1, x_kinem2, run_over_only_n_evts=n, exclusive=exclusive)    
        y_vals = self.apply_mask_get_data(y_kinem1, y_kinem2, run_over_only_n_evts=n, exclusive=exclusive)    
        
        # A special case to make comparison of (delta_phi vs. delta_theta) easy with (delta_phi vs. delta_eta).
        if (x_kinem1[:-1] == "delta_theta") and (y_kinem1[:-1] == "delta_phi"):
            x_vals *= -1

        #--- Make plots ---#
        if (ax is None):
            f, ax = plt.subplots(figsize=(12.8, 9.6))
        
        x_2D_bins, x_2D_bin_width = make_binning_array(x_bin_limits)
        y_2D_bins, y_2D_bin_width = make_binning_array(y_bin_limits) 
        
        # Plot 1: dphi vs. deta
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
        cuts = "Selection:\n" + r"{}".format(self.cuts)
        ax.set_title(title + cuts)#, fontsize=label_size)
    
        # Stats: 
        stat_text_x = 0.82
        stat_text_y = 0.9
        
        stats_ls = get_stats_2Dhist(x_vals, y_vals)
        leg_label = make_stats_legend_for_2dhist(stats_ls)
        ax.text(stat_text_x, stat_text_y, leg_label, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
        
        newcmp = change_cmap_bkg_to_white('rainbow')
        bin_vals, x_bin_edges, y_bin_edges, im = ax.hist2d(x_vals, y_vals, bins=[x_2D_bins, y_2D_bins], cmap=newcmp)
        plt.colorbar(im, ax=ax)
        
        # Save plots, if you want. 
        file_name  = "2Dplot_dphi_vs_detaANDdtheta"
        file_name += "__%.2f_eta_%.2f" % (self.eta_min, self.eta_max)
        file_name += "__%d_%s_%d" % (self.pT_min, self.p_str, self.pT_max)
        
        save_plots_to_outpath(save_plot, outpath, file_name, save_as_png, verbose)


    def apply_mask_get_data(self, kinem1, kinem2, run_over_only_n_evts=-1, exclusive=True):
        """
        Return the kinematic muon data which pass selections of this kinematic bin.
        Example: Retrieve all "delta_R1" and "delta_R2" values after masks have been applied.

        Parameters
        ----------
        kinem1 : str
            A branch in the DataFrame or root file. Should correspond to muon1. 
                Examples: "pT1", "genLep_pt1"
        kinem2 : str
            A branch in the DataFrame or root file. Should correspond to muon2. 
                Examples: "pT2", "genLep_pt2"
        run_over_only_n_evts : int
            Number of events to run over. Use '-1' to specify all available events.
        exclusive : bool
            Means "only put muons which passed all selections in this plot".
            FIXME: It must be set to True for now...
        
        Returns
        -------
        kinem_vals : array 
            kinem1 values (muon1) woven together with kinem2 values (muon2).
            All values satisfy selection criteria for this kinematic bin.
        """
        # Quick check.
        if kinem1[:-1] != kinem2[:-1]:
            raise ValueError("Problem! `kinem1` and `kinem2` are not the same kind of kinematical variable (e.g. `pT1` and `pT2`).\nStopping now.")
        
        df = self.binned_df
        #--- Get data ---#
        if (exclusive):   
            # Keep only the muons which pass the kinem bin selection. 
            vals_muon1 = df[kinem1][self.mask_kinembin_lep1].values  # Muon 1 passes kinem bin selection.
            vals_muon2 = df[kinem2][self.mask_kinembin_lep2].values  # Muon 2 passes kinem bin selection.
                        
            if run_over_only_n_evts == -1:
                # Running over all events. No need to slice.
                kinem_vals = np.append(vals_muon1, vals_muon2)
            else:
                n = run_over_only_n_evts
                # Weave muon1 and muon2 values together in systematic way.
                # Otherwise muon1 values would mostly be selected when doing a slice, like [:n].
                weave_ls = weave_lists(vals_muon1, vals_muon2)
                kinem_vals = np.array(weave_ls)[:n]

        else:
            #--------# DEPRECATED FOR NOW.
            # Keep both muons from each event in which AT LEAST ONE muon passed the kinem bin selection. 
            raise RuntimeError("Stopping now. This section hasn't been fully developed.\nSet exclusive=True.")
            
            kinem_vals = np.append(self.binned_df['delta_eta1'][:n].values, self.binned_df['delta_eta2'][:n].values)
            y_vals = np.append(self.binned_df['delta_phi1'][:n].values, self.binned_df['delta_phi2'][:n].values)

            x2_vals = np.append(self.binned_df['delta_theta1'][:n].values, self.binned_df['delta_theta2'][:n].values)
            x2_vals = x2_vals * -1
            y_vals = np.append(self.binned_df['delta_phi1'][:n].values, self.binned_df['delta_phi2'][:n].values)
            #--------#
        
        combined_key = "{} and {}".format(kinem1, kinem2)
        self.kinem_vals_after_selection[combined_key] = kinem_vals
        
        return kinem_vals
    

    def plot_kinem_genrec_comparison(self,
                              kinem_gen, kinem_rec, 
                              x_range_ls=[-1, 1],
                              bin_limits=[-0.5,0.5,0.5], 
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
            x_range_ls : 2-element list
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

            lep = kinem_gen[-1]  # Get last character of kinematic. Example: pT1 --> 1. Should be a string!
            
            #--- Get data ---#
            if ("1" in kinem_gen) and ("1" in kinem_rec): 
                data_rec = df[kinem_rec][self.mask_kinembin_lep1].values  # Muon 1 passes kinem bin selection.
                data_gen = df[kinem_gen][self.mask_kinembin_lep1].values  # Muon 1 passes kinem bin selection.
            elif ("2" in kinem_gen) and ("2" in kinem_rec):
                data_rec = df[kinem_rec][self.mask_kinembin_lep2].values  # Muon 1 passes kinem bin selection.
                data_gen = df[kinem_gen][self.mask_kinembin_lep2].values  # Muon 1 passes kinem bin selection.
            else:
                err_msg = "\n    Either kinem_gen or kinem_rec does not end with a '1' or '2', or they are not the same as each other.\nStopping now."
                raise ValueError(err_msg)

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
            ax.set_title(r"Selection: {}".format(self.cuts), fontsize=textsize_title)
            
            plt.minorticks_on()

            x_min = x_range_ls[0]
            x_max = x_range_ls[1]
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

    def plot_1D_kinematics(self, lep, kinem, x_limits=[0,0], bin_limits=[0,0,0], ax=None, x_label="", y_label="", title="", y_max=-1, log_scale=False, iter_gaus=(False, 0)):
        """
        Make a histogram of a kinematical variable in the DF.  
        FIXME: 
            - Currently only works with kinem variables that end with `1` or `2`. 
            [ ] Eventually plot other quantities like massZ, etc.
        
        Parameters
        ----------
        lep : int
            Either `1` or `2`. Indicates which lepton you are referring to. 
        kinem : str
            The full name of the kinematical branch/column. 
            E.g. `delta_eta1` or `delta_R2`
        x_limits : 2-element list
            The x_min and x_max to show along x_axis, for viewing purposes: [x_min, x_max]
        bin_limits : 3-element list
            [first_bin_left_edge, last_bin_right_edge, bin_width]   
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
        df = self.binned_df
        
        if ax is None:
            # Axes doesn't exist yet. Make it.
            fig, ax = plt.subplots(figsize=(18,13.5))
            
        if bin_limits == [0,0,0]:
            # No bin limits specified, so use default binning for this kinematical variable.
            bin_limits = label_LaTeX_dict[kinem]["default_bin_limits"]
             
        if x_limits == [0,0]:
            # No x-limits specified, so use default x-limits for this kinematical variable.
            x_limits = label_LaTeX_dict[kinem]["default_x_limits"]
            
        unit = label_LaTeX_dict[kinem]["units"]
        
        #--- Get data ---#
        if lep == 1: 
            data = df[kinem][self.mask_kinembin_lep1].values  # Muon 1 passes kinem bin selection.
        elif lep == 2:
            data = df[kinem][self.mask_kinembin_lep2].values  # Muon 1 passes kinem bin selection.
        else:
            err_msg = "\n    You  must specify lep = 1 or 2.\nStopping now."
            raise ValueError(err_msg)
            
        x_bins, binw = make_binning_array(bin_limits)

        if len(x_label) == 0:
            x_label = label_LaTeX_dict[kinem]["label"]
            if len(unit) > 0:
                x_label += " [{}]".format(unit)
            
        if len(y_label) == 0:
            # User didn't specify label, so make it.
            y_label = hist_y_label(binw, unit)
            
        if len(title) == 0:
#             title = "Selection: " + self.cuts
            title = r"Selection: {}".format(self.cuts)
        
        ax, bin_vals, bin_edges, stats = make_1D_dist(ax, data, x_limits, x_bins, x_label, y_label, title, y_max=-1, log_scale=False)    
        
        if (iter_gaus[0]):
            # Do iterative fitting procedure.
            stats_dict, ax = iterative_fit_gaus(iter_gaus[1], bin_edges, bin_vals, first_mean=stats[1], first_stdev=stats[3], ax=ax)
            self.stats_dict = stats_dict
         
#             iterations = iter_gaus[1]
#             msg = "Performing {} iterative Gaussian fits".format(iterations)
#             if (iterations == 1):
#                 msg = msg.replace("fits", "fit") 
#             print(msg)

#             bin_centers = shift_binning_array(bin_edges)
            
#             count = 0
#             popt = np.zeros(3)
#             while count < iterations:
#                 count += 1
                
# #                 fit_ls = iterative_fit_gaus(bin_centers, bin_vals, fit_range_start=[0,0], iterations=1)
#                 if count == 1:
#                     # First fit: use original histogram's mean and stdev to choose a fit range.
#                     this_mean  = stats[1]
#                     this_stdev = stats[3]
#                 else:
#                     # Otherwise use the last fit's optimized parameters.
#                     this_mean  = popt[1]
#                     this_stdev = popt[2]
#                 this_x_min = this_mean - 2*this_stdev
#                 this_x_max = this_mean + 2*this_stdev
                
#                 mask = get_subset_mask(bin_centers, x_min=this_x_min, x_max=this_x_max)
                
#                 new_bin_centers = bin_centers[mask]
#                 new_bin_vals = bin_vals[mask]
                
#                 popt, popt_err, pcov = fit_with_gaussian(new_bin_centers, new_bin_vals, guess_params=[1,this_mean,this_stdev])
                
#                 # Get the y-vals of the Gaussian fit for plotting
#                 gaus_y_vals = gaussian(new_bin_centers, *popt)

#                 leg_label_fit = make_stats_legend_for_gaus_fit(popt, popt_err)
#                 leg_label_fit = leg_label_fit.replace("Fit", "Fit {}:".format(count))

#                 ax.plot(new_bin_centers, gaus_y_vals, color=color_dict[count+1], label=leg_label_fit, linestyle='-', marker="")
#                 ax.legend()
            
            
    def count_in_pT_eta_region_exclusive(self):
        """
        Return the number of muons in the DF subset which pass the eta and pT cuts. 
        'Exclusive' because only those muons which pass such cuts are counted 
        (as opposed to counting both muons from event just because one pass the cuts).
        """
        pass
#         return n_passed = 