# This file is essentially all old classes from the d0_Analyzer.ipynb. 
# Perhaps some of the class structure is salvageable. 
# It may not even be worth your time. 
# The last time I touched this file was well before 2020-02-20

class HistProxy:
    
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
                 test_in_tier2,
                 make_big_pdf):

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
        this_d0_str = str(self.d0_bin_bounds[0])
        next_d0_str = str(self.d0_bin_bounds[1])
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
        elif bin_style in 'eta_with_pT_cut'
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
        elif (self.make_big_pdf):
            
            if bin_style in 'eta':
                  pdf_title = ("deltapT_dist_MC%s_mu%s__%s_eta_%s__%s_d0%s_%s_increm%s" % (year, charge, this_eta_str, next_eta_str, str(d0_min), str(bspv), str(d0_max), str(d0_bin_width)))
                pdf_title = ("deltapT_dist_MC%s_mu%s_%s__%s_pT_%s__%s_eta_%s__%s_d0%s_%s_increm%s" % (year, charge, bspv, pT_low_str, pT_high_str, this_eta_str, next_eta_str, str(d0_min), str(bspv), str(d0_max), str(d0_bin_width)))
            elif bin_style in 'pT':
                if len(this_pT_str) < 4:  # Turn: 5p0 --> 05p0, for plot-ordering purposes.
                    this_pT_str = '0'+this_pT_str  
                pdf_title = ("deltapT_dist_MC%s_mu%s_%s__%s_pT_%s__%s_d0%s_%s_increm%s" % (year, charge, bspv, this_pT_str, next_pT_str, str(d0_min), str(bspv), str(d0_max), str(d0_bin_width)))
            
            pdf_title = make_str_title_friendly(pdf_title)
            fullpath = os.path.join(outpath_plots_deltapT_dist, pdf_title)
            
            if this_d0_str in str(d0_min):
                c.Print(fullpath + '.pdf[')
            c.Print(fullpath + '.pdf')
            if next_d0_str in str(d0_max):
                c.Print(fullpath + '.pdf]')

        # Extract mean, RMS, and store for later.
        #--- FIXME! Implement a fit with CBxBW + exp
        #--- FIXME! Implement a fit with Voigtian
        
            
        if (verbose):
#             print "deltapT_mean_list:", deltapT_mean_list
#             print "deltapT_mean_err_list:", deltapT_mean_err_list
            if bin_style in 'eta':
                print "Completed %s<pT<%s, %s<eta<%s, %s<d0%s<%s \n" % (pT_low_str, pT_high_str, this_eta_str, next_eta_str, this_d0_str, bspv, next_d0_str)
                print "Here's the important stored info:"
                print "self.c", self.c
                print "self.h", self.h
                print "self.tree", self.tree
                print "self.hist_name", self.hist_name
                print "self.cuts_per_event", self.cuts_per_event
                print "self.h_mean", self.h_mean
                print "self.h_mean_err", self.h_mean_err
                print "self.h_stdev", self.h_stdev
                print "self.h_stdev_err", self.h_stdev_err
            if bin_style in 'pT':
                print "Completed %s<pT<%s, %s<d0%s<%s \n" % (this_pT_str, next_pT_str, this_d0_str, bspv, next_d0_str)
                
        if (make_plots_deltapT_vs_d0):
            return mean, mean_err


class HistInfo:
    
    def __init__(self):
        pass
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
                 test_in_tier2,
                 make_big_pdf):

        self.year      = year
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
      
#         self.draw_hist()
  
    def draw_hist(self):
      
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
        this_d0_str = str(self.d0_bin_bounds[0])
        next_d0_str = str(self.d0_bin_bounds[1])
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
#         elif bin_style in 'eta_with_pT_cut'
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
      
        # Extract mean, RMS, and store for later.
        #--- FIXME! Implement a fit with CBxBW + exp
        #--- FIXME! Implement a fit with Voigtian
      
          
        if (verbose):
#             print "deltapT_mean_list:", deltapT_mean_list
#             print "deltapT_mean_err_list:", deltapT_mean_err_list
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
