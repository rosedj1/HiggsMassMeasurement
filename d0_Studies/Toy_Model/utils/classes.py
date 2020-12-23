# FIXME: graph_dXOverX() may have an issue when called multiple times
# while do_fit = True.

import numpy as np
import ROOT as r
from array import array



class MuonTrack:
    """
    This class represents a muon trajectory.
    """

    def __init__(self, pT, hits_ls, yerr_ls, shift_yvals=(True, 0.01)):
        """Plot and fit the muon track.
        
        hits_ls : list
            The x-vals of hits in Tracker layers (cm).
        yerr_ls : list
            The uncertainties on the y-vals (um).
        shift_yvals : 2-tuple (bool, float)
            The random Gaussian shift in y-vals of each hit.
            The float is the sigma of the Gaussian (cm)
        """
        self.pT_true = pT
        self.hits_ls = hits_ls       # x-vals of Tracker/Pixel hits.
        self.n_hits = len(hits_ls)
        self.xerr_ls = None
        self.yerr_ls = np.array(yerr_ls) / 10000. # Uncertainty on y-vals. Convert to cm.
        self.smear = shift_yvals[0]  # 
        self.sigma = shift_yvals[1]  # sigma from Gaussian used to shift y-val.

        # Since we know true pT, we know what true track should look like.
        self.a_true = self.convert_pT_to_a(pT)
        self.y_vals = self.find_parabolic_y_vals(self.a_true, hits_ls)
        self.y_vals_smear = self.smear_hits(sigma=self.sigma)

        # Beam spot info.
        self.constrain_to_BS = False
        
        self.a_fit = None
        self.a_fit_err = None
        self.d0 = None      # Same as p0 from fit.
        self.d0_err = None 

        self.a_fit_BS = None
        self.a_fit_BS_err = None
        self.d0_BS = None      # Same as p0 from fit.
        self.d0_BS_err = None 
        self.pT_fit = None
        
        self.dpTOverpT = None
        self.daOvera = None
        
    def convert_pT_to_a(self, pT, q=-1, B=3.8):
        """Calculates the 'p2' parameter for a parabola (also called 'a') based on pT: 
            pT = q*B*R
            R = radius of curvature = 1/(2*a)
            
            ==> a = |q|*B / (2*pT)
            
        q : charge of the particle in units of 'e' (example: -1 e, +1 e, )
            # Recall that in HEP, the conversion is:
            # 1 e [unitless] = 1.602E-19 = (0.3 / |c|) * GeV/J
            # Here, |c| is simply the number 3E8. 
        B : magnetic field of particle 
        
        Returns `a` in cm^(-1).
        """
        # Factor of 100 to go from m^(-1) to cm^(-1).
        return 0.3 * abs(q) * B / float(2.0 * pT * 100.)
    
    def convert_a_to_pT(self, a, q=-1, B=3.8):
        """
        Converts `a` in pT, doing the reverse of convert_pT_to_a().
        Assumes a is in cm^(-1).
        Returns pT in GeV/c.
        """
        return 0.3 * abs(q) * B / float(2.0 * a * 100.)
    
    def calc_daOvera(self):
        """Return (a_fit - a_true) / a_true."""
        return float(self.a_fit - self.a_true) / self.a_true
    
    def calc_dpTOverpT(self):
        """Return (pT_fit - pT_true) / pT_true."""
        return float(self.pT_fit - self.pT_true) / self.pT_true
        
    def find_parabolic_y_vals(self, a, x_arr):
        """Return the y-values, assuming a parabolic trajectory."""
        x_arr = np.array(x_arr)
        return a * x_arr**2
    
    def add_BS(self, pos=(0,0), bs_unc_y=10):
        """Add the BS info to start of list of hits, y-vals, and y-errors.
        
        pos : 2-tuple
            Location of BS. 
            (x-val, y-val) in cm
        bs_unc_y : float
            The y-uncertainty on the BS in um.
        """
        self.constrain_to_BS = True
        bs_xval, bs_yval = pos[0], pos[1]
        self.hits_ls = [bs_xval] + self.hits_ls
        self.y_vals = [bs_yval] + self.y_vals
        self.yerr_ls = [bs_unc_y] + self.yerr_ls

    def plot_trajectory(self, leg=None, color=1):
        """
        Draw a TGraph with the muon track as it hits the Pixel/Strip layers.
        Return (TGraph, TLegend) tuple.

        Parameters
        ----------
        sigma : float
            The uncertainty (cm) on all y-vals.
        """
        # Add to the legend.
        if leg is None:
            leg = r.TLegend(0.15, 0.50, 0.45, 0.65)
            # leg = r.TLegend()
        leg_text = r"p_{T}^{gen} = %.0f GeV" % self.pT_true
            
        n_pts = len(self.hits_ls)
        
        # Use y = a*x^2 to model the trajectory of the muon.
        # For high pT this is reasonable. 
        # Radius of curvature goes like: R = 1 / (2*a)
        x_vals = np.array(self.hits_ls)
        y_vals = self.y_vals
        # y_vals = self.find_parabolic_y_vals(self.a_true, x_vals)
        # self.y_vals = y_vals.copy()  # self.y_vals will get smear later. 
        
        # Then convert them to array.arrays afterward.
        x_arr = array('f', x_vals)
        
        graph_title = r" Muon trajectory for %s" % leg_text
        if self.smear:
            toppiece = graph_title
            botpiece = r"y-vals shifted by G(#mu=0, #sigma=%.4g #mum)" % (self.sigma * 1E4)
            graph_title  = r"#splitline{%s}{%s}" % (toppiece, botpiece)
            y_arr = array('f', self.y_vals_smear)
        else:
            y_arr = array('f', y_vals)

        # Make the graph.
        xerr_arr = array('f', np.zeros_like(self.yerr_ls)) if self.xerr_ls is None else array('f', self.xerr_ls)
        yerr_arr = array('f', self.yerr_ls)
        gr = r.TGraphErrors(n_pts, x_arr, y_arr, xerr_arr, yerr_arr)
        # Clean up the title a bit.
        toppiece = graph_title
        str_err_ls = ', '.join([r"%.0f"%x for x in yerr_arr])
        botpiece = r"#sigma_{y} = [%s] #mum" % str_err_ls
        graph_title  = r"#splitline{%s}{%s}" % (toppiece, botpiece)
        # Pretty up the graph.
        gr.SetLineColor(color)
        gr.SetLineWidth(1)
        gr.SetMarkerColor(color)
        gr.SetMarkerStyle( 21 )
        gr.SetTitle(graph_title)
        gr.GetXaxis().SetTitle('Transverse Pixel/Strip Positions (x) (cm)')
        gr.GetYaxis().SetTitle('y (cm)')
        gr.GetYaxis().SetTitleOffset(1.5)
        return (gr, leg)

    def smear_hits(self, sigma=0.01):
        """
        Return the y-value of each hit in vals after adding a small Gaussian 
        uncertainty. 
        
        sigma should be given in cm.
        By default, sigma of Gaus is 100 um = 0.01 cm (since 
        positions are measured in cm).
        """
        smears = np.random.normal(loc=0, scale=sigma, size=self.n_hits)
        return self.y_vals + smears

    def fit_hits_pol2(self, graph, color=2):
        """Return the pol2 fit parameters for a set of (x,y) coordinates.
        
        Parameters
        ----------
        graph : TGraph-like
            Graph to be fit.
            
        Returns
        -------
        fit_func : TF1
            The fitting function whose best-fit parameters are obtained from `graph`.
        """
        x_min = min(self.hits_ls) 
        x_max = max(self.hits_ls)
        fit_func = r.TF1(f'f1_{graph.GetName()}', '[0]+[1]*x+[2]*x^2', x_min, x_max)
        fit_func.SetLineColor(color)
        fit_func.SetLineWidth(2)
        fit_func.SetLineStyle(2)
        # Fit it onto a histogram `h1`:
        result = graph.Fit(fit_func,'S')
        # The option 'S' saves the fit results into a pointer.
        r.gStyle.SetOptFit(111)
        
        self.set_fit_vals(fit_func)
        
        return (fit_func, result)
    
    def set_fit_vals(self, fit_func):
        """Check that the best-fit params are reasonable and set attributes."""
        a_fit = fit_func.GetParameter(2)
        a_fit_err = fit_func.GetParError(2)
        d0 = fit_func.GetParameter(0)
        d0_err = fit_func.GetParError(0)
        assert a_fit is not None
        assert a_fit_err is not None
        
        self.pT_fit = self.convert_a_to_pT(self.a_fit)
        self.dpTOverpT = self.calc_dpTOverpT()
        self.daOvera = self.calc_daOvera()

        # May need to set: self.a_fit, self.a_fit_err, self.d0

        return ()



class BiasPlotter:
    def __init__(self, pT, hits_ls, smear=True):
        """Toss n_toys to test the theoretical bias of da/a or dpT/pT vs. d0."""
        self.pT = pT
        self.hits_ls = hits_ls
        self.smear = smear
        
        # Attribute which will eventually be filled.
        self.n_toys = None
        self.d0_ls = None
        self.daOvera_ls = None
        self.dpTOverpT_ls = None

    def toss_toys(self, n_toys):
        """Make n_toys worth of muon trajectories with pT (self.pT)."""
        self.n_toys = n_toys
        d0_ls = []
        daOvera_ls = []
        dpTOverpT_ls = []
        
        for n in range(n_toys):
            hitplotorg = HitPlotOrg(pT=self.pT, hits_ls=self.hits_ls)
            gr_smear, leg_smear = hitplotorg.plot_hit_trajectory(smear=self.smear)
            fit_func = hitplotorg.fit_hits_pol2(gr_smear)
            # Append new values. 
            d0_ls.append(hitplotorg.d0)
            daOvera_ls.append(hitplotorg.daOvera)
            dpTOverpT_ls.append(hitplotorg.dpTOverpT)
            
        self.d0_ls = d0_ls
        self.daOvera_ls = daOvera_ls
        self.dpTOverpT_ls = dpTOverpT_ls
    
    def check_kinem_name(self, kinem):
        """Make sure kinem is either 'dpTOverpT' or 'daOvera'."""
        assert kinem in ["dpTOverpT", "daOvera"]
        
    def make_label(self, kinem):
        """Return the LaTeX string of kinem."""
        self.check_kinem_name(kinem)
        if kinem == "dpTOverpT":
            return r"-#Deltap_{T}/p_{T}"
        elif kinem == "daOvera":
            return r"#Deltaa/a"
        
    def get_kinem_vals(self, kinem):
        """Return the values corresponding to kinem."""
        self.check_kinem_name(kinem)
        if kinem == "dpTOverpT":
            return np.array(self.dpTOverpT_ls) * -1.0
        elif kinem == "daOvera":
            return np.array(self.daOvera_ls)
        
    def graph_dXOverX(self, kinem="dpTOverpT", color=4, draw=True, do_fit=True):
        """Returns a plot of dX/X vs. d0."""
        err_msg = "You have to toss some toys first! Do: `biasplotter.toss_toys(int)`"
        assert self.n_toys is not None, err_msg
        # print "Making {} plot after tossing {} toys...".format(kinem, self.n_toys)
        print(f"Making {kinem} plot after tossing {self.n_toys} toys...")
        
        y_label = self.make_label(kinem)
        x_vals = np.array(self.d0_ls)
        y_vals = self.get_kinem_vals(kinem)
        graph_title = r"Bias of %s vs. d_{0} #bf{(p_{T}= %s GeV, n_{toys}= %d)}" % (y_label, self.pT, self.n_toys)
        
        # Then convert them to array.arrays afterward.
        x_arr = array('f', x_vals)
        y_arr = array('f', y_vals)
        
        # Make the graph.
        gr = r.TGraph(self.n_toys, x_arr, y_arr)
        gr.SetMarkerStyle(20)
        gr.SetMarkerSize(0.2)
#         gr.Draw("AP")
        gr.SetLineColor(color)
        gr.SetLineWidth(1)
        gr.SetMarkerColor(color)
        gr.SetTitle(graph_title)
        gr.GetXaxis().SetTitle(r"d_{0} [cm] (p0 from fit)")
        gr.GetYaxis().SetTitle(y_label)
        gr.GetXaxis().SetLimits(-0.05, 0.05)
        gr.GetYaxis().SetRangeUser(-0.2, 0.2)
#         gr.GetYaxis().SetTitleOffset(1.5)
        # Must do 2 Draw() calls to set TGraphs properly...
        if (draw):
            gr.Draw("AP")
        
        if (do_fit):
            # FIXME: Instead of self.fit_line = ...,
            # maybe try something like gr_dict[kinem] = ...
            # since this gets called twice for the same gr.
            self.fit_line = self.fit_graph_with_line(gr)
        return gr

    def fit_graph_with_line(self, graph):
        """Return a linear fit function and after fitting it to graph."""
        x_min = min(self.d0_ls) 
        x_max = max(self.d0_ls)
        # x = r.RooRealVar("x", "d_{0} (impact parameter)", x_min, x_max, "cm")  # (name, title, min, max, units)
        # interc = r.RooRealVar("p0", "intercept", -999, 999)
        # slope = r.RooRealVar("p1", "slope", -9999, 9999)
        # linefit = r.RooFormula("linefit", "@0+@1*@2", r.RooArgList(interc, slope, x))
        fit_func = r.TF1('f1', '[0]+[1]*x', x_min, x_max)
        fit_func.SetLineColor(2)
        fit_func.SetLineWidth(1)
        fit_func.SetLineStyle(2)
        # Fit it onto a histogram `h1`:
        graph.Fit(fit_func,'S')
        # The option 'S' saves the fit results into a pointer.
        r.gStyle.SetOptFit(111)
        fit_func.Draw("same")
        return fit_func
    
    def make_and_fill_TH2F(self, kinem="dpTOverpT", x_bins=100, y_bins=100, x_lim=None, y_lim=None, draw=True):
        """Make a 2D hist of dX/X vs. d0.
        
        x_lim : 2-elem list, [x_min, x_max]
            Used to determine range of x-axis.
        """
        y_label = self.make_label(kinem)
        title = r"Bias of %s vs. d_{0} #bf{(p_{T}= %s GeV, n_{toys}= %d)}" % (y_label, self.pT, self.n_toys)
        
        x_vals = np.array(self.d0_ls)
        y_vals = self.get_kinem_vals(kinem)
        
        def decide_limits(user_limits, vals):
            """If User specifies limits, then use those. Otherwise take limits from vals."""
            if user_limits is None:
                min_ = vals.min()
                max_ = vals.max()
            else: 
                min_ = user_limits[0]
                max_ = user_limits[1]
            return (min_, max_)
        x_min, x_max = decide_limits(x_lim, x_vals)
        y_min, y_max = decide_limits(y_lim, y_vals)
        
        # Make TH2F.
        h2 = r.TH2F("h2_%s_%s" % (kinem, self.pT), title, 
              x_bins, x_min, x_max,  # num_bins_x, x_min, x_max
              y_bins, y_min, y_max   # num_bins_y, y_min, y_max
             )
        
        # Fill TH2F.
        assert len(x_vals) == len(y_vals)
        for x,y in zip(x_vals, y_vals):
            h2.Fill(x, y, 1)
        
        # Decorate TH2F.
        h2.GetXaxis().SetTitle(r"d_{0} [cm] (p0 from fit)")
        h2.GetYaxis().SetTitle(y_label)
        
        r.gStyle.SetOptStat(0)   # Don't show the stats box.
        r.gStyle.SetPalette(55)  # Change the color map.
        # More color maps: https://root.cern.ch/doc/master/classTColor.html#C06
        h2.SetContour(200)
        #  h2.GetZaxis().SetRangeUser(0, 15)  # Restrict color bar to have limits: [0, 15].
        
        if (draw):
            h2.Draw("colz")
        return h2