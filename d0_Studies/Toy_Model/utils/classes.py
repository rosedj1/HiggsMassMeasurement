import numpy as np
import ROOT as r
from array import array

class HitPlotOrg:
    """
    This class represents a muon trajectory.
    """

    def __init__(self, pT, position_ls):
        """Plot and fit the muon track."""
        self.pT_true = pT
        self.position_ls = position_ls 
        self.n_hits = len(position_ls)
        
        self.a_true = self.convert_pT_to_a(pT)
        self.a_fit = None
        self.a_fit_err = None
        self.d0 = None
        self.pT_fit = None
        self.y_vals = None
        self.y_vals_smear = None
        
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
        return 0.3 * abs(q) * B / float(2.0 * pT * 100)
    
    def convert_a_to_pT(self, a, q=-1, B=3.8):
        """
        Converts `a` in pT, doing the reverse of convert_pT_to_a().
        Assumes a is in cm^(-1).
        Returns pT in GeV/c.
        """
        return 0.3 * abs(q) * B / float(2.0 * a * 100)
    
    def calc_daOvera(self):
        """Return (a_fit - a_true) / a_true."""
        return float(self.a_fit - self.a_true) / self.a_true
    
    def calc_dpTOverpT(self):
        """Return (pT_fit - pT_true) / pT_true."""
        return float(self.pT_fit - self.pT_true) / self.pT_true
        
    def find_parabolic_y_vals(self, a, x_arr):
        """Return the y-values, assuming a parabolic trajectory."""
        return a * x_arr**2
    
    def plot_hit_trajectory(self, leg=None, color=1, smear=False, sigma=0.01):
        """
        Draw the muon track as it hits the Pixel and Strip layers.
        Return (TGraph, TLegend) tuple.
        """
        # Add to the legend.
        if leg is None:
            leg = r.TLegend()
        leg_text = r"p_{T} = %.0f GeV" % self.pT_true
            
        n_pts = len(self.position_ls)
        
        # Use y = a*x^2 to model the trajectory of the muon.
        # For high pT this is reasonable. 
        # Radius of curvature goes like: R = 1 / (2*a)
        x_vals = np.array(self.position_ls)
        y_vals = self.find_parabolic_y_vals(self.a_true, x_vals)
        self.y_vals = y_vals.copy()  # self.y_vals will get smear later. 
        
        # Then convert them to array.arrays afterward.
        x_arr = array('f', x_vals)
        
        graph_title = "Hits along muon path"
        if (smear):
            graph_title += " (y-val smeared)"
            leg_text += " (y-val smeared)"
            self.y_vals_smear = self.smear_hits(sigma=sigma)
            y_arr = array('f', self.y_vals_smear)
        else:
            y_arr = array('f', y_vals)

        # Make the graph.
        gr = r.TGraph(n_pts, x_arr, y_arr)
        gr.SetLineColor(color)
        gr.SetLineWidth(1)
        gr.SetMarkerColor(color)
        gr.SetMarkerStyle( 21 )
        gr.SetTitle(graph_title)
        gr.GetXaxis().SetTitle('Transverse Pixel/Strip Positions [cm]')
        gr.GetYaxis().SetTitle('Distance from x-axis [cm]')
        gr.GetYaxis().SetTitleOffset(1.5)
        
        leg.AddEntry(gr, leg_text, "lp")
#         leg.SetTextColor(color)
        return gr, leg

    def smear_hits(self, sigma=0.01):
        """
        Return the y-value of each hit after adding a small Gaussian uncertainty.
        By default, sigma of Gaus is 100 um = 0.01 cm (since positions are measured in cm).
        """
        smears = np.random.normal(loc=0, scale=sigma, size=self.n_hits)
#         rng = np.random.default_rng(1)
#         smears = rng.normal(loc=0, scale=sigma, size=self.n_hits)
        return self.y_vals + smears

    def fit_hits_pol2(self, graph):
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
        x_min = min(self.position_ls) 
        x_max = max(self.position_ls)
        fit_func = r.TF1('f1', '[0]+[1]*x+[2]*x^2', x_min, x_max)
        fit_func.SetLineColor(1)
        fit_func.SetLineWidth(2)
        fit_func.SetLineStyle(2)
        # Fit it onto a histogram `h1`:
        graph.Fit(fit_func,'S')
        # The option 'S' saves the fit results into a pointer.
        r.gStyle.SetOptFit(111)
        
        self.set_fit_vals(fit_func)
        
        return fit_func
    
    def set_fit_vals(self, fit_func):
        """Check that the best-fit params are reasonable and set attributes."""
        self.a_fit = fit_func.GetParameter(2)
        self.a_fit_err = fit_func.GetParError(2)
        self.d0 = fit_func.GetParameter(0)
        assert self.a_fit is not None
        assert self.a_fit_err is not None
        
        self.pT_fit = self.convert_a_to_pT(self.a_fit)
        self.dpTOverpT = self.calc_dpTOverpT()
        self.daOvera = self.calc_daOvera()

class BiasPlotter:
    def __init__(self, pT, position_ls, smear=True):
        """Toss n_toys to test the theoretical bias of da/a or dpT/pT vs. d0."""
        self.pT = pT
        self.position_ls = position_ls
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
            hitplotorg = HitPlotOrg(pT=self.pT, position_ls=self.position_ls)
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
        print "Making {} plot after tossing {} toys...".format(kinem, self.n_toys)
        # print(f"Making {kinem} plot after tossing {self.n_toys} toys...")
        
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
            self.fit_line = self.fit_graph_with_line(gr)
        return gr

    def fit_graph_with_line(self, graph):
        """Return a linear fit function and after fitting it to graph."""
        x_min = min(self.d0_ls) 
        x_max = max(self.d0_ls)
        x = r.RooRealVar("x", "d_{0} (impact parameter)", x_min, x_max, "cm")  # (name, title, min, max, units)
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