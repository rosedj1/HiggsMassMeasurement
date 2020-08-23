"""
PURPOSE: This script is used to make plots of simulated muon track hits
    in the Pixel and Strip Tracker layers. 
SYNTAX: python this_script.py
NOTES:
AUTHOR: Jake Rosenzweig
CREATED: 2020-08-20
UPDATED: 2020-08-21
"""
from ROOT import TGraph, TCanvas, TLegend, TF1, TMultiGraph, gROOT
import numpy as np
from array import array

from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid

new_pixel_pos = [2.9, 6.8, 10.5, 16.0, 25.5, 33.9, 41.9, 49.8, 60.8, 69.2, 78.0, 86.8, 96.5, 108]
pT_ls = [5, 10, 20, 40, 50, 75, 100]

class HitPlotOrg:
    def __init__(self, pT, position_ls):
        """Plot and fit the muon track."""
        self.pT_true = pT
        self.position_ls = position_ls 
        
        self.a_true = self.convert_pT_to_a(pT)
        self.a_fit = None
        self.pT_calc = None
        
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
        
    def find_parabolic_y_vals(self, a, x_arr):
        """Return the y-values, assuming a parabolic trajectory."""
        return a * x_arr**2
    
    def draw_hits(self, leg=None, color=1):
        """
        Draw the muon track as it hits the Pixel and Strip layers.
        Return a TGraph object.
        """
        n_pts = len(self.position_ls)
        
        # Use y = a*x^2 to model the trajectory of the muon.
        # For high pT this is reasonable. 
        # Radius of curvature goes like: R = 1 / (2*a)
        x_vals = np.array(self.position_ls)
        y_vals = self.find_parabolic_y_vals(self.a_true, x_vals)
        
        # Then convert them to array.arrays afterward.
        x_arr = array('f', x_vals)
        y_arr = array('f', y_vals)

        # Make the graph.
        gr = TGraph(n_pts, x_arr, y_arr)
        gr.SetLineColor(color)
        gr.SetLineWidth(1)
        gr.SetMarkerColor(color)
        gr.SetMarkerStyle( 21 )
        gr.SetTitle('Hits along muon path')
        gr.GetXaxis().SetTitle('Transverse Pixel/Strip Positions [cm]')
        gr.GetYaxis().SetTitle('Distance from x-axis [cm]')
        gr.GetYaxis().SetTitleOffset(1.5)
        
        # Add to the legend.
        if leg is None:
            leg = TLegend()
        text = r"p_{T} = %.0f GeV" % self.pT_true
        leg.AddEntry(gr, text, "lp")
#         leg.SetTextColor(color)
        return gr, leg

#--- Script functions ---#
def plot_muon_hits_all_pTs(pT_ls, position_ls):
    """
    Show one plot with all muon hits for each pT and Tracker position specified.
    
    Parameters
    ----------
    pT_ls : list
        The muon pTs whose tracks will be plotted.
    position_ls : list
        The hits along the Tracker as a radial distance from the interaction point.
        
    Returns
    -------
    mg : TMultiGraph
        Contains all the TGraphs, each of which is a muon trajectory.
    leg : TLegend
        The legend.
    c1 : TCanvas
        The drawn-on canvas.
    """
    c1 = TCanvas()
    leg = TLegend(0.2, 0.7, 0.4, 0.9)
    mg = TMultiGraph()

    for count,pT in enumerate(pT_ls, 1):
        hitplotorg = HitPlotOrg(pT=pT, position_ls=position_ls)
        gr, leg = hitplotorg.draw_hits(leg=leg, color=count)
        mg.Add(gr, "CP")
    mg_title = "Hits along muon path; Transverse Pixel/Strip Positions [cm];Distance from x-axis [cm]"
    mg.SetTitle(mg_title)
    mg.Draw("A")
    mg.GetXaxis().SetRangeUser(0, 120)
    # mg.GetYaxis().SetRangeUser(0, 120)
    leg.SetTextSize(0.02)
    leg.Draw("same")
    c1.Draw()
    return mg, leg, c1

def set_properties():
    gROOT.SetBatch(True)
    tdrStyle = setTDRStyle()
    tdrGrid(tdrStyle, gridOn=True)

if __name__ == '__main__':
    set_properties()
    mg, leg, c1 = plot_muon_hits_all_pTs(pT_ls, new_pixel_pos)
    c1.SaveAs("/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/toy_model/muon_hits_thru_Tracker.pdf")
    
    # fit_func = TF1('f1', '[0]+[1]*x+[2]*x^2', 0, new_pixel_pos[-1])