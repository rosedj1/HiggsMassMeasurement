""" Theoretical Muon Trajectory Analyzer
PURPOSE: This script is used to make plots of simulated muon track hits
    in the Pixel and Strip Tracker layers. 
SYNTAX: python this_script.py
NOTES: This code produces 4 kinds of plots: 
    (1) A one-page PDF which shows the trajectories of muons with 
    varying pTs.
    (2) An X-page PDF is made. On each page, a plot is made of the 
    true muon track and of a y-val-smeared muon track (reco track). 
    (3) An X-page PDF is made. On each page, a delta_a/a vs. d0
    scatterplot is made on the left side and a delta_pT/pT vs. d0 
    on the right side for each pT muon.
    (4) Instead of the scatterplot like in (3), instead a TH2F is made
    to detect an anomalous 'clumping" of points.
AUTHOR: Jake Rosenzweig
CREATED: 2020-08-20
UPDATED: 2020-08-27
TO DO: 
    [ ] Figure out how to suppress fit stats output.
    [ ] Perform a pol2 fit in make_pdf_smear_and_nonsmear_traj().
"""
import os
import ROOT as r
# Package imports.
from Toy_Model.utils.classes import HitPlotOrg, BiasPlotter
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid

#--- User Parameters ---#
pdf_outpath = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots"
# Processes to run. Structure: (bool_to_run_this, output_path)
do_muon_traj_cumul_plot = (0, os.path.join(pdf_outpath, "muon_trajectories.pdf"))
do_smear_and_nonsmear_plots = (0, os.path.join(pdf_outpath, "smear_and_nonsmear_traj_eachpT.pdf"))
do_daOvera_dpTOverpT_plots = (1, os.path.join(pdf_outpath, "deltaa_and_deltapT_vs_d0_scatterplots_test01.pdf"))
do_th2f_scatterplots = (0, os.path.join(pdf_outpath, "h2d_withbestfitlines.pdf"))

n_toys = 100
pT_ls = [5, 10, 20, 30, 40, 50, 75, 100]
new_pixel_pos = [2.9, 6.8, 10.5, 16.0, 25.5, 33.9, 41.9, 49.8, 60.8, 69.2, 78.0, 86.8, 96.5, 108.0]

#------------------------#
#--- Script functions ---#
#------------------------#
def plot_muon_hits_all_pTs(pT_ls, position_ls, smear=False):
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
    c1 = r.TCanvas()
    leg = r.TLegend(0.2, 0.7, 0.4, 0.9)
    mg = r.TMultiGraph()

    for count,pT in enumerate(pT_ls, 1):
        hitplotorg = HitPlotOrg(pT=pT, position_ls=position_ls)
        gr, leg = hitplotorg.plot_hit_trajectory(leg=leg, color=count, smear=smear, sigma=0.01)
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

def plot_smear_and_nonsmear_yvals(pT, position_ls, canv=None):
    """
    Show one plot with all muon hits for each pT and Tracker position specified.
    
    Parameters
    ----------
    pT : list
        The muon pT whose track will be plotted.
    position_ls : list
        The hits along the Tracker as a radial distance from the interaction point.
        
    Returns
    -------
    mg : TMultiGraph
        Contains all the TGraphs, each of which is a muon trajectory.
    leg : TLegend
        The legend.
    canv : TCanvas
        The drawn-on canvas.
    """
    if canv is None:
        canv = r.TCanvas()
    leg = r.TLegend(0.2, 0.8, 0.55, 0.9)
    mg = r.TMultiGraph()

    hitplotorg1 = HitPlotOrg(pT=pT, position_ls=position_ls)
    hitplotorg2 = HitPlotOrg(pT=pT, position_ls=position_ls)
    gr, leg = hitplotorg1.plot_hit_trajectory(leg=leg, color=2, smear=False, sigma=0.01)  
    gr_smear, leg_smear = hitplotorg2.plot_hit_trajectory(leg=leg, color=4, smear=True, sigma=0.01)
    gr.SetMarkerStyle(2)
    gr_smear.SetMarkerStyle(5)
    mg.Add(gr, "CP")
    mg.Add(gr_smear, "CP")
    x_label = r"Transverse Pixel/Strip Positions [cm]"
    y_label = r"Distance from x-axis [cm]"
    title = r"p_{T}= %s GeV: Before and after y-val smearing (#sigma_{y}=100 #bf{#mu}m)" % pT
    mg_title = r"%s;%s;%s" % (title, x_label, y_label)
    mg.SetTitle(mg_title)
    mg.Draw("A")
    mg.GetXaxis().SetRangeUser(0, 120)
    mg.GetYaxis().SetRangeUser(-2, 14)
    leg.SetTextSize(0.02)
    leg.Draw("same")
#     canv.Draw()
    return mg, leg

def set_properties():
    """Adjust plot properties by calling on setTDRStyle()."""
    r.gROOT.SetBatch(True)
    tdrStyle = setTDRStyle()
    tdrGrid(tdrStyle, gridOn=True)

def make_pdf_smear_and_nonsmear_traj(outfile_path, pT_ls, position_ls):
    """
    Make multi-page PDF. Each page shows smear vs. non-smear trajectory.
    
    TO DO: [ ] Perform a pol2 fit.
    """
    mg_ls = []
    leg_ls = []

    c1 = r.TCanvas()
    c1.Print(outfile_path + "[")
    for pT in pT_ls:
        mg, leg = plot_smear_and_nonsmear_yvals(pT, position_ls, canv=c1)
        mg_ls.append(mg)
        leg_ls.append(leg)
        c1.Draw()
        c1.Print(outfile_path)

    c1.Print(outfile_path + "]")

def make_pdf_daOvera_dpTOverpT_plots(outfile_path, pT_ls, position_ls, n_toys):
    """Make multi-page PDF. Each page has a delta_a/a plot on left and delta_pT/pT on right."""
    c1 = r.TCanvas()
    biasplotter_ls = []
    gr_daOvera_ls = []
    gr_dpTOverpT_ls = []
    c1.Print(outfile_path + "[")
    c1.Divide(2,1)
    for pT in pT_ls:
        
        biasplotter_ls.append( BiasPlotter(pT=pT, position_ls=position_ls) )
        biasplotter_ls[-1].toss_toys(n_toys=n_toys)

        c1.cd(1)
        gr_daOvera_ls.append( biasplotter_ls[-1].graph_dXOverX(kinem="daOvera", color=4, draw=True, do_fit=True) )
        c1.cd(2)
        gr_dpTOverpT_ls.append( biasplotter_ls[-1].graph_dXOverX(kinem="dpTOverpT", color=4, draw=True, do_fit=True) )
        c1.Draw()
        c1.Print(outfile_path)

    c1.Print(outfile_path + "]")

def make_TH2F_plots(outfile_path, pT_ls, position_ls, n_toys, x_lim=[-0.05, 0.05], y_lim=[-0.2, 0.2]):
    """Make multi-page PDF. Each page shows TH2F of da/a on left and TH2F of dpT/pT on right."""
    c1 = r.TCanvas()
    biasplotter_ls = []
    gr_daOvera_ls = []
    gr_dpTOverpT_ls = []
    h2_ls = []
    c1.Print(outfile_path + "[")
    c1.Divide(2,1)

    for pT in pT_ls:
        biasplotter_ls.append( BiasPlotter(pT=pT, position_ls=position_ls) )
        biasplotter_ls[-1].toss_toys(n_toys=n_toys)

        c1.cd(1)
        gr_daOvera_ls.append( biasplotter_ls[-1].graph_dXOverX(kinem="daOvera", color=4, draw=False, do_fit=True) )
        h2_ls.append( biasplotter_ls[-1].make_and_fill_TH2F(kinem="daOvera", x_bins=100, y_bins=100, x_lim=x_lim, y_lim=y_lim, draw=True) )
        biasplotter_ls[-1].fit_line.Draw("same")
        
        c1.cd(2)
        gr_daOvera_ls.append( biasplotter_ls[-1].graph_dXOverX(kinem="dpTOverpT", color=4, draw=False, do_fit=True) )
        h2_ls.append( biasplotter_ls[-1].make_and_fill_TH2F(kinem="dpTOverpT", x_bins=100, y_bins=100, x_lim=x_lim, y_lim=y_lim, draw=True) )
        biasplotter_ls[-1].fit_line.Draw("same")
        
        c1.Draw()
        c1.Print(outfile_path)

        c1.Print(outfile_path + "]")

if __name__ == '__main__':
    set_properties()

    if (do_muon_traj_cumul_plot[0]):
        mg, leg, c1 = plot_muon_hits_all_pTs(pT_ls, new_pixel_pos, smear=False)
        c1.SaveAs(do_muon_traj_cumul_plot[1])

    if (do_smear_and_nonsmear_plots[0]):
        # NOTE: Finish this function by fitting with pol2.
        outfile = do_smear_and_nonsmear_plots[1]
        make_pdf_smear_and_nonsmear_traj(outfile, pT_ls, new_pixel_pos)

    if (do_daOvera_dpTOverpT_plots[0]):
        outfile = do_daOvera_dpTOverpT_plots[1].replace(".pdf", f"_{n_toys}toys.pdf")
        make_pdf_daOvera_dpTOverpT_plots(outfile, pT_ls, new_pixel_pos, n_toys)

    if (do_th2f_scatterplots[0]):
        outfile = do_th2f_scatterplots[1].replace(".pdf", f"_{n_toys}toys.pdf")
        make_TH2F_plots(outfile, pT_ls, new_pixel_pos, n_toys, x_lim=[-0.05, 0.05], y_lim=[-0.2, 0.2])