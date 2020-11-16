"""A module to manage Printers in ROOT (like TCanvas)."""
import ROOT as r
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid, fixOverlay

class CanvasPrinter:
    """A class to manage the printing (drawing) of plots onto a canvas."""

    def __init__(self, show_plots=False, gridOn=True, canv=None):
        r.gROOT.SetBatch(not show_plots)  # SetBatch(True) means don't show plots.
        self.make_plots_pretty(gridOn)
        # self.canv = r.TCanvas("c1", "c1", 550, 600) if canv is None else canv
        self.canv = r.TCanvas() if canv is None else canv

    def make_plots_pretty(self, gridOn=True):
        """Activate a TStyle which sets a consistent style for all plots."""
        tdrStyle = setTDRStyle()
        tdrGrid(tdrStyle, gridOn=gridOn)
        # return tdrStyle

    def make_pdf_of_plots(self, plot_ls, outpdf_path):#, use_TDR_style=True):
        """For each plot in plot_ls, draw the plots to outpdf_path."""
        print(f"Drawing plots to:\n{outpdf_path}")
        self.canv.Print(outpdf_path + "[")
        for plot in plot_ls:
            if isinstance(plot, r.RooPlot):
                plot.Draw()
            if isinstance(plot, r.TH1):
                plot.Draw("hist")
            if isinstance(plot, r.TH2):
                # Don't show under-overflow bin stats matrix.
                r.gStyle.SetOptStat("iRMe")
                plot.Draw("colz")
                # self.canv.Update()
            if isinstance(plot, r.TGraph):
                plot.Draw("AP")
            self.canv.Print(outpdf_path)
            r.gPad.Update()
            # Ensure original style in case of any changes e.g. TH2 above.
            self.make_plots_pretty()
        self.canv.Print(outpdf_path + "]")