"""A module to manage Printers in ROOT (like TCanvas)."""
import ROOT as r
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid, fixOverlay

class CanvasPrinter:
    """A class to manage the printing (drawing) of plots onto a canvas."""

    def __init__(self, show_plots=False, show_statsbox=True, canv=None):
        r.gROOT.SetBatch(not show_plots)  # SetBatch(True) means don't show plots.
        self.style = self.make_plots_pretty(show_statsbox=show_statsbox)
        # self.canv = r.TCanvas("c1", "c1", 550, 600) if canv is None else canv
        self.canv = r.TCanvas() if canv is None else canv

    def make_plots_pretty(self, show_statsbox=True,
                                pad_top_margin=0.10,
                                pad_bottom_margin=0.10,
                                pad_left_margin=0.15,
                                pad_right_margin=0.05):
        """Activate a TStyle which sets a consistent style for all plots."""
        tdrStyle = setTDRStyle(show_statsbox=show_statsbox,
                pad_top_margin=pad_top_margin,
                pad_bottom_margin=pad_bottom_margin,
                pad_left_margin=pad_left_margin,
                pad_right_margin=pad_right_margin)
        tdrGrid(tdrStyle, gridOn=False)
        return tdrStyle

    def make_pdf_of_plots(self, plot_ls, outpdf_path):#, use_TDR_style=True):
        """For each plot in plot_ls, draw the plots to outpdf_path."""
        print(f"Drawing plots to:\n{outpdf_path}")
        self.canv.Print(outpdf_path + "[")
        self.draw_plots(plot_ls, outpdf_path=outpdf_path)
        self.canv.Print(outpdf_path + "]")

    def draw_plots(self, plot_ls, outpdf_path, show_statsbox=True, add_plots_to_pdf=True, canv=None):
        """For each plot in plot_ls, draw the plots.
        
        If add_plots_to_pdf, then draw each plot on its own page in the
        growing pdf. The pdf will be saved at outpdf_path.
        """
        for plot in plot_ls:
            if isinstance(plot, r.TH1):
                plot.Draw("hist")
            if isinstance(plot, r.RooPlot):
                plot.Draw()
            if isinstance(plot, r.TH2):
                # Don't show under-overflow bin stats matrix.
                r.gStyle.SetOptStat("iRMe")
                plot.Draw("colz")
                # self.canv.Update()
            if isinstance(plot, r.TGraph):
                plot.Draw("AP")
            if add_plots_to_pdf:
                assert len(outpdf_path) > 0, "You must specify an outpdf_path."
                self.canv.Print(outpdf_path) if canv is None else canv.Print(outpdf_path)
            r.gPad.Update()
            # Ensure original style in case of any changes e.g. TH2 above.
            self.make_plots_pretty(show_statsbox=show_statsbox)