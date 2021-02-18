import ROOT as r
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid

def make_ratio_pads():
    """Create an upper pad (for hists), and a lower pad (for a ratio plot)."""
    ptop = r.TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0)
    pbot = r.TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25)
    ptop.SetBottomMargin(0)
    pbot.SetTopMargin(0)
    pbot.SetBottomMargin(0.25)
    # ptop.Draw()
    # pbot.Draw()
    # ptop.cd()
    # ptop.SetTicks(1,1)
    # pbot.cd()
    # pbot.SetTicks(1,1)
    return (ptop, pbot)

def Root_Hist_GetLastBinRightEdge(hist):
    """Returns the right-most bin edge of a ROOT.TH1 object."""
    n_bins = hist.GetNbinsX()
    edge = hist.GetBinLowEdge(n_bins) + hist.GetBinWidth(n_bins)
    return edge

def make_new_xframe(data_x_min, data_x_max, x_lim=None, x_label="", units="", n_bins=100, title=""):
    """Make and return 4-tuple of new RooFit frame info.
    
    Parameters
    ----------
    data_x_min : float
        The min val of the 1-dim dataset.
    data_x_max : float
        The min val of the 1-dim dataset.
    x_lim : 2-elem list
        The range of the plot's x-axis.
        [x_axis_min, x_axis_max]
    x_label : str
        Description of the x variable.
    units : str
        Units of the x variable.
    n_bins : int
        Number of bins for x variable.
    title : str
        Title to be put onto plot once x variable is plotted.
    """
    # A frame was not provided, so make and decorate one.
    x_min = data_x_min if x_lim is None else x_lim[0]
    x_max = data_x_max if x_lim is None else x_lim[1]
    # Changing RooRealVar name: # x_roofit = r.RooRealVar("x_roofit", x_label, x_min, x_max, units)  # (name, title, min, max, units)
    x_roofit = r.RooRealVar("x_var", x_label, x_min, x_max, units)  # (name, title, min, max, units)
    # x_roofit.setRange("test_range_again", x_min, x_max)
    x_roofit.setBins(n_bins)
    # xframe = x_roofit.frame(r.RooFit.Title(title)) #r.RooFit.Range(x_min, x_max)
    xframe = x_roofit.frame(r.RooFit.Title(title), r.RooFit.Range(x_min, x_max))
    return (x_min, x_max, x_roofit, xframe)

def set_plot_styles(gridOn=True):
    """Make your plots consistently pretty."""
    # r.gROOT.SetBatch(r.kTRUE)
    # Plotting info.
    tdrStyle = setTDRStyle()
    tdrGrid(tdrStyle, gridOn=gridOn)