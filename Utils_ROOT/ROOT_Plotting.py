import ROOT as r

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

def make_TPave(w, h, topright_corner_pos=(0.33, 0.62)):#, when="Before"):
    """Return a TPaveText that shows either before or after corr stats.
    
    Parameters
    ----------
    topright_corner_pos : 2-tup
        (x_max, y_max) coordinate of pave.
    when : str
        Should be either 'Before' or 'After'
    """
    # color = r.kBlue if when.lower() in "before" else r.kRed
    x_right, y_top = topright_corner_pos[0], topright_corner_pos[1]
    assert x_right >= w
    assert y_top >= h
    x_left = x_right - w
    y_bot = y_top - h

    pave = r.TPaveText(x_left, y_bot, x_right, y_top, "NDC")
    # pave = r.TPaveText(0.13, text_y_min-0.11, 0.33, text_y_min, "NDC")
    pave.SetFillColor(0)
    pave.SetBorderSize(1) # Use 0 for no border.
    pave.SetTextAlign(12) # 22 is centered vert and horiz.
    pave.SetTextSize(0.016)
    # pave.SetTextColor(1)
    pave.SetFillStyle(1001)  # Solid fill.
    # pave.AddText(f"Used %i Gaus fit iterations per fit" % iters)  # Accommodates LaTeX!
    # pave.SetTextColor(color)
    # pave.AddText(r"%s p_{T} corrections:" % when)
    # pave.AddText(f"  #mu = {mean.getVal():.4g} #pm {mean.getError():.4g}")
    # pave.AddText(f"  #sigma = {sigma.getVal():.4g} #pm {sigma.getError():.4g}")
    # pave.AddText(r"  #chi^{2}/ndf = %.4g" % xframe.chiSquare())
    return pave