import numpy as np
import ROOT
from array import array

def make_TH1F(internal_name, title=None, xlabel=None, n_bins=100, x_min=0, x_max=10, units=None):
    """A function to quickly make TH1F. Return the TH1F.
    
    NOTE: 
    - Don't worry about parentheses around units.
      It will be taken care of.
    - Accommodates ROOT LaTeX (using '#' of course).
      Use raw strings: r'your string'.

    Parameters
    ----------
    internal_name : str
        asdfasdfasdf
    title : str
        The title of the histogram.
    """
    h = ROOT.TH1F(internal_name, "", n_bins, x_min, x_max)
    bin_w = (x_max - x_min) / float(n_bins)
    ylabel = r"Events / (%s)" % bin_w
    if units is not None:
        xlabel += r" (%s)" % units
        ylabel = ylabel.rstrip(")") + r" %s)" % units
    h.SetTitle(xlabel) if title is None else h.SetTitle(title)
    h.SetXTitle(xlabel)
    h.SetYTitle(ylabel)
    return h

def make_TGraphErrors(x, y, x_err=None, y_err=None,
                      use_binwidth_xerrs=False, title="",
                      x_label="", x_min=0, x_max=10, x_units=None,
                      y_label="", y_min=0, y_max=10, y_units=None,
                      line_color=1, line_width=2,
                      marker_color=1, marker_style=21, marker_size=0.5):
    """A function to quickly make a TGraphErrors. Return the TGraphErrors.
    
    NOTE: Don't pass in x_err yet. Not ready.
    """
    assert x_err is None
    # Make sure the lengths are all the same
    n_pts = len(x)
    assert all([len(arr) == n_pts for arr in [y, x_err, y_err] if arr is not None])
    # Convert to arrays since that's what ROOT likes.
    yerr_tmp = np.zeros_like(x) if y_err is None else y_err
    xerr_tmp = np.zeros_like(x) if x_err is None else x_err
    x_arr = array('f', x)
    y_arr = array('f', y)
    xerr_arr = array('f', xerr_tmp)
    yerr_arr = array('f', yerr_tmp)
    # FIXME: Implement this later:
    # if (x_err is None) and not use_binwidth_xerrs:
    #     x_err = 
    # if use_binwidth_xerrs:
    #     low_err_arr, high_err_arr = calc_x_err_bins_from_bin_edges(binedge_ls)
    # ROOT.TGraphAsymmErrors()
    gr = ROOT.TGraphErrors(n_pts, x_arr, y_arr, xerr_arr, yerr_arr)
    gr.SetLineColor(line_color)
    gr.SetLineWidth(line_width)
    gr.SetMarkerColor(marker_color)
    gr.SetMarkerStyle(marker_style)
    gr.SetMarkerSize(marker_size)
    gr.SetTitle(title)
    def add_units(label, units):
        return label + r" (%s)" % units if len(units) > 0 else label
    x_label = add_units(x_label, x_units)
    y_label = add_units(y_label, y_units)
    gr.GetXaxis().SetTitle(x_label)
    gr.GetYaxis().SetTitle(y_label)
    gr.GetXaxis().SetLimits(x_min, x_max)  # SetRangeUser() doesn't work for x-axis! SetLimits() instead.
    gr.GetYaxis().SetRangeUser(y_min, y_max)
    return gr

def make_TLegend(x_dim=(0,1), y_dim=(0,1), screenshot_dim=None, buffer_dim=None):
    """Return a TLegend, given its x and y coord.
        
    Parameters
    ----------
    x_dim : 2-tup
        (x_min, x_max) of TLegend position as a fraction of the TCanvas.
        So x_min = 0.2 means that the TLegend will begin at 20% of the
        canvas width.
    y_dim : 2-tup
        (y_min, y_max) of TLegend position as a fraction of the TCanvas.
    screenshot_dim : 2-tup
        (width, height) of full screenshot window, in pixels.
        Put the window around the entire canvas.
    buffer_dim : 4-tup
        (left_offset, right_offset, bottom_offset, top_offset),
        measured in pixels.
        Use another screenshot window to measure the offset from all four
        edges of the canvas to where you want your legend to be.
    """
    if screenshot_dim is not None:
        assert buffer_dim is not None
        canv_width = screenshot_dim[0]
        canv_height = screenshot_dim[1]
        left_offset = buffer_dim[0]
        right_offset = buffer_dim[1]
        bot_offset = buffer_dim[2]
        top_offset = buffer_dim[3]
        x_min, x_max, y_min, y_max = get_normcoord_from_screenshot(canv_width, canv_height,
                                  left_offset, right_offset, bot_offset, top_offset)
    else:
        x_min, x_max = x_dim[0], x_dim[1]
        y_min, y_max = y_dim[0], y_dim[1]
    return ROOT.TLegend(x_min, y_min, x_max, y_max)

def make_TMultiGraph_and_Legend(gr_ls=[], leg_txt_ls=[], y_min=None, y_max=None):
    """Return a (TMultiGraph, TLegend) with all the TGraphs from gr_ls added."""
    mg = ROOT.TMultiGraph()
    for gr in gr_ls:
        mg.Add(gr, "p")
    x_title = gr_ls[0].GetXaxis().GetTitle()
    y_title = gr_ls[0].GetYaxis().GetTitle()
    title = gr_ls[0].GetTitle()
    all_titles = "%s;%s;%s" % (title, x_title, y_title)
    mg.SetTitle(all_titles)

    leg = make_TLegend(x_dim=(0,1), y_dim=(0,1),
                    screenshot_dim=(878,872), buffer_dim=(176,438,610,131))
    # leg_text = r"%.2f < #deltam_{4#mu}/m_{4#mu} < %.2f%%" % (relmin, relmax)
    leg.SetTextSize(0.02)
    for gr, txt in zip(gr_ls, leg_txt_ls):
        leg.AddEntry(gr, txt, "lpfe")
    # leg.SetLineWidth(3)
    leg.SetBorderSize(1)
    if (y_min is not None) and (y_max is not None):
        mg.SetMinimum(y_min) # Change y-axis limits.
        mg.SetMaximum(y_max)
    return (mg, leg)

def get_normcoord_from_screenshot(canv_width, canv_height,
                                  left_offset, right_offset, bot_offset, top_offset):
    """Return a 4-tup of normalized coordinates, given screenshot dimensions.

    NOTE: 
        MacOS has a nice screenshot feature (Cmd + Shift + 5).
        It shows the dimensions of the screenshot window.
        You can use this window as a ruler of sorts, in units of pixels.
        In this function, use the window to measure the dimensions of your
        canvas and then to measure the offsets from the canvas edges.
        This will get you the canvas-normalized dimensions of your object.

    Parameters
    ----------
    canv_height : float
    canv_width : float
    left_offset : float
    right_offset : float
    bot_offset : float
    top_offset : float
    """
    # A fraction of the total width.
    x_min = left_offset / float(canv_width)
    x_max = 1 - (right_offset / float(canv_width))
    y_min = bot_offset / float(canv_height)
    y_max = 1 - (top_offset / float(canv_height))
    return (x_min, x_max, y_min, y_max)

# def add_branch_pTcorr_fromadhocmethod():
#     """Copy  a '.root' file 