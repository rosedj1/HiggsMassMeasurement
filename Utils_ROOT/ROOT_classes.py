import numpy as np
import ROOT
from array import array

#-----------#
#--- TH1 ---#
#-----------#
def make_TH1F(internal_name, title=None, n_bins=100,
              xlabel=None, x_min=0, x_max=10,
              ylabel=None, units=None):
    """A function to quickly make TH1F. Return the TH1F.
    
    NOTE: 
    - Don't worry about parentheses around units.
      It will be taken care of.
    - Accommodates ROOT LaTeX (use '#' instead of '\').
      Use raw strings: r'your string'.

    Parameters
    ----------
    internal_name : str
        ROOT's internal name for this histogram.
    title : str
        The title of the histogram.
    """
    h = ROOT.TH1F(internal_name, "", n_bins, x_min, x_max)
    bin_w = (x_max - x_min) / float(n_bins)
    ylabel = r"Events / (%s)" % bin_w if ylabel is None else ylabel
    xlabel = "" if xlabel is None else xlabel
    if units is not None:
        xlabel += r" (%s)" % units
        ylabel = ylabel.rstrip(")") + r" %s)" % units
    h.SetTitle(xlabel) if title is None else h.SetTitle(title)
    h.SetXTitle(xlabel)
    h.SetYTitle(ylabel)
    return h

#-----------#
#--- TH2 ---#
#-----------#
def make_TH2F(internal_name, title=None, 
              n_binsx=5, x_label="", x_units=None, x_min=0, x_max=10,
              n_binsy=5, y_label="", y_units=None, y_min=0, y_max=10,
              z_min=None, z_max=None, z_label_size=None,
              n_contour=100):
    """A function to quickly make TH2F. Return the TH2F.
    
    TODO: Update docstring.
    
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
    n_binsx : int or list
        Number of bins along the x-axis.
        If a list is provided, then the elements are the bin edges (which can
        be non-uniform!).
    x_min : float
        If n_binsx is type int, then the first bin will start at this x-val.
        If n_binsx is array-like, then this parameter is ignored.
    x_max : float
        If n_binsx is type int, then the last bin will end at this x-val.
        If n_binsx is array-like, then this parameter is ignored.
    n_binsy : int or list
        Number of bins along the y-axis.
        If a list is provided, then the elements are the bin edges (which can
        be non-uniform!).
    n_contour : int
        Number of contour lines in color label.
    """
    if isinstance(n_binsx, int) and isinstance(n_binsy, int):
        h2 = ROOT.TH2F(internal_name, "",
                       n_binsx, x_min, x_max,
                       n_binsy, y_min, y_max
        )
    elif isinstance(n_binsx, list) and isinstance(n_binsy, list):
        nx_bin_edges = np.array(n_binsx, dtype=float)
        ny_bin_edges = np.array(n_binsy, dtype=float)
        h2 = ROOT.TH2F(internal_name, "",
                       len(nx_bin_edges)-1, nx_bin_edges,
                       len(ny_bin_edges)-1, ny_bin_edges)
    else:
        raise TypeError(
            f"n_binsx and n_binsy must either both be `int` or `list`.\n"
            f"You provided type `{type(n_binsx)}`"
            )
    h2.Sumw2()
    # bin_wx = (x_max - x_min) / float(n_binsx)
    # ylabel = r"Events / (%s)" % bin_w
    x_label_withunits = x_label
    y_label_withunits = y_label
    if x_units is not None:
        x_label_withunits += r" (%s)" % x_units
    if y_units is not None:
        y_label_withunits += r" (%s)" % y_units
    h2.SetTitle(f"{y_label} vs. {x_label}") if title is None else h2.SetTitle(title)
    h2.SetXTitle(x_label_withunits)
    h2.SetYTitle(y_label_withunits)
    if z_min is not None and z_max is not None:
        h2.GetZaxis().SetRangeUser(z_min, z_max)
    if z_label_size is not None:
        h2.GetZaxis().SetLabelSize(z_label_size)
    h2.SetContour(n_contour)
    return h2

def fill_TH2F(h2, x_vals, y_vals, z_vals, zerr_vals=None):
    """Fill a 2-D hist with (x,y) entries (coordinates).
    
    FIXME:
    [ ] Think this function through to decide if necessary.
    - Seems like a dict {(x,y) : val} would be a more intuitive way to fill.
    
    Parameters
    ----------
    h2 : ROOT.TH2F
        The 2-D hist to be filled with len(x_coord_ls) entries.
    x_vals : list or array-like
        The x-coordinate of each entry.
        The first element pairs with the first in y_vals.
    y_vals : list or array-like
        The y-coordinate of each entry.
        The first element pairs with the first in x_vals.
    z_vals : list or array-like
        The number to be displayed in cell (x,y).
    """
    assert len(x_vals) == len(y_vals)
    if zerr_vals is not None:
        # Make a TH2 of errors.
        h2_err = h2.Clone()
    for x, y, z in zip(x_vals, y_vals, z_vals):
        h2.Fill(x, y, z)

def set_TH2F_errs(h2, h2_err):
    """Set the values stored in h2_err as the errors in h2."""
    assert h2.GetNbinsX() == h2_err.GetNbinsX()
    assert h2.GetNbinsY() == h2_err.GetNbinsY()
    for binx in range(1, h2.GetNbinsX() + 1):
        for biny in range(1, h2.GetNbinsY() + 1):
            # Similar cells (e.g. (2,3)) in both hists correspond to each other.
            err = h2_err.GetBinContent(binx, biny)
            h2.SetBinError(binx, biny, err)
#--------------#
#--- TGraph ---#
#--------------#
def make_TGraphErrors(x, y, x_err=None, y_err=None,
                      use_binwidth_xerrs=False, title="",
                      x_label="", x_min=0, x_max=10, x_units=None,
                      y_label="", y_min=0, y_max=10, y_units=None,
                      line_color=1, line_width=2,
                      marker_color=1, marker_style=21, marker_size=0.5):
    """A function to quickly make a TGraphErrors. Return the TGraphErrors.
    
    NOTE: Don't pass in x_err yet. Not ready.
    
    Parameters
    ----------
    x : list or array-like
        The x values to be plotted.
    y : list or array-like
        The y values to be plotted.
    """
    assert x_err is None, "FIXME: cannot implement `x_err` yet."
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
        return label + r" (%s)" % units if units is not None else label
    x_label = add_units(x_label, x_units)
    y_label = add_units(y_label, y_units)
    gr.GetXaxis().SetTitle(x_label)
    gr.GetYaxis().SetTitle(y_label)
    gr.GetXaxis().SetLimits(x_min, x_max)  # SetRangeUser() doesn't work for x-axis! SetLimits() instead.
    gr.GetYaxis().SetRangeUser(y_min, y_max)
    return gr

def get_coord_ls_from_TGraph(gr, axis):
    """Return a list of x- or y-coordinates of gr.
    
    Parameters
    ----------
    gr : TGraph
    axis : str
        "x" or "y"
    """
    n_pts = gr.GetN()
    if axis.lower() in "x":
        pt_ls = [gr.GetPointX(pt) for pt in range(n_pts)]
    elif axis.lower() in "y":
        pt_ls = [gr.GetPointY(pt) for pt in range(n_pts)]
    else:
        raise ValueError(f"Axis type ({type(axis)}) should be in ['x', 'y'].")
    return pt_ls

def get_xcoord_ls_from_TGraph(gr):
    """Return a list of x-coordinates of gr (TGraph)."""
    n_pts = gr.GetN()
    return [gr.GetPointX(pt) for pt in range(n_pts)]

#---------------#
#--- TLegend ---#
#---------------#
def make_TLegend(x_dim=(0.7, 0.9), y_dim=(0.7, 0.9), screenshot_dim=None, buffer_dim=None):
    """Return a TLegend, given its x and y coord.
    
    NOTE:
    - screenshot_dim and buffer_dim take precedence over x_dim and y_dim.
    If you want to use x_dim and y_dim, then set screenshot_dim and buffer_dim
    to `None`.

    - Some good numbers for screenshot_dim and buffer_dim:
    screenshot_dim=(878,872), buffer_dim=(176,438,610,131)

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

def make_TMultiGraph_and_Legend(gr_ls=[], leg_txt_ls=[], y_min=None, y_max=None,
                                x_dim=(0.7, 0.9), y_dim=(0.7, 0.9),
                                screenshot_dim=None, buffer_dim=None):
    """Return a (TMultiGraph, TLegend) with all the TGraphs from gr_ls added.
    
    NOTE:
    - screenshot_dim and buffer_dim take precedence over x_dim and y_dim.
    If you want to use x_dim and y_dim, then set screenshot_dim and buffer_dim
    to `None`.

    Parameters
    ----------
    gr_ls : list
        All the TGraphs to be put into the TMultiGraph.
    leg_txt_ls : list
        All the legend texts to be associated with the TGraphs.
        The order should coincide with the elements of gr_ls.
    y_min : float
        Set the minimum of the y-axis.
    y_max : float
        Set the maximum of the y-axis.
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
    mg = ROOT.TMultiGraph()
    for gr in gr_ls:
        mg.Add(gr, "p")
    x_title = gr_ls[0].GetXaxis().GetTitle()
    y_title = gr_ls[0].GetYaxis().GetTitle()
    title = gr_ls[0].GetTitle()
    all_titles = "%s;%s;%s" % (title, x_title, y_title)
    mg.SetTitle(all_titles)

    leg = make_TLegend(x_dim=x_dim, y_dim=y_dim,
                    screenshot_dim=screenshot_dim, buffer_dim=buffer_dim)
    # leg_text = r"%.2f < #deltam_{4#mu}/m_{4#mu} < %.2f%%" % (relmin, relmax)
    leg.SetTextSize(0.02)
    for gr, txt in zip(gr_ls, leg_txt_ls):
        leg.AddEntry(gr, txt, "lpfe")
    # leg.SetLineWidth(3)
    leg.SetBorderSize(1)
    if y_min is not None:
        mg.SetMinimum(y_min) # Change y-axis limits.
    if y_max is not None:
        mg.SetMaximum(y_max)
    return (mg, leg)

#-------------#
#--- TTree ---#
#-------------#
def add_branch_to_TTree(tree, br):
    """Add a branch (`br`) to `tree` and return the pointer to `br`.
    
    NOTE:
    - Only works with floats/doubles.
    - How to add a value to the TTree:
        ptr[0] = some_val
        tree.Fill()

    - If saving the TTree to '.root', remember to open the new TFile first!

    Parameters
    ----------
    tree : ROOT.TTree
    br : str
    """
    ptr = array('f', [0.])
    tree.Branch(br, ptr, "%s/F" % br)
    return ptr

def add_branch_pTcorr_fromadhocmethod(tree, pT_corr_factor_dict, outpath_root):
    """
    FIXME: Is this function done?
    
    Add a branch to a clone of a TTree.
    The branch is filled with pT corrected values, using ad hoc correction
    factors from pT_corr_factor_dict.
    The cloned tree is written to a new '.root' file (outpath_root).
    """
    outfile_root = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/"
    newfile = ROOT.TFile(outfile_path, "recreate")
    treeclone = tree.CloneTree(0)
    ptr = array('f', [0.])
    treeclone.Branch("pTcorr_adhoc", ptr, "m4mu/F")

def get_normcoord_from_screenshot(canv_width, canv_height,
                                  left_offset=None, right_offset=None, bot_offset=None, top_offset=None,
                                  obj_width=None, obj_height=None):
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
    obj_width : float
        If specified, will override any val passed into right_offset.
    obj_height : float
        If specified, will override any val passed into top_offset.

    Returns
    -------
    (x_min, x_max, y_min, y_max)
    """
    # A fraction of the total width.
    x_min = left_offset / float(canv_width)
    if obj_width is not None:
        msg = "Specify either `obj_width` or `right_offset`."
        assert right_offset is None, msg
        x_max = x_min + float(obj_width / canv_width)
    else:
        x_max = 1 - (right_offset / float(canv_width))
    # y dimension:
    y_min = bot_offset / float(canv_height)
    if obj_height is not None:
        msg = "Specify either `obj_height` or `top_offset`."
        assert top_offset is None, msg
        y_max = y_min + float(obj_height / canv_height)
    else:
        y_max = 1 - (top_offset / float(canv_height))
    return (x_min, x_max, y_min, y_max)