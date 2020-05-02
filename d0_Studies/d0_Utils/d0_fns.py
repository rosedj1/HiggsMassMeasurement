import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.backends.backend_pdf import PdfPages

#from d0_Utils.d0_cls import HistInfo    # Not sure why I can't import this...

def get_subset_mask(x_vals, x_min, x_max):
    """
    Return a mask of an array such that: x_min <= element <= x_max
    Useful for selecting fit ranges along the 'x-axis'.
    
    Apply the mask: x_vals[mask]
    
    Parameters
    ----------
    x_vals : array
        The original array from which the subset will be made.
    x_min : float
        The minimum value in x_vals to be kept.
    x_max : float
        The maximum value in x_vals to be kept.
        
    Returns
    -------
    mask : bool array
        An array of booleans. If an element is 'True' then: x_min <= element <= x_max.
        The mask will be the same length as x_vals. 
    """
    mask = (x_min <= x_vals) & (x_vals <= x_max)
    
    return mask

def combine_cut_list_ROOT(cut_list):
    """Used for ROOT-style cuts while doing tree.Draw()"""
    total_cuts = ''
    for cut in cut_list:
        total_cuts += cut + ' && '
    cut_str = total_cuts.rstrip(' && ')
    return cut_str

def combine_cut_list(cut_list):
    """Take the cuts from a list of cuts and string them together into a single string."""
    total_cuts = ''
    for cut in cut_list:
        total_cuts += cut + '\n'
    cut_str = total_cuts.rstrip('\n')
    return cut_str

def calc_num_bins(bin_min, bin_max, bin_width):
    """
    Calculates the number of bins (int) given: bin_min, bin_max, bin_width

    Parameters
    ----------
    bin_min, bin_max, bin_width
    """
    return int(round( (bin_max-bin_min)/bin_width ))

def calc_bin_widths(bin_edges):
    """
    Return an array which contains the bin widths of each bin. 
    Bins can be of unequal widths. 
    
    Parameters
    ----------
    bin_edges : array
        The boundaries of all the bins. They can be variable bin widths.
        E.g., np.array([first_bin_left_edge, first_bin_right_edge, ..., last_bin_left_edge, last_bin_right_edge])
    
    Returns
    -------
    bin_width_arr : array
        The width of each bin. len(bin_width_arr) = len(bin_edges) - 1
    """
    # Create left and right edges of bins and stagger them. Then take the difference.
    left_edges = bin_edges[:-1]
    right_edges = bin_edges[1:]
    bin_width_arr = right_edges - left_edges
    return bin_width_arr

def calc_x_err_bins(x_val_center_list):
    """
    Returns lists of the "x-errors", i.e. the half-distances between neighboring bins.
    These may be symmetrical or asymmetrical. 
    
    FIXME?: Possibly incorporate calc_bin_widths() in here?
    """
    high_err_list = []
    for x_center in range(len(x_val_center_list)-1):
        high_err = float( x_val_center_list[x_center+1] - x_val_center_list[x_center] ) / 2
        high_err_list.append(high_err)
        
    low_err_list = [high_err_list[0]] + high_err_list
    high_err_list = high_err_list + [high_err_list[-1]]
    
    return low_err_list, high_err_list

def calc_ymin_for_legend(n_graphs, text_height=0.042):
    """
    ROOT specific! Allows for dynamic expansion of bottom axis of TLegend. 
    There is probably an in-built function, but I want to write my own!
    
    Parameters
    ----------
    n_obj : int
        Number of objects (like graphs, histos) that will show up in your legend.
    text_height : float between 0 and 1
        Essentially sets the font. Was 0.033333 but you'll have to play with it. 
    """
    delta_y = text_height*n_graphs
    return 0.9 - delta_y # When you do TCanvas(), it puts the top axis at y = 0.9.

            
def make_kinem_subplot(lep, ax, data, x_limits, x_bins, x_label, y_label, y_max=-1, log_scale=False):
    """
    MAY BE DEPRECATED.
    Draw a kinematic distribution (e.g. eta1, gen_phi2, etc.) to a subplot in a figure.
    
    Parameters
    ----------
    lep : int
        Either `1` or `2`. Indicates which lepton you are referring to. 
    ax : axis object
        The axes to which plot will be drawn.
    data : array-like
        Data to be histogrammed. 
    x_limits : 2-element list
        A list of the x-min and x-max to be plotted.
    x_bins : array-like 
        Array of bin edges. Should be of length = len(data)+1. <-- I don't think that's true...
    x_label : str
        Label for x-axis.
    y_label : str
        Label for y-axis.
    label_legend : str
        Label for legend.
    y_max : float
        Max on y-axis. If y_max <= 0, then matplotlib will choose y_max.
    log_scale : bool
        If True, sets the y-axis to log scale. 
    """
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    ax.set_xlim(x_limits)
    ax.grid(False)
    
    mod_data, n_underflow, n_overflow = add_underoverflow_entries(data, x_limits[0], x_limits[1])
    
    if y_max > 0: ax.set_ylim([0,y_max])
    if (log_scale): ax.set_yscale('log')
        
    stats = get_stats_1Dhist(data)
    label_legend = make_stats_legend_for_1dhist(stats)
    bin_vals, bin_edges, _ = ax.hist(mod_data, bins=x_bins, label=label_legend, histtype='step', color='b')
    ax.legend(loc='upper right', fontsize=7)
    return
  
# def make_1D_dist(ax, data, x_limits, x_bins, x_label, y_label, title, y_max=-1, log_scale=False):
#    """
#    Draw a kinematic distribution (e.g. eta1, gen_phi2, etc.) to an axes object in a figure.
#    
#    Parameters
#    ----------
#    ax : axes object
#        The axes to which plot will be drawn.
#    data : array-like
#        Data to be histogrammed. 
#    x_limits : 2-element list
#        A list of the x-axis range to be plotted: 
#        x_limits=[x-min, x-max]
#    x_bins : array-like 
#        Array of bin edges. Should be of length = len(n_bins)+1.
#    x_label : str
#        Label for x-axis.
#    y_label : str
#        Label for y-axis.
#    title : str
#        Label for title.
#    y_max : float
#        Max on y-axis. If y_max <= 0, then matplotlib will choose y_max.
#    log_scale : bool
#        If True, sets the y-axis to log scale. 
#    """
#    textsize_legend = 9
#    textsize_axislabels = 12
#    textsize_title = 12
#            
#    ax.set_xlabel(x_label, fontsize=textsize_axislabels)
#    ax.set_ylabel(y_label, fontsize=textsize_axislabels)
#    ax.set_title(title, fontsize=textsize_title)
#    
#    ax.set_xlim(x_limits)
#    ax.grid(False)
#    
#    mod_data, n_underflow, n_overflow = add_underoverflow_entries(data, x_limits[0], x_limits[1])
#    
#    if y_max > 0: ax.set_ylim([0,y_max])
#    if (log_scale): ax.set_yscale('log')
#        
#    stats = get_stats_1Dhist(data)
#    label_legend = make_stats_legend_for_1dhist(stats)
#    bin_vals, bin_edges, _ = ax.hist(mod_data, bins=x_bins, label=label_legend, histtype='step', color='b')
#    ax.legend(loc='upper right', framealpha=0.9, fontsize=textsize_legend)
#    
#    return ax, stats
    
def print_header_message(msg):
    n = len(msg)
    octothorpes = (n+12)*'#'  # Who names their variable 'octothorpes'? Honestly?
    buff = 5*'#'
    print(octothorpes)
    print(buff, msg, buff)
    print(octothorpes)
    
def make_binning_array(lim_ls):
    """
    Turn a list of [bin_min, bin_max, bin_width] into an array of length = (bin_max - bin_min)/bin_width,
    whose first element is bin_min, last element is bin_max and has spacing bin_width.
    
    Returns
    -------
    bin_edges : array
        Contains equally-spaced bin edges, where:
            bin_edges[0] = left edge of first bin
            bin_edges[-1] = right edge of last bin
    binw : float
        The width of a bin. All are equal width.
    """
    ls_min = lim_ls[0]
    ls_max = lim_ls[1]
    binw = lim_ls[2]
    
    bin_edges = np.arange(ls_min, ls_max+0.5*binw, binw)
    
    return bin_edges, binw

def shift_binning_array(bin_arr):
    """
    Returns an array of the centers of the bins of bin_arr, excluding the very last bin.
    
    Example: 
        bin_arr = np.array([1, 2, 3, 4])
        bin_arr_shifted = np.array([1.5, 2.5, 3.5])
        
    FIXME: Somewhere in my DarkZ code is this function which can accommodate an asymmetrical bin_arr.
    """
    bin_width = bin_arr[1] - bin_arr[0]
    
    return bin_arr[:-1] + 0.5 * bin_width

def event_counter(start, stop, step, df, branch):
    """
    I THOUGHT I WAS CLEVER WHEN WRITING THIS FUNCTION, BUT THE MUCH MORE EFFICIENT WAY IS:
    
    vals, bin_edges = np.histogram(df[branch], bins=np.arange(start, stop * step/2, step))
    return vals
    
    ORIGINAL:
    Counts the number of events in each bin.
    This function provides a way to make sure that your histogram is binning the correct number of events per bin.
    
    Parameters
    ----------
    start : float 
        The left edge of the very first bin to begin counting.
    stop : float 
        The right edge of the very last bin to stop counting. 
    step : float
        The bin width. 
    df : pandas.DataFrame
        The DataFrame that holds branch.
    branch : str
        The DataFrame column name that holds the data to be binned.
    """
    # vals = np.arange(start, stop, step)
    # for v in vals:
    #     print(f"n_evts with {v:.1f} < eta < {v+step : .1f}:", sum((v < df[branch]) & (df[branch] < v+step) ) )
    vals, bin_edges = np.histogram(df[branch], bins=np.arange(start, stop * step/2, step))
    return vals