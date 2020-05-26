import os
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.backends.backend_pdf import PdfPages
from PyUtils.Utils_Physics import perc_diff

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
    
    mod_data = account_for_underoverflow_entries(data, x_limits[0], x_limits[1], x_bins)
    
    if y_max > 0: ax.set_ylim([0,y_max])
    if (log_scale): ax.set_yscale('log')
        
    stats = get_stats_1Dhist(data)
    label_legend = make_stats_legend_for_1dhist(stats)
    bin_vals, bin_edges, _ = ax.hist(mod_data, bins=x_bins, label=label_legend, histtype='step', color='b')
    ax.legend(loc='upper right', fontsize=7)
    return
  
def account_for_underoverflow_entries(data, x_min, x_max, bin_edges):
    """
    Add the entries in the underflow and overflow bins, 
    which would normally not be plotted on a histogram.
    Underflow entries will be added to first bin and 
    overflow entries are added to last bin, within the x_limits. 

    NOTES: 
        - It is not a good idea to run statistics on mod_data
          since it has been modified.
        - Under-overflow entries are only accounted for when
          'zooming' in on the x-axis, within the given bin_edges.

    Parameters
    ----------
    data : array-like
        The data to be binned in a histogram. 
    x_min : int or float
        The minimum value to be shown on the x-axis. 
    x_max : int or float
        The maximum value to be shown on the x-axis. 
    bin_edges : list or array-like
        [left_edge_first_bin, ..., right_edge_last_bin]

    Returns
    -------
    mod_data : array-like
        The modified data in which data values < x_min are converted to x_min and
        data_values > x_max are converted to x_max.
    """
    last_bin_width = bin_edges[-1] - bin_edges[-2]
    mod_data = data

    if x_min > bin_edges[0]:
        # Convert any value < x_min to x_min (underflow entries).
        mod_data = np.clip(mod_data, x_min, max(mod_data))
    if x_max < bin_edges[-1]:
        # Convert any value > x_max to x_max (overflow entries).

        # Overflow entries will not be seen on plot,
        # since x_max ends the last bin shown on plot.
        # However, clipped values are in the next bin (shown off screen).
        # So put overflow entries into last bin shown on plot by going
        # back 1 bin.
        mod_data = np.clip(mod_data, min(mod_data), x_max - last_bin_width)

    return mod_data

def print_header_message(msg, pad_char="-", n_center_pad_chars=5):
    pad_char_long = pad_char * (len(msg) + n_center_pad_chars*2 + 2)
    pad_char_short = pad_char * n_center_pad_chars
    banner = "#{}#".format(pad_char_long)
    middle = "#{} {} {}#".format(pad_char_short, msg, pad_char_short)
    print(banner)
    print(middle)
    print(banner)
    
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

def centers_of_binning_array(bin_arr, decimals=0):
    """
    Returns an array of the centers between adjacent bin edges.
    Also works with bin_arr that has unequal bin widths. 
    
    Example: 
        bin_arr = np.array([1, 2, 4, 8])
        bin_arr_shifted = np.array([1.5, 3, 6])
        
    NOTE: Was formerly called shift_binning_array().
    
    Parameters
    ----------
    bin_arr : list or array-like
        Edges of bins: [left_edge_first_bin, ..., right_edge_last_bin]
    decimals : int, optional
        Number of decimal places to round bin centers.
        
    Returns
    -------
    bin_centers_arr : array
        Centers between adjacent bin edges.
        Will be of length = len(bin_arr) - 1
    """
    bin_arr = np.array(bin_arr)
    bin_edges_low = bin_arr[:-1]
    bin_edges_high = bin_arr[1:]
    
    bin_centers_arr = (bin_edges_high + bin_edges_low) / 2.
    if decimals != 0:
        bin_centers_arr = np.round(bin_centers_arr, decimals=decimals)
    
    return bin_centers_arr

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

def integrate_hist_over_region(bin_edges, bin_vals):
    """
    Calculate the total area of a histogram over the specified bin_edges.
    Performs a simple width * height calculation for each bin area, then sums up all areas. 
    
    Parameters
    ----------
    bin_edges : array
        The boundaries of all the bins. They can be variable bin widths.
        E.g., np.array([first_bin_left_edge, first_bin_right_edge, ..., last_bin_left_edge, last_bin_right_edge])
    bin_vals : array
        An array which contains the entries in each bin.
        len(bin_vals) = len(bin_edges) - 1
        
    Returns
    -------
    total_area : float
        Usually the total entries of all bins being summed.
    first_bin_left_edge : float
        Where the summing started along the x-axis, in x-axis units.
    last_bin_right_edge : float
        Where the summing ended along the x-axis, in x-axis units.
    """
    first_bin_left_edge = bin_edges[0]
    last_bin_right_edge = bin_edges[-1]
    
    # Area calculation
    bin_width_arr = calc_bin_widths(bin_edges)
    total_area = np.sum(bin_width_arr * bin_vals)
    
    return total_area, first_bin_left_edge, last_bin_right_edge

def find_equal_hist_divisions(bin_edges, bin_vals, K, verbose=False):
    """
    Return the bin edges which divide the histogram into divisions with equal entries per division. 
    
    Algorithm example:
        Want to split histogram of 100 entries up into 3 divisions. 
        Calculate hypothetical entries_per_division: 100/3 = 33.333.
        Start from first bin and keep adding entries in each bin until you exceed entries_per_division. 
        Record the right_end_bin_edge where this happens.
        Check to see if maybe stopping at the previous bin would have been a better decision
        (this is done by a percent difference calculation).
        Do only 2 scans over the first two divisions.
        Third division is found by deduction.
    Remember that entries is defined as bin_width * bin_height!

    Parameters
    ----------
    bin_edges : array
        Edges of each bin, starting from the left-most edge to the right-most edge. 
            bin_edges[0] = left edge of first bin
            bin_edges[-1] = right edge of last bin
    bin_vals : array
        "Height" of each bin (NOT necessarily the entries in each bin). 
        Again: bin_width * bin_height = bin_entries
        
    K : int
        Number of divisions to split histogram.
        Each division will contain approximately the same number of entries. 
        
    Returns
    -------
    bin_div_ls : list
        A list of the edges of each division along the x-axis. 
        The number of entries between divisions should have ~same number of entries.
    bin_stats_dict : dict
        Contains the better percent difference of the number of entries found in each division
        compared to the ideal number of entries.
    Notes:
    Should be the case that: len(bin_edges) == len(bin_vals) + 1
    """
    # TOTAL ENTRIES IS SUM(bin_width * bin_height)
    entries_total, first_bin, last_bin = integrate_hist_over_region(bin_edges, bin_vals)
    entries_per_div = entries_total / float(K)
    bin_div_ls = [first_bin]  # This is what we are after ultimately. 
    bin_stats_dict = {}
    
    if (verbose): 
        print("[INFO] Performing {} divisions on {} total entries.".format(K, entries_total))
    
    start_elem = 0
    # Always do K - 1 loops, except when K = 1. 
    div_ls = [1] if K == 1 else list(range(1,K))
    
    tmp_division_entries = 0
    for count in div_ls:
        if (verbose): print("[INFO]  Beginning division {}...".format(count))
        
        end_elem = start_elem + 1
        while True:
            # Start at the beginning of this division.
            # Go bin by bin until you reach or surpass entries_per_div.
            division_bin_edges = bin_edges[start_elem:end_elem+1]  # Extra 1 because of Python's exclusion.
            division_bin_vals = bin_vals[start_elem:end_elem]
            division_entries, first_bin_left_edge, last_bin_right_edge = integrate_hist_over_region(division_bin_edges, division_bin_vals)
            if division_entries >= entries_per_div: 
                break
            else:
                end_elem += 1
                
        # Check if previous bin would have been a better choice. 
        prev_division_bin_edges = bin_edges[start_elem:end_elem]
        prev_division_bin_vals = bin_vals[start_elem:end_elem-1]
        prev_division_entries, prev_first_bin_left_edge, prev_last_bin_right_edge = integrate_hist_over_region(prev_division_bin_edges, prev_division_bin_vals)
        if (verbose):
            print("[INFO]    Integration complete:")
            print("[INFO]      division_bin_edges:",division_bin_edges)
            print("[INFO]      division_bin_vals:",division_bin_vals)
            print("[INFO]    Found {:.2f} entries when looking for {:.2f}. Checking previous bin...".format(division_entries, entries_per_div))
            print("[INFO]    Found {:.2f} entries when looking for {:.2f}.".format(prev_division_entries, entries_per_div))

        div_perc_diff = perc_diff(division_entries, entries_per_div)
        prev_div_perc_diff = perc_diff(prev_division_entries, entries_per_div)
        if abs(div_perc_diff) <= abs(prev_div_perc_diff):
            best_div_entries = division_entries
            best_perc_diff = div_perc_diff
            best_bin_edge = last_bin_right_edge
            best_div_bin_edges = division_bin_edges
            best_div_bin_vals = division_bin_vals
            # Keep end_elem where it is.
        else:
            best_div_entries = prev_division_entries
            best_perc_diff = prev_div_perc_diff
            best_bin_edge = prev_last_bin_right_edge
            best_div_bin_edges = prev_division_bin_edges
            best_div_bin_vals = prev_division_bin_vals
            # Need to decrement end_elem. 
            end_elem -= 1
                
        if (verbose): 
            print("[INFO]    best_div_entries:   {:.2f}".format(best_div_entries)) 
            print("[INFO]    best_perc_diff:     {:.4f}%".format( best_perc_diff))
            print("[INFO]    best_bin_edge:      {}".format( best_bin_edge))
            print("[INFO]    best_div_bin_edges: {}".format( best_div_bin_edges))
            print("[INFO]    best_div_bin_vals:  {}".format( best_div_bin_vals), "\n")
                
        # Found the best element, sum, and perc diff which correspond to this division bin. 
        bin_stats_dict["Division{}".format(count)] = best_perc_diff
        bin_div = bin_edges[end_elem]  # Looks right.
        bin_div_ls.append(bin_div)
        # Start the next division at the element that we stopped at last. 
        # So don't change the counting of elem.
        start_elem = end_elem
        
        tmp_division_entries += best_div_entries
    # All divisions, except last, performed successfully.
    # Since last division isn't performed, append final bin.
    if K != 1:
        bin_div_ls.append(bin_edges[-1])
    
    if len(set(bin_div_ls)) != len(bin_div_ls):
        err_msg = "[ERROR] The same bin edge was found multiple times for some reason.\n"
        err_msg += "Most likely the value of K ({}) was too large. Try fewer divisions.".format(K)
        raise RuntimeError(err_msg)

    if K != 1: 
    # Deduce last division. 
        final_bin_edges = bin_edges[end_elem:]
        excess_entries, first_bin_here, last_bin_here = integrate_hist_over_region(final_bin_edges, bin_vals[end_elem:])
        perc_diff_excess = perc_diff(excess_entries, entries_per_div)
        tmp_division_entries += excess_entries
        bin_stats_dict["Division{}".format(K)] = perc_diff_excess
        if (verbose):
            print("[INFO]  Division {} was analyzed by deduction...".format(K)) 
            print("[INFO]    first_bin_last_division={}, last_bin_last_division={}".format(first_bin_here,last_bin_here))
            print("[INFO]    final bin edges: {}".format(final_bin_edges))
            print("[INFO]    Excess entries in this division: {} ({:.4f}% different from number searched for)".format(excess_entries, perc_diff_excess))
            
    if (verbose):
        print("[INFO]    Total entries found after full analysis: {}".format(tmp_division_entries))
        print("[INFO]    Final bin division list is:\n{}".format(bin_div_ls))

    return bin_div_ls, bin_stats_dict

# def find_equal_hist_regions_unbinned(vals_arr, r, round_to_n_decimals=2, algo=("normal", -1), verbose=False):
def find_equal_hist_regions_unbinned(vals_arr, r, algo=("normal", -1), verbose=False):
    """
    Return the "bin edges" (really just specific values of vals_arr) 
    which divide the array into regions with equal number of values per region. 
    
    Algorithm example:
        Want to split values of 100 entries up into 3 regions. 
        Calculate hypothetical entries_per_region: 100/3 = 33.333.
        Start from first val and simply count until you exceed entries_per_region. 
        Record the right_end_bin_edge where this happens.
        It's actually a little more complicated than that, 
        but you get the gist.

    Parameters
    ----------
    vals_arr : array
        Array of values which will be sorted and then split into `r` regions.
    r : int
        Number of regions to split histogram into.
        Each region will contain approximately the same number of entries. 
    DELETEround_to_n_decimals : int, optional
        Number of decimals to round bin edge values to. 
    algo : 2-tuple, optional
        algo[0] : str
            The kind of algorithm to run as described below.
        algo[1] : int 
            The minimum number of entries per region acceptable. 
            If `-1` then use the `r` specified.
            
        Possible algo's:
        "normal" : default, split vals_arr up into `r` regions. 
        "at_least" : request at least algo[1] entries per region, first trying `r` regions. 
                     If there are fewer than algo[1] entries per region, then try `r-1` regions.
                     Try recursively until to a minimum of 2 regions. Stops at 2 regions. 

    verbose : bool, optional
        Gives lots of juicy debugging details.
        
    Returns
    -------
    bin_reg_ls : list
        A list of the edges of each region along the x-axis. 
        The number of entries between regions should have ~same number of entries.
        If 
    r : int
        A possibly updated number of regions, if algo was set to something other than "normal".
        
    Notes:
        - I used to call them "divisions" but now I call them "regions"
    """
    entries_total = len(vals_arr)
    if entries_total == 0:
        return [None, None] , -1
    sorted_vals_arr = np.sort(vals_arr)
    # first_bin = round(sorted_vals_arr[0], round_to_n_decimals)
    # last_bin = round(sorted_vals_arr[-1], round_to_n_decimals)
    first_bin = sorted_vals_arr[0]
    last_bin = sorted_vals_arr[-1]
    bin_reg_ls = [first_bin]  # This is what we are after ultimately. 

    # Prepare the expectation of each region, based on sorted_vals_arr.
    entries_per_reg = float(entries_total) / float(r)
    
    mode = algo[0]
    find_at_least_per_reg = algo[1]
    if (mode in "at_least"):
        # Make sure User specified desired number of entries per region.
        assert find_at_least_per_reg > 0

        while entries_per_reg < find_at_least_per_reg:
            if (r == 2): 
                msg = "[WARNING] Could not find at least {} entries per region.".format(find_at_least_per_reg)
                print_header_message(msg)
                break
            # There are too few actual entries per region. 
            # Decrement the number of regions.  
            msg  = "  Found {:.2f} entries per region (using {} regions), ".format(entries_per_reg, r)
            msg += "but require at least {} entries per region.\n".format(find_at_least_per_reg)
            msg += "    Decrementing the number of regions from {} to {}".format(r, r-1)
            r -= 1
            print(msg)
            
            entries_per_reg = float(entries_total) / float(r)
        print("Splitting array into {} equal-entry regions, {:.2f} entries per region.".format(r, entries_per_reg))

    entries_per_reg_roundup = math.ceil(entries_per_reg)
    entries_per_reg_rounddown = math.floor(entries_per_reg)
    
    num_regions_with_more = entries_total % r
    num_regions_with_fewer = r - num_regions_with_more
    # Example: 
    #--- sorted_vals_arr = [1,1,2,2,3]
    #--- N=5, say r=3 
    #--- entries_per_reg=1.667, entries_per_reg_roundup=2, entries_per_reg_rounddown=1
    #--- num_regions_with_more=2, num_regions_with_fewer=1
    #--- make 3 regions: [1,1 | 2,2 | 3]  
    #--- two regions have 2 entries and 1 region has 1 entry
    
    if (verbose):
        print("entries_total:             {}".format(entries_total))
        print("num_regions_with_more:     {}".format(num_regions_with_more))
        print("num_regions_with_fewer:    {}".format(num_regions_with_fewer))
        print("first_bin_edge:            {:.4f}".format(first_bin))
        print("last_bin_edge:             {:.4f}".format(last_bin))
        print("entries_per_reg:           {:.2f}".format(entries_per_reg))
        print("entries_per_reg_roundup:   {:.2f}".format(entries_per_reg_roundup))
        print("entries_per_reg_rounddown: {:.2f}".format(entries_per_reg_rounddown))
        print("[INFO] Making {} regions on {} total entries.".format(r, entries_total))
        print("[INFO] Looking for {:.2f} entries per region.".format(entries_per_reg))
        
    def scan_arr_get_index(arr, start_elem, entries_to_scan):
        """
        Scan an array over a specified number of sequential elements
        and return the "bin edge" corresponding to the final element. 
        
        Parameters
        ----------
        arr : list or array-like
            Array to scan over. 
        start_elem : int
            The element in the array at which to start scanning. 
        entries_to_scan : int
            Number of sequential entries in array to begin scanning.
        
        Returns
        -------
        end_elem : int
            The final element in arr on which scanning has ended. 
            This is the element corresponding to bin_edge.
        bin_edge : anything
            The value corresponding to arr[end_elem].
        """
        end_elem = start_elem + entries_to_scan-1
        bin_edge = arr[end_elem]
        # bin_edge = round(bin_edge, round_to_n_decimals)
        return end_elem, bin_edge

    # Did brief testing and concluded that it doesn't matter 
    # whether you loop over the regions with more entries
    # or fewer entries first. There are so many statistics generally that 
    # it doesn't matter. 
    elem = 0
    for _ in range(num_regions_with_fewer):
        elem, bin_edge = scan_arr_get_index(sorted_vals_arr, elem, entries_per_reg_rounddown)
        bin_reg_ls.append(bin_edge)
        elem += 1  # Must not include this element again. 
    for _ in range(num_regions_with_more):
        elem, bin_edge = scan_arr_get_index(sorted_vals_arr, elem, entries_per_reg_roundup)
        bin_reg_ls.append(bin_edge)
        elem += 1  # Must not include this element again. 
              
    # I think this will only trigger if r == entries_total.
    # Can happen when the number of entries is VERY small. 
    if len(set(bin_reg_ls)) != len(bin_reg_ls):
        from collections import Counter
        c = Counter(bin_reg_ls)
        print(c)
        multiple = c.most_common(1)[0][0]
        msg  = "[WARNING] The same bin edge ({multiple}) was found multiple times.\n"
        msg += "Using an ad hoc solution: duplicating first bin and shifting it left a little."
        print(msg)
        bin_reg_ls[0] = bin_reg_ls[0] - 1E-14 * abs(bin_reg_ls[0])
        # raise RuntimeError(err_msg)

    if (verbose):
        print("[INFO] Final bin region list is:\n{}\n".format(bin_reg_ls))

    return bin_reg_ls, r

def collapse_eta_bin_edges(bin_ls, round_to_n_decimals=2):
    """
    Return a list of bin edges that where symmetric elements are averaged.  
    
    Concept:
        Originally, an eta bin edge list looks something like:
            bin_ls = [-2.4, -0.89, 0.001, 0.91, 2.4]
        But since eta is binned using abs(eta) we need to collapse
        this list to become:
            new_bin_ls = [0.00, 0.90, 2.4]
    
    Notes:
        len(new_arr) == len(orig_arr) / 2, rounded up.

    Parameters
    ----------
    bin_ls : list or array-like
        The bin edges to be combined.
    round_to_n_decimals : int
        Number of decimals to round bin edge values to. 
        
    Returns
    -------
    rev_ls : list
        A sorted list (ascending order) of the averaged, symmetric bin edges.
    """
    n = len(bin_ls)
    mid_elem = math.floor(n / 2)
    if (n % 2) == 0:
        # Stop when passing the midpoint of an even-length list.
        # We would go too far, so decrement. 
        mid_elem -= 1
    
    new_ls = []
    elem = 0
    while elem <= mid_elem:
        elem_mirror = -1*elem - 1
        this_num = bin_ls[elem]
        mirror_num = bin_ls[elem_mirror]

        tmp_arr = np.array([this_num, mirror_num])

        avg = np.mean(np.abs(tmp_arr))
        new_ls.append(round(avg, round_to_n_decimals))
        elem += 1
    rev_ls = list(reversed(new_ls))
        
    return rev_ls