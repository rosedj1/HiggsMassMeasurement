import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import sys
sys.path.append('/Users/Jake/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly

def change_cmap_bkg_to_white(colormap, n=256):
    """
    Make any matplotlib colormap have a white background (i.e. any bins with 0 entries are colored white).
    This makes it easy to distinguish between bins with 0 entries and bins with only a few entries.
    
    Parameters
    ----------
    colormap : Colormap or str, optional
        A `.colors.Colormap` instance.  If not set, use rc settings.
    n : int
        The number of divisions of the color bar. 
        Use a higher number for a smoother transition of colors.
    """
    colors = cm.get_cmap(colormap, n)
    newcolors = colors(np.linspace(0, 1, n))
    white = np.array([1, 1, 1, 1])    # White background (Red, Green, Blue, Alpha).
    newcolors[0, :] = white    # Only change bins with 0 entries.
    newcmp = ListedColormap(newcolors)
    
    return newcmp

def save_plots_to_outpath(save_plot=False, outpath="", file_name="DEFAULT_NAME", save_as_png=False, verbose=False):
    if (save_plot):
        file_name = make_str_title_friendly(file_name)
        
        fullpath = os.path.join(outpath, file_name)
        
        makeDirs(outpath)
        
        plt.savefig(fullpath + '.pdf')
        if (verbose): print(f"Figure saved at: {fullpath + '.pdf'}")
        if (save_as_png):
            plt.savefig(fullpath + '.png')
            if (verbose): print(f"Figure saved at: {fullpath + '.png'}")

def make_1D_dist(ax, data, x_limits, x_bins, x_label, y_label, title, y_max=-1, log_scale=False):
    """
    Draw a kinematic distribution (e.g. eta1, gen_phi2, etc.) to an axes object in a figure.
    
    Parameters
    ----------
    ax : axes object
        The axes to which plot will be drawn.
    data : array-like
        Data to be histogrammed. 
    x_limits : 2-element list
        A list of the x-axis range to be plotted: 
        x_limits=[x-min, x-max]
    x_bins : array-like 
        Array of bin edges. Should be of length = len(n_bins)+1.
    x_label : str
        Label for x-axis.
    y_label : str
        Label for y-axis.
    title : str
        Label for title.
    y_max : float
        Max on y-axis. If y_max <= 0, then matplotlib will choose y_max.
    log_scale : bool
        If True, sets the y-axis to log scale. 
    """
    textsize_legend = 9
    textsize_axislabels = 12
    textsize_title = 12
            
    ax.set_xlabel(x_label, fontsize=textsize_axislabels)
    ax.set_ylabel(y_label, fontsize=textsize_axislabels)
    ax.set_title(title, fontsize=textsize_title)
    
    ax.set_xlim(x_limits)
    ax.grid(False)
    
    mod_data, n_underflow, n_overflow = add_underoverflow_entries(data, x_limits[0], x_limits[1])
    
    if y_max > 0: ax.set_ylim([0,y_max])
    if (log_scale): ax.set_yscale('log')
        
    stats = get_stats_1Dhist(data)
    label_legend = make_stats_legend_for_1dhist(stats)
    bin_vals, bin_edges, _ = ax.hist(mod_data, bins=x_bins, label=label_legend, histtype='step', color='b')
    ax.legend(loc='upper right', framealpha=0.9, fontsize=textsize_legend)
    
    return ax, bin_vals, bin_edges, stats
                                
def add_underoverflow_entries(data, xmin, xmax):
    """
    Add the entries in the underflow and overflow bins, which would normally not be plotted on a histogram.
    NB: Run all statistics on data BEFORE binning and before adding under-overflow entries.
    
    Parameters
    ----------
    data : array-like
        The data, usually a numpy array or DataFrame series, which will get binned in a histogram. 
    xmin : float
        The left edge of the first bin. Any values in data that are less than xmin
        (which usually would not be plotted) are now put into the xmin bin. 
    xmax : float
        The right edge of the last bin. Overflow values fit into this final bin. 
        
    Returns
    -------
    mod_data : array-like
        The modified data with fake underflow and overflow entries appended. 
        This way the binned mod_data will show under-overflow entries.
    n_underflow : int
        The number of entries in the underflow bin.
    n_overflow : int
        The number of entries in the overflow bin.
    """
    n_underflow = len(data[data < xmin])
    n_overflow = len(data[data > xmax])
    
    # Create new data array where entries that would fall outside x-range
    # now sit at underflow/overflow bins.
    fake_under = n_underflow * [xmin]
    fake_over = n_overflow * [xmax]

    mod_data = data[ (xmin <= data) & (data <= xmax) ]
    mod_data = np.append(mod_data, fake_under + fake_over)

    n_data_init = len(data)
    n_data_final = len(mod_data)

    # Make sure that things make sense:
    if n_data_init != n_data_final:
        err_msg = (f"The initial number of entries in data ({n_data_init}) is different from "
                   f"the final number ({n_data_final}), after accounting for underflow/overflow bins.\n"
                   "Stopping now.")
        raise RuntimeError(err_msg)
                                
    return mod_data, n_underflow, n_overflow

def get_stats_1Dhist(data):
    """
    Return the statistics of array-like data.
    Particularly good for displaying the statistics on binned data in a histogram legend. 

    Parameters
    ----------
    data : array-like
        Your data (in the form of a list, numpy array, series, etc.) 

    Returns
    -------
    stats : list
        A 5-element list of the statistics of data.
        stats[0] -> n : int
            Number of entries in data.
        stats[1] -> mean : float
            Unweighted average of data. 
        stats[2] -> mean_err : float
            Standard error on the mean. 
        stats[3] -> stdev : float
            Standard deviation of data.
        stats[4] -> stdev : float 
            Standard error on the stdev. 
    """
    n = len(data)
    mean = np.mean(data)
    mean_err = np.abs(mean) / np.sqrt(n)
    stdev = np.std(data)
    stdev_err = stdev / np.sqrt(2*n)
    
    stats = [n, mean, mean_err, stdev, stdev_err]

    return stats

def get_stats_2Dhist(x_data, y_data):
    """
    Return the statistics of 2D histogrammed, array-like data.
    Particularly good for displaying the statistics on binned data in a histogram legend. 

    FIXME: does not do any mean or stdev calculation yet.
    
    Parameters
    ----------
    data : array-like
        Your data (in the form of a list, numpy array, series, etc.) 

    Returns
    -------
    stats : list
        A 5-element list of the statistics of data.
        stats[0] -> n : int
            Number of entries in 2D hist. Note that this is not 
        stats[1] -> mean : float
            Unweighted average of data. 
        stats[2] -> mean_err : float
            Standard error on the mean. 
        stats[3] -> stdev : float
            Standard deviation of data.
        stats[4] -> stdev : float 
            Standard error on the stdev. 
    """
    if len(x_data) != len(y_data):
        raise ValueError("The len(x_data) != len(y_data). Check your data going into this 2D histogram.")
    n = len(x_data)
#     mean = np.mean(data)
#     mean_err = abs(mean) / np.sqrt(n)
#     stdev = np.std(data)
#     stdev_err = stdev / np.sqrt(2*n)
    
#     stats = [n, mean, mean_err, stdev, stdev_err]
    stats = [n]

    return stats

def make_stats_legend_for_1dhist(stats_ls):
    """
    Create a legend label that displays the statistics of a 1D histogram. 
    
    Parameters
    ----------
    stats_ls : list
        A 5-element list of statistics from the 1D data used to make the histogram. 
        
        stats_ls[0] -> number of entries in data
        stats_ls[1] -> mean of data
        stats_ls[2] -> standard error on the mean of data
        stats_ls[3] -> standard deviation of data
        stats_ls[4] -> standard error on the stdev of data
    
    Returns
    -------
    leg_label : str
        A string of the data statistics, useful for making a legend label. 
    """
    n = stats_ls[0]
    mean = stats_ls[1]
    mean_err = stats_ls[2]
    stdev = stats_ls[3]
    stdev_err = stats_ls[4]

    leg_label  = "Entries = {}".format(n) + "\n"
    leg_label += "Mean = {:.2E}".format(mean) + r" $\pm$ " + "{:.2E}".format(mean_err) + "\n"
    leg_label += "Std Dev = {:.2E}".format(stdev) + r" $\pm$ " + "{:.2E}".format(stdev_err)

    return leg_label

def make_stats_legend_for_2dhist(stats_ls):
    """
    Create a legend label that displays the statistics of a 2D histogram. 
    
    FIXME: does not show any mean or stdev values yet.
    
    Parameters
    ----------
    stats_ls : list
        A 5-element list of statistics from the 1D data used to make the histogram. 
        
        stats_ls[0] -> number of entries in data
        stats_ls[1] -> mean of data
        stats_ls[2] -> standard error on the mean of data
        stats_ls[3] -> standard deviation of data
        stats_ls[4] -> standard error on the stdev of data
    
    Returns
    -------
    leg_label : str
        A string of the data statistics, useful for making a legend label. 
    """
    n = stats_ls[0]
#     mean = stats_ls[1]
#     mean_err = stats_ls[2]
#     stdev = stats_ls[3]
#     stdev_err = stats_ls[4]

    leg_label  = "Entries = {}".format(n)
#    leg_label += "Mean = {:.2E}".format(mean) + r" $\pm$ " + "{:.2E}".format(mean_err) + "\n"
#    leg_label += "Std Dev = {:.2E}".format(stdev) + r" $\pm$ " + "{:.2E}".format(stdev_err)

    return leg_label
                                
def make_stats_legend_for_gaus_fit(popt, popt_err):
    """
    Create a legend label that displays the statistics after fitting a 1D histogram with a Gaussian curve.
    
    Parameters
    ----------
    
    Returns
    -------
    leg_label : str
        A string of the fit statistics, useful for making a legend label. 
    """
    coeff = popt[0]
    mean = popt[1]
    stdev = popt[2]
    
    coeff_err = popt_err[0]
    mean_err = popt_err[1]
    stdev_err = popt_err[2]
    
    leg_label  = r"Fit $C$ = {:4.3E} $\pm$ {:4.3E}".format(coeff, coeff_err) + "\n"
    leg_label += r"Fit $\mu$ = {:4.3E} $\pm$ {:4.3E}".format(mean, mean_err) + "\n"
    leg_label += r"Fit $\sigma$ = {:4.3E} $\pm$ {:4.3E}".format(stdev, stdev_err)

    return leg_label
                                
def hist_y_label(bin_width, unit):
    """
    Return a y-label for a histogram: "Events / [bin_width units]".

    Parameters
    ----------
    bin_width : float
        Bin width of histogram. 
    unit : str
        The unit of whatever distribution you are making.
        
    Returns
    -------
    y_label : str
        The y-label for the distribution.
    """
    # Automatically label y-axis.
    if bin_width < 0.01:
        y_label = "Events / [{:.2E}]".format(bin_width)
    else: 
        y_label = "Events / [{:.2f}]".format(bin_width)

    if len(unit) > 0:
        y_label = y_label.replace("]", " {}]".format(unit))

    return y_label

def make_2by2_subplots_for_ratioplots(fig_width=28, fig_height=16):
    """
    Return a figure, 4 main axes, and 4 ratio axes, useful for making ratio plots. 
    The returned tuples can be intuitively sliced to get the axes you want:
        E.g. ax_ratio_tup[0][1] -> ax12_ratio
        
    Parameters
    ----------
    fig_width : float
        Width of figure.
    fig_height : float
        Height of figure.
        
    Returns
    -------
    fig : figure
        A figure object with 8 axes objects, but 4 main areas. 
        Kind of like subplot(22), but with more control.
    ax_tup : tuple
        A tuple of 4 axes objects. These are the "main" axes for plotting.
        Each main axes object (e.g. ax21) has an associated ratio axes (e.g. ax21_ratio) associated with it. 
    ax_ratio_tup : tuple
        A tuple of 4 axes objects. These are the ratio plot axes. 
        The ratio axes make it useful to make ratio plots of the corresponding main axes object. 
    """
    fig = plt.figure(figsize=(fig_width, fig_height))

    x_left = 0.04
    y_bottom = 0.07

    leftright_spacing = 0.04  # Spacing between plots.
    updown_spacing = 0.06  # Spacing between plots.
    
    buffer_on_right = 0.03
    buffer_on_top = 0.03
    
    magnification = 2.5  # How many times taller the main plot is relative to ratio plot.
    plot_width = (1 - x_left - leftright_spacing - buffer_on_right) * 0.5
    plot_ratio_height = (1 - y_bottom - updown_spacing - buffer_on_top) / (2 * (magnification + 1))
    plot_height = magnification * plot_ratio_height

    # You can derive these formulae pretty easily. Just make a diagram.
    x_right = x_left + plot_width + leftright_spacing
    y_middle = y_bottom + plot_ratio_height + plot_height + updown_spacing

    ax11       = fig.add_axes([x_left, y_middle+plot_ratio_height, plot_width, plot_height])  # [low_left_corner_x, low_left_corner_y, plot_width, height]
    ax11_ratio = fig.add_axes([x_left, y_middle,                   plot_width, plot_ratio_height])

    ax12       = fig.add_axes([x_right, y_middle+plot_ratio_height, plot_width, plot_height])  # [low_left_corner_x, low_left_corner_y, plot_width, height]
    ax12_ratio = fig.add_axes([x_right, y_middle,                   plot_width, plot_ratio_height])

    ax21       = fig.add_axes([x_left, y_bottom+plot_ratio_height, plot_width, plot_height])  # [low_left_corner_x, low_left_corner_y, plot_width, height]
    ax21_ratio = fig.add_axes([x_left, y_bottom,                   plot_width, plot_ratio_height])

    ax22       = fig.add_axes([x_right, y_bottom+plot_ratio_height, plot_width, plot_height])  # [low_left_corner_x, low_left_corner_y, plot_width, height]
    ax22_ratio = fig.add_axes([x_right, y_bottom,                   plot_width, plot_ratio_height])
    
    ax_tup = ((ax11, ax12), (ax21, ax22))
    ax_ratio_tup = ((ax11_ratio, ax12_ratio), (ax21_ratio, ax22_ratio))
    
    return fig, ax_tup, ax_ratio_tup