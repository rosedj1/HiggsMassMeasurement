#FIXME: NOT FIXED YET!

#---------------------------------------------------#
#----- Make 5 x 3 grid of dPhi vs. dEta plots. -----#
#---------------------------------------------------#

import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('/Users/Jake/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

# Neat tricks.
from matplotlib.backends.backend_pdf import PdfPages

# Local imports. 
#from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from PyUtils.Utils_Files import makeDirs 
#from PyUtils.Utils_Plotting import make_2by2_subplots_for_ratioplots
from d0_Utils.d0_cls import KinematicBin

# FIXME: cannot make these 2D lists of kbins easily. 
from d0_Studies.kinematic_bins import (kbin_ls_eta_INCLUSIVE_pT_LOW1,  kbin_ls_eta_INCLUSIVE_pT_LOW2, 
                                       kbin_ls_eta_INCLUSIVE_pT_MED1,  kbin_ls_eta_INCLUSIVE_pT_MED2, 
                                       kbin_ls_eta_INCLUSIVE_pT_HIGH1, kbin_ls_eta_INCLUSIVE_pT_HIGH2)

plt.rcParams.update({'figure.max_open_warning': 10})    # Related to saving memory and opening plots.
plt.style.use('cmsstyle_plot')

outfile = "/Users/Jake/Desktop/20200426/grid_of_2D_plots/dPhivsdEtaANDdTheta_TEST1.pdf"

n_evts_scan = 10000
n_evts_keep = 1000

infile_path_MC_2016 = '/Users/Jake/Desktop/MC_2016.h5'

df_MC_2016 = pd.read_hdf(infile_path_MC_2016)

def make_grid_of_2D_plots(x_kinem, y_kinem, kbin_2D_ls, x_bin_limits, y_bin_limits, pdf):
    """
    Create a grid of 2D plots. User can choose the x and y variables to be plotted. 
    The shape of kbin_2D_ls determines the shape of the plots on the figure. 
    
    Parameters
    ----------
    x_kinem : str
        The x kinematic variable.
    y_kinem : str
        The y kinematic variable.
    kbin_2D_ls : 2D array of KinemBinnedEtaPt objects
        Shape of this array determines shape of plots on the figure. 
        Example: [[kbin1, kbin2], 
                  [kbin3, kbin4]]
        would make a figure with 2x2 plots on it. 
    x_bin_limits : list
        x-axis bin limits [first_bin_left_edge, last_bin_right_edge, bin_width]
    y_bin_limits : list
        y-axis bin limits [first_bin_left_edge, last_bin_right_edge, bin_width]
    pdf : pdf save object
        PdfPages object. 
    """
    plt.style.use("grid_multiple_plots")
    
    rows = np.shape(kbin_2D_ls)[0]
    cols = np.shape(kbin_2D_ls)[1]

    f, ax = plt.subplots(rows, cols)#, figsize=(19.2, 14.4))

    # Flatten arrays down for easy zipping.
    kbin_flat = [kb for sublist in kbin_2D_ls for kb in sublist]
    ax_flat = ax.flatten()

    for kbin, ax in zip(kbin_flat, ax_flat):
        kbin.make_2D_plot(x_kinem, 
                          y_kinem, 
                          x_bin_limits=x_bin_limits,
                          y_bin_limits=x_bin_limits,
                          lep_selection_type="independent",
                          run_over_only_n_evts=-1,
                          title="",
                          ax=ax,
                          verbose=False)

    plt.tight_layout(w_pad=35)
    pdf.savefig()
    plt.close("all")

mydir = os.path.dirname(outfile)
makeDirs(mydir)

bin_limits = [-0.01, 0.01, 0.0001]

with PdfPages(outfile) as pdf:
    make_grid_of_2D_plots("delta_eta",   "delta_phi", kbin_ls_eta_INCLUSIVE_pT_LOW1,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_theta", "delta_phi", kbin_ls_eta_INCLUSIVE_pT_LOW1,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_eta",   "delta_phi", kbin_ls_eta_INCLUSIVE_pT_LOW2,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_theta", "delta_phi", kbin_ls_eta_INCLUSIVE_pT_LOW2,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_eta",   "delta_phi", kbin_ls_eta_INCLUSIVE_pT_MED1,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_theta", "delta_phi", kbin_ls_eta_INCLUSIVE_pT_MED1,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_eta",   "delta_phi", kbin_ls_eta_INCLUSIVE_pT_MED2,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_theta", "delta_phi", kbin_ls_eta_INCLUSIVE_pT_MED2,  bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_eta",   "delta_phi", kbin_ls_eta_INCLUSIVE_pT_HIGH1, bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_theta", "delta_phi", kbin_ls_eta_INCLUSIVE_pT_HIGH1, bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_eta",   "delta_phi", kbin_ls_eta_INCLUSIVE_pT_HIGH2, bin_limits, bin_limits, pdf)
    make_grid_of_2D_plots("delta_theta", "delta_phi", kbin_ls_eta_INCLUSIVE_pT_HIGH2, bin_limits, bin_limits, pdf)