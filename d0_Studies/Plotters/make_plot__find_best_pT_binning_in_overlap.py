# PURPOSE: Makes a PDF of 2 plots so that User can visually determine best pT binning:
#              plot 1 : Shows the pT bin locations vs. the number of divisions made on pT distribution. 
#              plot 2 : Shows the mean(abs(% diff.)) vs. number of divisions.
#          A KinematicBin object can make a pT hist of all events within that KinematicBin.
#          The hist is then divided up into regions on the condition that there are ~equal entries per region.
#          Where the divisions occur are the best pT bins to use for this eta region.
#          Splits pT distributions from KinematicBin objects into equal-entries divisions, to find best pT bins. 
# NOTES:   The workflow is to: 
#              (1) Make eta cut (barrel, overlap, endcap), 
#              (2) then figure out best pT bins to use per eta region
#              (3) Finally to figure out what q*d0 bins to use per (eta,pT) region.
#          % diff. = (num_entries_found_in_division - ideal_number_in_division) / ideal_number_in_division
#          For this script, User should check the eta region to run over. 
# SYNTAX:  python <script>.py
# AUTHOR:  Jake Rosenzweig
# DATE:    2020-05-04

import sys
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

# Local imports.
from d0_Utils.d0_cls import KinematicBin
from d0_Utils.d0_fns import find_equal_hist_divisions
from PyUtils.Utils_Files import makeDirs, check_overwrite


dataframe = '/Users/Jake/Desktop/MC_2016.h5'

outpath = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/find_best_binning__eta_pT_qd0/overlap_3Mevents_5to200GeV_RMS_zoom.pdf"

n_events_scan = 3000000

n_divisions_start = 1
n_divisions_end = 16

d0q_cut_ls = [-1, 1]
massZ_cut_ls=[60,120]
pT_cut_ls=[5,1000]
eta_cut_ls=[0.8, 1.0]

d0_type='BS'
dR_cut=0.008
use_ptotal_instead = False
verbose = True
overwrite = False

bin_limits = [5,200,1]
x_limits=[-50,250]
y_max = 70

#--------------------------------#
#---------- Automatons ----------#
#--------------------------------#
dir_ = os.path.dirname(outpath)
check_overwrite(outpath,overwrite)
makeDirs(dir_) 

plt.style.use("cmsstyle_plot")

df = pd.read_hdf(dataframe)

kb = KinematicBin(df_original=df,
                  n_evts=n_events_scan,
                  massZ_cut_ls=massZ_cut_ls,
                  eta_cut_ls=eta_cut_ls, 
                  pT_cut_ls=pT_cut_ls, 
                  d0q_cut_ls=d0q_cut_ls,
                  d0_type=d0_type,
                  dR_cut=dR_cut,
                  use_ptotal_instead=use_ptotal_instead, 
                  verbose=verbose)

with PdfPages(outpath) as pdf:

    # Must make distribution to divide it up into regions.
    kb.plot_1D_kinematics(kinem="pT1", lep_selection_type="independent", 
                          x_limits=x_limits, bin_limits=bin_limits, 
                          run_over_only_n_evts=-1, 
                          ax=None, x_label="", y_label="", title="", y_max=-1, log_scale=False, 
                          iter_gaus=(False, 0))
    pdf.savefig()
    plt.close("all")
    
    bin_vals = kb.stats_dict['pT1']['bin_vals']
    bin_edges = kb.stats_dict['pT1']['bin_edges']

    f, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=False, figsize=(9.6,9.6))
    eta_cuts = kb.cut_dict['eta']

    division_ls = list(range(n_divisions_start, n_divisions_end+1))
    mean_perc_diff_ls = []

    for divisions in division_ls:
        pT_bins, this_dict = find_equal_hist_divisions(bin_edges, bin_vals, divisions, verbose=verbose)

        abs_perc_diff_ls = [abs(perc) for perc in this_dict.values()] 
#        mean_perc_diff = np.mean(abs_perc_diff_ls)
        mean_perc_diff = np.sqrt(np.mean(np.power(abs_perc_diff_ls, 2)))
        mean_perc_diff_ls.append(mean_perc_diff)

        x_vals = np.ones(len(pT_bins)) * divisions
        abs_perc_diff_ls.insert(0,0)

        # mfc = markerface color
        ax1.errorbar(x=x_vals, y=pT_bins, yerr=abs_perc_diff_ls, linestyle="", markersize=2, mfc='k')
        if y_max != -1:
            ax1.set_ylim([0,y_max])

    ax2.plot(division_ls, mean_perc_diff_ls)

    ax1.set_title(r"%s" % eta_cuts)
    ax1.set_ylabel(r"$p_T$ bins [GeV]")
    ax2.set_xlabel(r"Equal-entry divisions in $p_T$ distribution")
    ax2.set_ylabel(r"RMS(% diff.)")

    pdf.savefig()
    plt.close("all")

#with 
