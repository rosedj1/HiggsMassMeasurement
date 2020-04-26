import sys
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.path.append('/Users/Jake/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

# Neat tricks.
#from itertools import chain
from matplotlib.backends.backend_pdf import PdfPages
#from matplotlib.patches import Rectangle
#from IPython.display import display
#from IPython.core.display import display, HTML
#display(HTML("<style>.container { width:100% !important; }</style>"))

# Local imports. 
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
#from PyUtils.Utils_Plotting import (change_cmap_bkg_to_white, save_plots_to_outpath, make_1D_dist, get_stats_1Dhist, 
#                                    get_stats_2Dhist, hist_y_label, make_2by2_subplots_for_ratioplots,
#                                    add_underoverflow_entries, make_stats_legend_for_1dhist, make_stats_legend_for_2dhist, 
#                                    make_stats_legend_for_gaus_fit)
from PyUtils.Utils_Plotting import make_2by2_subplots_for_ratioplots
#from PyUtils.Utils_Physics import theta2pseudorap, pseudorap2theta, calc_dR, calc_dphi
#from PyUtils.Utils_StatsAndFits import gaussian, fit_with_gaussian, iterative_fit_gaus
#from PyUtils.Utils_Collection_Helpers import weave_lists
from d0_Utils.d0_cls import KinematicBin
#from d0_Utils.d0_fns import (make_binning_array, shift_binning_array, get_subset_mask, make_kinem_subplot)
#from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict

plt.rcParams.update({'figure.max_open_warning': 10})    # Related to saving memory and opening plots.

plt.style.use('cmsstyle_plot')
# plt.style.use("grid_multiple_plots")



### Make all kinematic plots into 1 PDF with iterative fits.

# def make_all_kinem_plots(kbin, iterations):
#     with PdfPages("/Users/Jake/Desktop/20200415/kinematics_and_fits/practice1.pdf") as pdf:

outpath = "/Users/Jake/Desktop/20200426/test8.pdf"

make_this_dir = os.path.dirname(outpath)
makeDirs(make_this_dir)

n_evts_scan = 100000
n_evts_keep = 10000

infile_path_MC_2016 = '/Users/Jake/Desktop/MC_2016.h5'

df_MC_2016 = pd.read_hdf(infile_path_MC_2016)

kbin = KinematicBin(df_MC_2016, 
                        n_evts=n_evts_scan, 
                        massZ_cut_ls=[60,120],
                        eta_cut_ls=[0.0, 0.3], 
                        pT_cut_ls=[5, 100], 
                        d0q_cut_ls=[-0.002, 0.002],
                        d0_type='BS',
                        dR_cut=0.02,
                        use_ptotal_instead=False, 
                        verbose=True)


with PdfPages(outpath) as pdf:

    # Rec vs. Gen plots require special fig and ax. 
    fig, ax_tup, ax_ratio_tup = make_2by2_subplots_for_ratioplots()
    kbin.plot_kinem_genrec_comparison("genLep_eta1","eta1", lep_selection_type='1', x_limits=[-2.5,2.5], bin_limits=[-2.4, 2.4, 0.05], run_over_only_n_evts=-1, ax=ax_tup[0][0], ax_ratio=ax_ratio_tup[0][0], log_scale=False)    
    kbin.plot_kinem_genrec_comparison("genLep_eta2","eta2", lep_selection_type='2', x_limits=[-2.5,2.5], bin_limits=[-2.4, 2.4, 0.05], run_over_only_n_evts=-1, ax=ax_tup[0][1], ax_ratio=ax_ratio_tup[0][1], log_scale=False) 
    kbin.plot_kinem_genrec_comparison("genLep_phi1","phi1", lep_selection_type='1', x_limits=[-np.pi-0.2, np.pi+0.2], bin_limits=[-np.pi, np.pi, 0.05], run_over_only_n_evts=-1, ax=ax_tup[1][0], ax_ratio=ax_ratio_tup[1][0], log_scale=False)
    kbin.plot_kinem_genrec_comparison("genLep_phi2","phi2", lep_selection_type='2', x_limits=[-np.pi-0.2, np.pi+0.2], bin_limits=[-np.pi, np.pi, 0.05], run_over_only_n_evts=-1, ax=ax_tup[1][1], ax_ratio=ax_ratio_tup[1][1], log_scale=False)
    pdf.savefig()
    plt.close("all")
    
    fig, ax_tup, ax_ratio_tup = make_2by2_subplots_for_ratioplots()
    kbin.plot_kinem_genrec_comparison("genLep_pt1","pT1", lep_selection_type='1', x_limits=[0,150], bin_limits=[0, 120, 1], run_over_only_n_evts=-1, ax=ax_tup[0][0], ax_ratio=ax_ratio_tup[0][0], log_scale=False)
    kbin.plot_kinem_genrec_comparison("genLep_pt2","pT2", lep_selection_type='2', x_limits=[0,150], bin_limits=[0, 120, 1], run_over_only_n_evts=-1, ax=ax_tup[0][1], ax_ratio=ax_ratio_tup[0][1], log_scale=False)
    kbin.plot_kinem_genrec_comparison("genLep_eta1","eta1", lep_selection_type='1', x_limits=[-2.5,2.5], bin_limits=[-2.4, 2.4, 0.05], run_over_only_n_evts=-1, ax=ax_tup[1][0], ax_ratio=ax_ratio_tup[1][0], log_scale=False)    
    kbin.plot_kinem_genrec_comparison("genLep_eta2","eta2", lep_selection_type='2', x_limits=[-2.5,2.5], bin_limits=[-2.4, 2.4, 0.05], run_over_only_n_evts=-1, ax=ax_tup[1][1], ax_ratio=ax_ratio_tup[1][1], log_scale=False)    
    pdf.savefig()
    plt.close("all")
    
    # Regular kinematic distributions.
    fig = plt.figure(figsize=(28,16))
    ax = plt.subplot(221)
    kbin.plot_1D_kinematics(kinem="delta_eta1", lep_selection_type='independent', x_limits=[-0.015, 0.015], bin_limits=[-0.012, 0.012, 0.0001], run_over_only_n_evts=-1, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(222)
    kbin.plot_1D_kinematics(kinem="delta_eta2", lep_selection_type='independent', x_limits=[-0.015, 0.015], bin_limits=[-0.012, 0.012, 0.0001], run_over_only_n_evts=-1, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(223)
    kbin.plot_1D_kinematics(kinem="delta_theta1", lep_selection_type='independent', x_limits=[-0.002, 0.002], bin_limits=[-0.0018, 0.0018, 0.00002], run_over_only_n_evts=-1, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(224)
    kbin.plot_1D_kinematics(kinem="delta_theta2", lep_selection_type='independent', x_limits=[-0.002, 0.002], bin_limits=[-0.0018, 0.0018, 0.00002], run_over_only_n_evts=-1, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    pdf.savefig()
    plt.close("all")
    
    fig = plt.figure(figsize=(28,16))
    ax = plt.subplot(221)
    kbin.plot_1D_kinematics(kinem="delta_phi1", lep_selection_type='independent', x_limits=[-0.003, 0.003], bin_limits=[-0.003, 0.003, 0.00005], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))    
    ax = plt.subplot(222)
    kbin.plot_1D_kinematics(kinem="delta_phi2", lep_selection_type='independent', x_limits=[-0.003, 0.003], bin_limits=[-0.003, 0.003, 0.00005], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(223)
    kbin.plot_1D_kinematics(kinem="delta_R1", lep_selection_type='independent', x_limits=[-0.004, 0.008], bin_limits=[0, 0.008, 0.00005], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(False, 5))
    ax = plt.subplot(224)
    kbin.plot_1D_kinematics(kinem="delta_R2", lep_selection_type='independent', x_limits=[-0.004, 0.008], bin_limits=[0, 0.008, 0.00005], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(False, 5))
    pdf.savefig()
    plt.close("all")

    #--- Momentum plots ---#
    fig = plt.figure(figsize=(28,16))
    ax = plt.subplot(221)
    kbin.plot_1D_kinematics(kinem="delta_pT1", lep_selection_type='both', x_limits=[-20, 20], bin_limits=[-20, 20, 0.2], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(222)
    kbin.plot_1D_kinematics(kinem="delta_pT2", lep_selection_type='both', x_limits=[-20, 20], bin_limits=[-20, 20, 0.2], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(223)
    kbin.plot_1D_kinematics(kinem="massZ", lep_selection_type='both', x_limits=[0,0], bin_limits=[0,0,0], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(224)
    kbin.plot_1D_kinematics(kinem="massZErr", lep_selection_type='both', x_limits=[0,0], bin_limits=[0,0,0], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(False, 4))
    pdf.savefig()
    plt.close("all")
    
    fig = plt.figure(figsize=(28,16))
    ax = plt.subplot(221)
    kbin.plot_1D_kinematics(kinem="delta_pToverpT1", lep_selection_type='1', x_limits=[-0.3, 0.3], bin_limits=[-0.3, 0.3, 0.004], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 3))
    ax = plt.subplot(222)
    kbin.plot_1D_kinematics(kinem="delta_pToverpT2", lep_selection_type='2', x_limits=[-0.3, 0.3], bin_limits=[-0.3, 0.3, 0.004], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(223)
    kbin.plot_1D_kinematics(kinem="delta_pToverRecpT1", lep_selection_type='both', x_limits=[-0.3, 0.3], bin_limits=[-0.3, 0.3, 0.004], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(224)
    kbin.plot_1D_kinematics(kinem="delta_pToverRecpT2", lep_selection_type='either', x_limits=[-0.3, 0.3], bin_limits=[-0.3, 0.3, 0.004], run_over_only_n_evts=-1, ax=ax, log_scale=False, iter_gaus=(True, 4))
    pdf.savefig()
    plt.close("all")

    #--- d0 Plots ---#
    fig = plt.figure(figsize=(28,16))
    ax = plt.subplot(221)
    kbin.plot_1D_kinematics(kinem="d0BSq1", lep_selection_type='1', x_limits=[-0.003, 0.003], bin_limits=[-0.003, 0.003, 0.00005], run_over_only_n_evts=n_evts_keep, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(222)
    kbin.plot_1D_kinematics(kinem="d0BSq2", lep_selection_type='2', x_limits=[-0.003, 0.003], bin_limits=[-0.003, 0.003, 0.00005], run_over_only_n_evts=n_evts_keep, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(223)
    kbin.plot_1D_kinematics(kinem="d0BSq1", lep_selection_type='both', x_limits=[-0.003, 0.003], bin_limits=[-0.003, 0.003, 0.00005], run_over_only_n_evts=n_evts_keep, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    ax = plt.subplot(224)
    kbin.plot_1D_kinematics(kinem="d0BSq2", lep_selection_type='either', x_limits=[-0.003, 0.003], bin_limits=[-0.003, 0.003, 0.00005], run_over_only_n_evts=n_evts_keep, ax=ax, y_max=-1, log_scale=False, iter_gaus=(True, 4))
    pdf.savefig()
    plt.close("all")
    
    fig = plt.figure(figsize=(28,16))
    ax = plt.subplot(221)
    kbin.make_2D_plot("delta_eta","delta_phi",lep_selection_type='1', x_bin_limits=[-0.003, 0.003, 0.00003],y_bin_limits=[-0.003, 0.003, 0.00003],run_over_only_n_evts=-1, title="", exclusive=True,save_plot=False,save_as_png=False,verbose=False,outpath=outpath, ax=ax)
    ax = plt.subplot(222)
    kbin.make_2D_plot("delta_eta","delta_phi",lep_selection_type='2', x_bin_limits=[-0.003, 0.003, 0.00003],y_bin_limits=[-0.003, 0.003, 0.00003],run_over_only_n_evts=-1, title="", exclusive=True,save_plot=False,save_as_png=False,verbose=False,outpath=outpath, ax=ax)
    ax = plt.subplot(223)
    kbin.make_2D_plot("delta_eta","delta_phi",lep_selection_type='independent', x_bin_limits=[-0.003, 0.003, 0.00003],y_bin_limits=[-0.003, 0.003, 0.00003],run_over_only_n_evts=-1, title="", exclusive=True,save_plot=False,save_as_png=False,verbose=False,outpath=outpath, ax=ax)
    ax = plt.subplot(224)
    kbin.make_2D_plot("delta_eta","delta_phi",lep_selection_type='both', x_bin_limits=[-0.003, 0.003, 0.00003],y_bin_limits=[-0.003, 0.003, 0.00003],run_over_only_n_evts=-1, title="", exclusive=True,save_plot=False,save_as_png=False,verbose=False,outpath=outpath, ax=ax)
    plt.tight_layout()
    pdf.savefig()
    plt.close("all")