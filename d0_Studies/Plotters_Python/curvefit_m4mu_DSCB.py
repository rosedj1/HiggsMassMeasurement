"""
PURPOSE:
    This script unpacks the m4mu binned data from a TH1 in a root file.
    It then takes those data and bins them in a Python hist, using the same 
    bin values, bin edges, etc. 
    Next, a double-sided Crystal Ball fit is performed.
    Finally the histogram and fit are saved into a PDF.
SYNTAX:
    python <script>.py
AUTHOR:
    Jake Rosenzweig
UPDATED:
    2020-06-17
"""
import os
import ROOT

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from matplotlib import gridspec

# Local imports.
from Utils_Python.Utils_Plotting import make_1D_dist
from Utils_Python.Utils_StatsAndFits import (crystal_ball_doublesided_func, fit_with_crystal_ball_doublesided,
                                             exp_gaus_exp_func, fit_with_exp_gaus_exp, prop_err_x_div_y)
from Utils_Python.Utils_Files import make_dirs, check_overwrite
from d0_Utils.d0_fns import centers_of_binning_array, make_binning_array, calc_x_err_bins_from_bin_edges

# Make sure sexy Times New Roman font loads correctly.
import matplotlib.font_manager as font_manager; font_manager._rebuild()
#-----------------------#
#----- User Params -----#
#-----------------------#
infile_path = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/m4mu_withcorr_MC_2017_ggF_fullstats.root"
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/DSCB_fit_m4mu/"
filename = "MC_2017_m4mu_DSCBfit_beforeandafter_pTcorr_test03"
overwrite = False
x_limits = [105, 140]
x_bin_info = [105, 140, 0.25]  # min, max, bin_width
fit_with_DSCB = True
fit_with_ExpGausExp = False

#----------------------------#
#----- Script Functions -----#
#----------------------------#
def prep_area(outdir_plots, filename, overwrite):
    fullpath_pdf = os.path.join(outdir_plots, filename + ".pdf")
    check_overwrite(fullpath_pdf, overwrite)  # See if files already exist.
    make_dirs(outdir_plots)  # Make directories, if need be.
    return fullpath_pdf

def get_data_from_TTree(infile_path):
    print(f"Opening root file:\n{infile_path}")
    infile = ROOT.TFile.Open(infile_path, "read")
    tree = infile.Get("tree")

    print(f"...Unpacking m4mu vals from TTree...")
    data = []
    for evt in tree:
        data.append(evt.m4mu)
    return data

def get_corranduncorr_data_from_TTree(infile_path):
    print(f"Opening root file:\n{infile_path}")
    infile = ROOT.TFile.Open(infile_path, "read")
    tree = infile.Get("tree")

    print(f"...Unpacking m4mu vals from TTree...")
    data = []
    data_corr = []
    for evt in tree:
        data.append(evt.m4mu)
        data_corr.append(evt.m4mu_corr)
    return data, data_corr

def make_fig_ax():
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1], hspace=0.07)
    ax_upper = fig.add_subplot(gs[0])
    ax_upper.xaxis.set_major_formatter(plt.NullFormatter())
    
    ax_lower = fig.add_subplot(gs[1])
    return fig, ax_upper, ax_lower

def draw_m4mu_dist(data, x_limits, x_bin_info, ax, color, extra_leg_text=None, leg_loc=None):
    x_bin_edges, binw = make_binning_array(x_bin_info)

    x_label = ""
    y_label = r"Events / [%.2f GeV]" % binw
    title = r"Effect of $p_T$ corrections on Higgs Width (Production Mode: ggF)"
    print(f"Making Python histogram.")
    
    ax, bin_vals, bin_edges, stats = make_1D_dist(
                                        ax, data, 
                                        x_limits = x_limits, 
                                        x_bin_edges = x_bin_edges,
                                        x_label = x_label, 
                                        y_label = y_label, 
                                        title = title, 
                                        extra_leg_text=extra_leg_text,
                                        y_max=-1, log_scale=False,
                                        color=color, leg_loc=leg_loc)
    return ax, bin_vals, bin_edges, stats

def get_DSCB_params(bin_centers, bin_vals):
    print(f"Performing DSCB fit.")
    popt, popt_err, pcov = fit_with_crystal_ball_doublesided(
                                        bin_centers, bin_vals,
                                        guess_params=[2500,1.3,2,1.3,2,125,1]
                                        )
    return popt, popt_err, pcov
    
def do_ExpGausExp_fit(bin_centers, bin_vals):
    print(f"Performing Exp-Gaus-Exp fit.")
    popt, popt_err, pcov = fit_with_exp_gaus_exp(
                                        bin_centers, bin_vals,
                                        guess_params=[2500,1,1,125,1]
                                        )
    return popt, popt_err, pcov

def make_leg_text_DSCB_fit_params(popt, popt_err, extra_text=None):
    """Returns the legend text which contains the best-fit params of a 
    double-sided Crystal Ball fit."""
    # Remember that popt[0] is the best-fit coefficient (scaling factor).
    leg_fit_text  = r"DSCB fit" 
    if extra_text is not None:
        leg_fit_text += extra_text
    leg_fit_text += ":"
    leg_fit_text += "\n" + r"$\alpha_{L}$ = %.4f $\pm$ %.4f" % (popt[1], popt_err[1])
    leg_fit_text += "\n" + r"$\alpha_{R}$ = %.4f $\pm$ %.4f" % (popt[3], popt_err[3])
    leg_fit_text += "\n" + r"$n_{L}$ = %.4f $\pm$ %.4f" % (popt[2], popt_err[2])
    leg_fit_text += "\n" + r"$n_{R}$ = %.4f $\pm$ %.4f" % (popt[4], popt_err[4])
    leg_fit_text += "\n" + r"$\mu$ = %.4f $\pm$ %.4f GeV" % (popt[5], popt_err[5])
    leg_fit_text += "\n" + r"$\sigma$ = %.4f $\pm$ %.4f GeV" % (popt[6], popt_err[6])

    return leg_fit_text

def make_leg_text_EGE_fit_params(popt, popt_err):
    """Returns the legend text which contains the best-fit params of a 
    ExpGausExp fit."""
    leg_fit_text  = r"ExpGausExp fit:"
    leg_fit_text += "\n" + r"$k_{L}$ = %.4f $\pm$ %.4f" % (popt[1], popt_err[1])
    leg_fit_text += "\n" + r"$k_{R}$ = %.4f $\pm$ %.4f" % (popt[2], popt_err[2])
    leg_fit_text += "\n" + r"$\mu$ = %.4f $\pm$ %.4f GeV" % (popt[3], popt_err[3])
    leg_fit_text += "\n" + r"$\sigma$ = %.4f $\pm$ %.4f GeV" % (popt[4], popt_err[4])

    return leg_fit_text

def plot_DSCB_fit(ax, bin_edges, bin_vals, color, x_limits=None, linewidth=1.25, extra_leg_text=None):
    bin_centers = centers_of_binning_array(bin_edges)
    fit_x_vals = np.linspace(bin_edges[0], bin_edges[-1], 1000)

    popt, popt_err, pcov = get_DSCB_params(bin_centers, bin_vals)
    leg_fit_text = make_leg_text_DSCB_fit_params(popt, popt_err, extra_text=extra_leg_text)

    y_vals_fit = crystal_ball_doublesided_func(fit_x_vals, *popt)

    ax.plot(fit_x_vals, y_vals_fit, c=color, label=leg_fit_text, linewidth=linewidth, marker="None")
    if x_limits is not None:
        ax.set_xlim(x_limits)
    ax.legend()
    return ax

def plot_EGE_fit(ax, bin_edges, bin_vals, color):
    bin_centers = centers_of_binning_array(bin_edges)
    fit_x_vals = np.linspace(bin_edges[0], bin_edges[-1], 1000)

    popt, popt_err, pcov = do_ExpGausExp_fit(bin_centers, bin_vals)
    leg_fit_text = make_leg_text_EGE_fit_params(popt, popt_err)

    y_vals_fit = exp_gaus_exp_func(fit_x_vals, *popt)

    ax.plot(fit_x_vals, y_vals_fit, c=color, label=leg_fit_text, linewidth=1, marker="None")
    ax.legend()
    return ax

def do_DSCB_and_EGE_fits(ax, bin_edges, bin_vals):
    # bin_centers = centers_of_binning_array(bin_edges)
    # fit_x_vals = np.linspace(bin_edges[0], bin_edges[-1], 1000)
    plot_DSCB_fit(ax, bin_edges, bin_vals, color="r")
        # popt_dscb, popt_err_dscb, pcov_dscb = get_DSCB_params(bin_centers, bin_vals)
        # leg_fit_text_dscb = make_leg_text_DSCB_fit_params(popt_dscb, popt_err_dscb)

        # y_vals_fit_dscb = crystal_ball_doublesided_func(fit_x_vals, *popt_dscb)

        # ax.plot(fit_x_vals, y_vals_fit_dscb, c='r', label=leg_fit_text_dscb, linewidth=1.25, marker="None")
        # ax.legend()
    plot_EGE_fit(ax, bin_edges, bin_vals, color="g")
    # popt_ege, popt_err_ege, pcov_ege = do_ExpGausExp_fit(bin_centers, bin_vals)
    # leg_fit_text_ege = make_leg_text_EGE_fit_params(popt_ege, popt_err_ege)

    # y_vals_fit_ege = exp_gaus_exp_func(fit_x_vals, *popt_ege)

    # ax.plot(fit_x_vals, y_vals_fit_ege, c='g', label=leg_fit_text_ege, linewidth=1, marker="None")
    # ax.legend()
    return ax

def decorate_ratio_plot(ax, x_limits):
    ax.set_xlabel(r"$m_{4\mu}$ [GeV]")
    ax.set_ylabel(r"$\frac{\mathrm{distrib.^{after corr.}}}{\mathrm{distrib.}}$")
    ax.axhline(y=1, xmin=-0.1, xmax=1.1, color='r', lw=0.7)
    ax.xaxis.set_major_formatter(plt.NullFormatter())
    ax.set_xlim(x_limits)
    ax.set_ylim([-0.2,2.2])
    return ax

def plot_ratio(ax, bin_edges, numer, denom, x_limits, color=None):
    ratio, ratio_err = prop_err_x_div_y(numer, denom, np.sqrt(numer), np.sqrt(denom))

    x_err_low, x_err_high = calc_x_err_bins_from_bin_edges(bin_edges)
    x_vals = centers_of_binning_array(bin_edges)

    color = "blue" if color is None else color
    markerstyle = "none"
    capsize = 1
    errlw = markererrw = 0.50
    
    ax.errorbar(x_vals, ratio, xerr=[x_err_low, x_err_high], yerr=ratio_err, 
                        capsize=capsize, fmt=markerstyle, 
                        mfc=color, mec=color, ecolor=color, 
                        elinewidth=errlw, mew=markererrw)

    ax = decorate_ratio_plot(ax, x_limits)
    return ax

if __name__ == "__main__":
    # Pre-analysis.
    fullpath_pdf = prep_area(outdir_plots, filename, overwrite)
    plt.style.use('/ufrc/avery/rosedj1/HiggsMassMeasurement/Utils_Python/Plot_Styles_Matplotlib/cmsstyle_plot.mplstyle')
    
    # data = get_data_from_TTree(infile_path)
    data, data_corr = get_corranduncorr_data_from_TTree(infile_path)

    fig, ax_upper, ax_lower = make_fig_ax()
    ax_upper, bin_vals,      bin_edges, stats      = draw_m4mu_dist(
                                                            data,      x_limits, x_bin_info, ax_upper, 
                                                            color="black", leg_loc="upper right")
    ax_upper, bin_vals_corr, bin_edges, stats_corr = draw_m4mu_dist(data_corr, x_limits, x_bin_info, ax_upper, 
                                                            color="blue", extra_leg_text=r" after corr.",leg_loc="lower right")

    ax_upper = plot_DSCB_fit(ax_upper, bin_edges, bin_vals, x_limits=x_limits, color="r", linewidth=1.10)
    ax_upper = plot_DSCB_fit(ax_upper, bin_edges, bin_vals_corr, x_limits=x_limits,
                             color="g", linewidth=0.95, extra_leg_text=r" after corr.")

    ax_lower = plot_ratio(ax_lower, bin_edges, numer=bin_vals_corr, denom=bin_vals, x_limits=x_limits, color="black")
    # if (fit_with_DSCB) and (fit_with_ExpGausExp):
    #     ax_upper = do_DSCB_and_EGE_fits(ax_upper, bin_edges, bin_vals)

    plt.savefig(fullpath_pdf)
    print(f"PDF created: {fullpath_pdf}")