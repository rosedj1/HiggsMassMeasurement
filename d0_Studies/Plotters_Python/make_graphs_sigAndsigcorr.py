"""
# This script 
"""
import pickle
import os
from matplotlib import gridspec
from matplotlib.backends.backend_pdf import PdfPages

import numpy as np
import matplotlib.pyplot as plt

# Local imports.
from Utils_Python.Utils_StatsAndFits import Andrey_prop_err_on_dsigoversig
from Utils_Python.Utils_Files import check_overwrite, make_str_title_friendly
from Utils_Python.printing import print_header_message
from d0_Utils.d0_fns import calc_x_err_bins_from_bin_edges
#----- Main -----#
inpath_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/20200709_MC2017_combinesamples_applycorr_fullstats_2p00sigmas__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"
with open(inpath_pkl, "rb") as f:
    dpToverpT_comb_stats_dict = pickle.load(f)

outdir_pdf = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/sigma_itergaus_pTcorr/"
filename = "20200709_finalplots_sigmacorrafterpTcorr"
year = "2017"
overwrite = False
set_log_x = True
y_limits_mainplot = [0.005, 0.053]

eta_ls = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.25, 1.5, 1.75, 2.0, 2.1, 2.2, 2.3, 2.4]
pT_ls  = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]

color_sigma = "blue"
color_sigma_corr = "red"
markerstyle = "none"
capsize = 1
errlw = markererrw = 0.50

# Include all 
sample_str = r"$J/\psi$, DY"
#----- Automatons -----#
plt.style.use("/ufrc/avery/rosedj1/HiggsMassMeasurement/Utils_Python/Plot_Styles_Matplotlib/cmsstyle_plot.mplstyle")

suffix   = (
    f"__{min(eta_ls)}eta{max(eta_ls)}"
    f"__{min(pT_ls):.1f}pT{max(pT_ls):.1f}_GeV"
)
suffix = make_str_title_friendly(suffix)

outpath_pdf = os.path.join(outdir_pdf, f"MC{year}_{filename}{suffix}.pdf")
check_overwrite(outpath_pdf, overwrite)

def get_stats_info(eta_bin, pT_bin, stats_dict, stat, hist_type):
    eta_min = eta_bin[0]
    eta_max = eta_bin[1]
    pT_min = pT_bin[0]
    pT_max = pT_bin[1]
    
    key = f"h_{eta_min}eta{eta_max}_{pT_min}pT{pT_max}_combined_{hist_type}"
    return stats_dict[key][stat][-1]

str_sigma = "\sigma_{\mathrm{\Delta p_{T}/p_{T}}}"
str_sigma_corr = "\sigma_{\mathrm{\Delta p_{T}/p_{T}}}^{\mathrm{corr.}}"

bins = np.array(pT_ls)
x = (bins[1:] + bins[:-1]) / 2.

# Retrieve data.
with PdfPages(outpath_pdf) as pdf:
    fig = plt.figure()
    for eta_bin in zip(eta_ls[:-1], eta_ls[1:]):
        sigma_ls = []
        sigma_corr_ls = []
        sigma_err_ls = []
        sigma_corr_err_ls = []

        print(f"Running over eta_bin={eta_bin}")
        for pT_bin in zip(pT_ls[:-1], pT_ls[1:]):
            sig = get_stats_info(eta_bin=eta_bin, pT_bin=pT_bin, 
                        stats_dict=dpToverpT_comb_stats_dict, stat='std_ls', hist_type="dpToverpT")
            sig_err = get_stats_info(eta_bin=eta_bin, pT_bin=pT_bin, 
                        stats_dict=dpToverpT_comb_stats_dict, stat='std_err_ls', hist_type="dpToverpT")
            sig_corr = get_stats_info(eta_bin=eta_bin, pT_bin=pT_bin, 
                        stats_dict=dpToverpT_comb_stats_dict, stat='std_ls', hist_type="dpToverpTcorr")
            sig_corr_err = get_stats_info(eta_bin=eta_bin, pT_bin=pT_bin, 
                        stats_dict=dpToverpT_comb_stats_dict, stat='std_err_ls', hist_type="dpToverpTcorr")
            
            sigma_ls.append(sig)
            sigma_corr_ls.append(sig_corr)
            sigma_err_ls.append(sig_err)
            sigma_corr_err_ls.append(sig_corr_err)
        # End pT_bin loop. Data collected.

        sigma_arr = np.array(sigma_ls)
        sigma_corr_arr = np.array(sigma_corr_ls)
        sigma_err_arr = np.array(sigma_err_ls)
        sigma_corr_err_arr = np.array(sigma_corr_err_ls)

        # Design the plot.
        gs = gridspec.GridSpec(2, 1, height_ratios=[2.5, 1], hspace=0.11) # Allow for axes of various sizes
        ax = fig.add_subplot(gs[0]) # Build first axis with aspect ratio gs[0]

        x_min = bins.min()
        x_max = bins.max()
        x_err_low, x_err_high = calc_x_err_bins_from_bin_edges(bins)

        if (set_log_x):
            ax.set_xscale('log')
        ax.set_xlim([x_min, x_max])
        if y_limits_mainplot is None:
            # Auto-set the y limits.
            y_lim = [0, max(sigma_arr.max(), sigma_corr_arr.max())+0.01]
        else:
            y_lim = y_limits_mainplot
        ax.set_ylim(y_lim)
        ax.errorbar(x, sigma_arr, xerr=[x_err_low, x_err_high], yerr=sigma_err_arr, label=r"$%s$" % str_sigma, 
                    color=color_sigma, mfc=color_sigma, ecolor=color_sigma,
                    capsize=capsize, fmt=markerstyle, elinewidth=errlw, mew=markererrw)

        ax.errorbar(x, sigma_corr_arr, xerr=[x_err_low, x_err_high], yerr=sigma_corr_err_arr, label=r"$%s$" % str_sigma_corr, 
                    color=color_sigma_corr, mfc=color_sigma_corr, ecolor=color_sigma_corr, 
                    capsize=capsize, fmt=markerstyle, elinewidth=errlw, mew=markererrw)
        #                    mec=color_sigma_corr,

        ax.set_ylabel(r'$%s$' % str_sigma)
        
        title_str  = r"Improvement of $%s$ from $\mathrm{p_{T,\mu}}$ corrections:" % str_sigma + "\n" 
        title_str += r"$ %.1f < \left|\eta\right| < %.1f$  " % (eta_bin[0], eta_bin[1])
        title_str += r"(%s %s MC)" % (year, sample_str)
        
        ax.set_title(title_str)
        ax.legend(loc="upper left", framealpha=1)

        # Propagate uncertainties.
        #----- Overriding ratio_numer_errs with Andrey's formula. -----#
        # [X] Confirm that this formula is the same as mine.
        # It's not the same, but his uncertainties are smaller. Use his. 
        ratio = (sigma_corr_arr - sigma_arr) / sigma_arr
        # ratio = np.true_divide(ratio_numer, ratio_denom, out=np.zeros_like(sigma_arr), where=ratio_denom!=0)
        ratio_err = Andrey_prop_err_on_dsigoversig(sigma_arr, sigma_err_arr, sigma_corr_arr, sigma_corr_err_arr)

        ax_ratio = fig.add_subplot(gs[1])
        ax_ratio.errorbar(x, ratio, xerr=[x_err_low, x_err_high], yerr=ratio_err, 
                        capsize=capsize, fmt=markerstyle, 
                        mfc=color_sigma_corr, mec=color_sigma_corr, ecolor=color_sigma_corr, 
                        elinewidth=errlw, mew=markererrw)
        #                     color=color_dict[count], 
        #                         , ms=ms, 
        #                         , ,
        if (set_log_x):
            ax_ratio.set_xscale('log')
        ax_ratio.set_xlim([x_min, x_max])
        ax_ratio.set_xlabel(r'$\mathrm{p_T}$ bin edges [GeV]')
        ax_ratio.set_ylabel(r'$\frac{\sigma^{\mathrm{corr.}} - \sigma}{\sigma}$')

        ax_ratio.set_ylim([-0.25,0.25])
        ax_ratio.axhline(y=0, xmin=-0.1, xmax=1.1, color='k', lw=0.7)#, marker=None)

        ax.xaxis.set_major_formatter(plt.NullFormatter())
        ax_ratio.grid(b=None)
        pdf.savefig()
        plt.clf()
    # End eta_bin loop.
print_header_message(f"PDF succesfully made at:\n{outpath_pdf}")