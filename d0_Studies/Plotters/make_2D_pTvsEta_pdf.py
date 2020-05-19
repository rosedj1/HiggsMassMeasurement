import os
import sys
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

import matplotlib.pyplot as plt
import pandas as pd
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly
from d0_Utils.d0_cls import KinematicBin
from matplotlib.backends.backend_pdf import PdfPages

plt.style.use('cmsstyle_plot')

n_evts_scan = 10000
n_evts_keep = 500

massZ_cut_ls = [60,120]
d0q_cut_ls = [-0.005, 0.005]

# eta_bin_edges = [0.00, 0.30] # barrel
# eta_bin_edges = [0.80, 1.10]  # overlap
eta_bin_edges = [2.10, 2.20, 2.40]  # endcap

# eta_bin_edges = [0.00, 0.10, 0.20, 0.30]    # barrel
# eta_bin_edges = [0.70, 0.80, 0.90, 1.00, 1.10, 1.20]    # overlap
# eta_bin_edges = [2.00, 2.10, 2.20, 2.30, 2.40]    # endcap

#--- Could be either pT or p. User specifies when initializing kbin. ---#
# p_bin_edges = [5, 20, 30, 40, 50, 60, 100]
# p_bin_edges = [5, 7, 10, 15, 20, 25, 30, 35, 40, 45, 50, 100]
# p_bin_edges = [5, 7, 10, 15, 20]
p_bin_edges = [20, 40]#, 60]
d0_type = 'BS'
dR_cut = 0.02
verbose = True
use_ptotal_instead = False  # p_total or pT

x_bin_limits = [-2.5, 2.5, 0.1] 
y_bin_limits = [0, 100, 1] 

outpath = "/Users/Jake/Desktop/20200427/2Dplots_pTvsEta/"
pdf_name_base = "test2"  # Does not need filler underscores.

df = pd.read_hdf('/Users/Jake/Desktop/MC_2016.h5')

#--- Main ---#
d0q_min = d0q_cut_ls[0]
d0q_max = d0q_cut_ls[1]

makeDirs(outpath)
all_kbin_ls = []
for k in range(len(eta_bin_edges)-1):
    this_eta = eta_bin_edges[k]
    next_eta = eta_bin_edges[k+1]
    
    kbin_ls = []
    for m in range(len(p_bin_edges)-1):
        this_p = p_bin_edges[m]
        next_p = p_bin_edges[m+1]
        
        kbin = KinematicBin(df, 
                            n_evts=n_evts_scan, 
                            massZ_cut_ls=massZ_cut_ls,
                            eta_cut_ls=[this_eta, next_eta], 
                            pT_cut_ls=[this_p, next_p], 
                            d0q_cut_ls=d0q_cut_ls,
                            d0_type=d0_type, 
                            dR_cut=dR_cut, 
                            use_ptotal_instead=use_ptotal_instead, 
                            verbose=verbose)
        kbin_ls.append(kbin)
        
    # Finished looping over p bins.
    evts_max = n_evts_keep

    # Save plots into one PDF:
    pT_min = min( [kb.pT_min for kb in kbin_ls] )
    pT_max = max( [kb.pT_max for kb in kbin_ls] )
    n_plots = len(kbin_ls)

    title_str_pT_min = f"{pT_min}" if pT_min < 10 else f"{pT_min}"  # For plot-ordering purposes.
    pdf_name = pdf_name_base + f"__{this_eta}_eta_{next_eta}"# + f"__{title_str_pT_min}_pT_{pT_max}"
    pdf_name += f"__{pT_min}_pT_{pT_max}"# + f"__{title_str_pT_min}_pT_{pT_max}"
    pdf_name += f"__{d0q_min}_d0q{d0_type}_{d0q_max}"
    pdf_name += f"__{n_plots}plots"
    pdf_name = make_str_title_friendly(pdf_name) + ".pdf"

    outfile = os.path.join(outpath, pdf_name)

    with PdfPages(outfile) as pdf:
        for kbinned in kbin_ls:
            kbinned.make_2D_plot(x_kinem="eta", y_kinem="pT", 
                                 x_bin_limits=x_bin_limits, y_bin_limits=y_bin_limits,
                                 lep_selection_type="independent",
                                 run_over_only_n_evts=evts_max, 
                                 title="",
                                 ax=None,
                                 verbose=verbose)

            pdf.savefig()  # saves the current figure into a pdf page
            plt.close()

    print("[INFO] PDF created at", outfile, "\n")

    plt.close('all')