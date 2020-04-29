import os
import sys
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')
from PyUtils.Utils_Files import makeDirs
from matplotlib.backends.backend_pdf import PdfPages

from d0_Utils.d0_cls import KinematicBin, KinBinOrganizer, GraphLine
from d0_Utils.d0_fns import make_binning_array, shift_binning_array, calc_x_err_bins
from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict

plt.style.use('cmsstyle_plot')
#--------------------------------------------------------------#
n_evts_scan = -1
n_evts_keep = -1

dataframe = '/Users/Jake/Desktop/MC_2016.h5'
outpath = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/20200427/plot_dpToverpT_vs_d0q__all_eta_bins.pdf"

# d0_bin_limits = [-0.006, 0.006, 0.002]
d0_bin_limits = [-0.006, 0.006, 0.001]
massZ_cut_ls=[60,120]
pT_cut_ls=[5,100]
eta_cut_ls=[
    [0.0, 0.4],
    [0.4, 0.8],
    [0.8, 1.2],
    [1.2, 1.6],
    [1.6, 2.0],
    [2.0, 2.4],
]
dR_cut=0.02

d0_type='BS'
fit_iters = 3
use_ptotal_instead=False
verbose = True
#---------- Automatons ----------#
df = pd.read_hdf(dataframe)

dir_ = os.path.dirname(outpath)
makeDirs(dir_) 

org_kbin_ls = []
with PdfPages(outpath) as pdf:
    for eta_reg_ls in eta_cut_ls: 
        org_kbin = KinBinOrganizer(df, n_evts_scan, massZ_cut_ls, eta_reg_ls, pT_cut_ls, d0_type, dR_cut, use_ptotal_instead, verbose)
        org_kbin.make_kbin_ls_over_d0_range(d0_bin_limits)
        org_kbin.plot_dpToverpT_for_kbin_ls(kinem="delta_pToverpT1", lep_selection_type='independent', 
                                           x_limits=[-0.3, 0.3], bin_limits=[-0.3, 0.3, 0.004], 
                                           run_over_only_n_evts=n_evts_keep, ax=None, y_max=-1, log_scale=False, 
                                           iter_gaus=(True, fit_iters),
                                           make_pdf=True,
                                           pdf_obj=pdf    # If make_pdf=True, must pass in pdf object here!
                                           )  
        org_kbin.get_dpToverpT_iter_gaus_fit_means()

        org_kbin_ls.append(org_kbin)
        plt.close("all")
    # Done looping over all eta regions.
    
    # Plot all graph lines on single axes. 
    graph_ls = []
#     f, ax = plt.subplots(figsize=(12.8, 9.6))
    f, ax = plt.subplots()
    for count,org_kb in enumerate(org_kbin_ls, 1):
        graph = GraphLine(org_kb.d0_bin_arr_shifted, org_kb.fit_mean_ls, org_kb.fit_mean_err_ls)
        graph.draw_graph("d0BSq1", "delta_pToverpT1", binning_type="eta", kbin_example=org_kb.kbin_ls[0], ax=ax, count=count)
        # Note: if you want stats to show up on plot when doing draw_graph, then make sure count=1 at some point.
        graph_ls.append(graph)
    # Done looping over KinBinOrganizers.
    plt.tight_layout()
    pdf.savefig()
    plt.close("all")