# PURPOSE: Make distributinos of pT^(rec) in whatever eta bins User specifies. 
# NOTES:   User can also specify 
# SYNTAX:  python <script>.py
# AUTHOR:  Jake Rosenzweig
# DATE:    2020-05-01

import os
import sys
import matplotlib.pyplot as plt
import pandas as pd

sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')
from PyUtils.Utils_Files import makeDirs, check_overwrite
from matplotlib.backends.backend_pdf import PdfPages
from d0_Utils.d0_cls import KinBinOrganizer

plt.style.use('cmsstyle_plot')
#--------------------------------------------------------------#
n_evts_scan = 3000000
n_evts_keep = 500000
#n_evts_scan = 10000
#n_evts_keep = 5000

dataframe = '/Users/Jake/Desktop/MC_2016.h5'
outpath = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/pT_distributions_eta_binning/test_07.pdf"
overwrite = True

d0_bin_limits = [-1, 1, 2]
massZ_cut_ls=[60,120]
pT_cut_ls=[5,1000]
eta_cut_ls=[
    [0.0, 0.2],
    [0.2, 0.4],
    [0.4, 0.6],
    [0.6, 0.8],
    [0.8, 1.0],
    [1.0, 1.2],
    [1.2, 1.4],
    [1.4, 1.6],
    [1.6, 1.8],
    [1.8, 1.9],
    [1.9, 2.0],
    [2.0, 2.1],
    [2.1, 2.2],
    [2.2, 2.3],
    [2.3, 2.4],
]

binning_type = "eta"
d0_type="BS"
dR_cut=0.008
use_ptotal_instead = False
verbose = True
iter_gaus = (False, 3)

x_limits = [-40, 250]
bin_limits = [5, 200, 1] 
#---------- Automatons ----------#
df = pd.read_hdf(dataframe)

dir_ = os.path.dirname(outpath)
check_overwrite(outpath,overwrite)
makeDirs(dir_) 

with PdfPages(outpath) as pdf:
    for eta_reg_ls in eta_cut_ls: 
        org_kbin = KinBinOrganizer(df, n_evts_scan, massZ_cut_ls, eta_reg_ls, pT_cut_ls, 
                                   d0_type, dR_cut, use_ptotal_instead, verbose)
        org_kbin.make_kbin_ls_over_d0_range(d0_bin_limits)
        org_kbin.plot_dpToverpT_for_kbin_ls(kinem="pT1", lep_selection_type='independent', 
                                           x_limits=x_limits, bin_limits=bin_limits,
                                           run_over_only_n_evts=n_evts_keep, ax=None, y_max=-1, log_scale=False, 
                                           iter_gaus=iter_gaus,
                                           make_pdf=True,
                                           pdf_obj=pdf    # If make_pdf=True, must pass in pdf object here!
                                           )  
        plt.close("all")
