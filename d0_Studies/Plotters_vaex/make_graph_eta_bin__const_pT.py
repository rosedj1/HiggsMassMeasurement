"""
# PURPOSE: 
#  asdf
# NOTES:
#  asdf   
# SYNTAX:  python script.py
# AUTHOR:  Jake Rosenzweig
# UPDATED: 2020-05-18
"""
import sys
import pickle

import matplotlib.pyplot as plt
import numpy as np
plt.style.use("cmsstyle_plot")
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')

from d0_Utils.d0_cls import KinBin3DOrganizer, GraphLineKinBin3D
from d0_Utils.d0_dicts import label_LaTeX_dict, color_dict
from PyUtils.Utils_StatsAndFits import prop_err_x_div_y
from scipy.optimize import curve_fit

fullpath_kinbin_ls_pkl = "/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/d0_studies/kinbin3D_pkls/fullscan_12reg__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl"

with open(fullpath_kinbin_ls_pkl, "rb") as f:
    kinbin_ls = pickle.load(f)

kb_org = KinBin3DOrganizer(kinbin_ls)

ls_of_kb_ls = kb_org.find_all_KinBin_ls_constant_range(const_range=[10.0, 14.0], const_bin="pT")

fig, ax = plt.subplots(1,1, sharey=False)
for count, kb_ls in enumerate(ls_of_kb_ls[:5], 1):

    qd0_avg_ls, dpToverpT_bestfitmean_ls, dpToverpT_bestfitmean_err_ls = kb_org.get_plotting_vals_from_KinBin_ls(kb_ls)

    line = GraphLineKinBin3D(x_vals=qd0_avg_ls, 
                  y_vals=dpToverpT_bestfitmean_ls, 
                  x_err_vals=None,
                  y_err_vals=dpToverpT_bestfitmean_err_ls)
    line.draw_graph(x_label="",
                     y_label="", 
                     title="", 
                     kbin_example=kb_ls[0], 
                     ax=ax, count=count,
                     scale_by_1divpT=False,
                     verbose=False,
                     constant_bin="pT")

plt.tight_layout()
plt.savefig("/Users/Jake/Desktop/graph_eta_inclus__pT_10p0to14p0_test03.pdf")





# kb_ls = kinbin_dict[2.3][10.0]

# qd0_avg_ls, dpToverpT_bestfitmean_ls, dpToverpT_bestfitmean_err_ls = kb_org.get_plotting_vals_from_KinBin_ls(kb_ls)

# line = GraphLineKinBin3D(x_vals=qd0_avg_ls, 
#                   y_vals=dpToverpT_bestfitmean_ls, 
#                   x_err_vals=None,
#                   y_err_vals=dpToverpT_bestfitmean_err_ls)

# fig, ax = plt.subplots(1,2, sharey=False, figsize=(12.8, 9.6))

# line.draw_graph(x_label="",
#                  y_label="", 
#                  title="", 
#                  kbin_example=kb_ls[0], 
#                  ax=ax[0], count=1,
#                  scale_by_1divpT=False,
#                  verbose=False,
#                  constant_bin="pT")

# line.draw_graph(x_label="",
#                  y_label="", 
#                  title="", 
#                  kbin_example=kb_ls[0], 
#                  ax=ax[1], count=2,
#                  scale_by_1divpT=True,
#                  verbose=False, 
#                 constant_bin="pT")

# plt.tight_layout()
# plt.savefig("/Users/Jake/Desktop/barrel_test06.pdf")