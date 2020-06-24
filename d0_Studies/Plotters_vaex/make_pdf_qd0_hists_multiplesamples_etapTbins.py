"""
# PURPOSE: 
#    This script makes 1 PDF of qd0 distributions. 
#    It can probably easily accomodate other kinematic variables, like dpT/pT.
# NOTES: 
#    
# 
"""
import time
import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt

from Utils_vaex.vaex_fns import vaex_apply_masks, prepare_vaex_df
from Samples.sample_info import Sample
from d0_Studies.KinBin_Info.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum,
                                                   bin_edges_pT_sevenfifths_to1000GeV_wholenum)
from d0_Utils.d0_dicts import label_LaTeX_dict
from d0_Utils.d0_fns import make_binning_array

from Utils_Python.Utils_Files import makeDirs, make_str_title_friendly, check_overwrite
from Utils_Python.Utils_Plotting import hist_y_label, make_1D_dist, ncolsrows_from_nplots

from matplotlib.backends.backend_pdf import PdfPages

#---------------------------#
#----- User Parameters -----#
#---------------------------#
outdir = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/hists_qd0/"
filename_base = "qd0dist_etapTbins"
sample_ls = [
    Sample("MC", "2017", "Jpsi", "hdf5"),
    Sample("MC", "2017", "DY", "hdf5"),
]   
overwrite = False

# Binning.
eta_ls = equal_entry_bin_edges_eta_mod1_wholenum
pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum
qd0_bin_info = [-0.020, 0.020, 0.0002]  # min, max, bin_width.
x_range = [-0.020, 0.020]

#-------------------------------------#
#----- Script-Specific Functions -----#
#-------------------------------------#
def check_sample_ls(sample_ls):
        yr_ls = []
        type_ls = []
        for smpl in sample_ls:
            assert smpl.name in ["DY", "Jpsi"]
            yr_ls.append(smpl.year)
            type_ls.append(smpl.data_type)
        assert len(set(yr_ls)) == 1
        assert len(set(type_ls)) == 1

def make_prefix(sample_ls):
    name_str = ""
    for smpl in sample_ls:
        name_str += smpl.name
    return f"{sample_ls[0].data_type}{sample_ls[0].year}{name_str}"

def make_suffix(eta_ls, pT_ls):
    suffix = (
        f"{min(eta_ls):.1f}eta{max(eta_ls):.1f}"
        f"_{min(pT_ls):.1f}pT{max(pT_ls):.1f}GeV"
        )
    suffix = make_str_title_friendly(suffix)
    return suffix

def make_full_filename(sample_ls, eta_ls, pT_ls, filename_base):
    prefix = make_prefix(sample_ls)
    suffix = make_suffix(eta_ls, pT_ls)
    return f"{prefix}_{filename_base}_{suffix}.pdf"

def prep_area(outdir, outfile_name, overwrite):
    outfile_path = os.path.join(outdir, outfile_name)
    check_overwrite(outfile_path, overwrite=overwrite)
    makeDirs(outdir)
    return outfile_path

def prep_plots(sample_ls, eta_ls, pT_ls):
    n_pages = len(eta_ls)-1
    n_plots = len(pT_ls) - 1
    rows, cols = ncolsrows_from_nplots(n_plots)
    print(f"Making a {n_pages}-page PDF.")
    print(f"Making {n_plots} plots per page ({rows} x {cols} grid).\n")
    print(f"|eta| regions:\n{np.round(eta_ls, decimals=2)}")
    print(f"pT regions:\n{np.round(pT_ls, decimals=2)}\n")
    for smpl in sample_ls:
        print(f"   ...Running over {smpl.data_type} {smpl.year} {smpl.name} sample...")
    return n_pages, n_plots, rows, cols

def make_cut_str(sample_ls, eta_range, pT_range, qd0_range):
    """Return a string of cuts to show on plot."""
    eta_min, eta_max = eta_range[0], eta_range[1]
    pT_min, pT_max = pT_range[0], pT_range[1]
    qd0_min, qd0_max = qd0_range[0], qd0_range[1]

    cuts  = r"$%.2f < \left| \eta^{\mathrm{REC}} \right| < %.2f$," % (eta_min, eta_max) + "\n"
    cuts += r"$%.1f < p_{T} < %.1f$ GeV," % (pT_min, pT_max) + "\n"
    cuts += r"$%.4f < q(\mu)*d_{0}^{\mathrm{%s}} < %.4f$, " % (qd0_min, "BS", qd0_max) + "\n"

    for smpl in sample_ls:
        # Sample-specific items.
        LaTeX_str = smpl.LaTeX_name
        m_inv_cut_min = smpl.inv_mass_cut_lim[0]
        m_inv_cut_max = smpl.inv_mass_cut_lim[1]
        dR_max = smpl.dR_cut
        cuts += r"$%.1f < m_{%s} < %.1f$ GeV,  " % (m_inv_cut_min, LaTeX_str, m_inv_cut_max)
        cuts += r"$\Delta R < %.3f$" % (dR_max) + "\n"
    return cuts

def make_sample_masks(sample_ls, eta_range, pT_range, qd0_range):
    """
    Create masks for eta, pT, qd0 cuts. 
    Also create masks for m_mumu and dR cuts based on sample name.
    """
    for smpl in sample_ls:
        massZ_minmax = smpl.inv_mass_cut_lim
        dR_cut = smpl.dR_cut
        smpl.mask = vaex_apply_masks(smpl.vdf_prepped, eta_range, pT_range, qd0_range, massZ_minmax, dR_cut)

def apply_masks(sample_ls, kinem):
    data = []
    for smpl in sample_ls:
        vals = smpl.vdf_prepped[kinem].values  # Get col as numpy.ndarray.
        vals_after_cut = vals[smpl.mask.values]  # Apply mask. 
        data.extend(vals_after_cut)
    return np.asarray(data)

def make_ax_labels(kinem, bin_width):
    """Make str labels for x- and y-axes."""
    x_label = label_LaTeX_dict[kinem]["independent_label"]
    x_units = label_LaTeX_dict[kinem]["units"]
    y_label = hist_y_label(bin_width, x_units)
    if len(x_units) > 0:
        x_label += " [{}]".format(x_units)
    return x_label, y_label

if __name__ == "__main__":
    # Prep work.
    check_sample_ls(sample_ls)
    outfile_name = make_full_filename(sample_ls, eta_ls, pT_ls, filename_base)
    outfile_path = prep_area(outdir, outfile_name, overwrite)

    # Plot info.
    plt.style.use("/ufrc/avery/rosedj1/HiggsMassMeasurement/Utils_Python/Plot_Styles_Matplotlib/grid_multiple_plots.mplstyle")
    n_pages, n_plots, rows, cols = prep_plots(sample_ls, eta_ls, pT_ls)
    qd0_range = qd0_bin_info[0:2]
    x_bins, bin_width = make_binning_array(qd0_bin_info)
    
    # Analysis.
    # Make only 1 PDF. Each page is a different eta cut.
    with PdfPages(outfile_path) as pdf:
        total_entries = 0
        # Loop over eta regions. 
        for page, (eta_min, eta_max) in enumerate(zip(eta_ls[:-1], eta_ls[1:]), 1):
            eta_range = [eta_min, eta_max]
            
            # Within this eta region, scan the pT regions. 
            t_start = time.perf_counter()
            f = plt.figure()
            for count, (pT_min,pT_max) in enumerate(zip(pT_ls[:-1], pT_ls[1:])):
                pT_range = [pT_min, pT_max]
                ax = plt.subplot(rows,cols,count+1)
                x_label, y_label = make_ax_labels("qd0BS1", bin_width)

                make_sample_masks(sample_ls, eta_range, pT_range, qd0_range)
                ax, bin_vals, bin_edges, stats = make_1D_dist(ax=ax, 
                                                              data=apply_masks(sample_ls, "qd0BS"),
                                                              x_limits=x_range,
                                                              x_bin_edges=x_bins, 
                                                              x_label=x_label, 
                                                              y_label=y_label,
                                                              title="",
                                                              y_max=-1,
                                                              log_scale=False, color=None, leg_loc=None, display="sci")
                
                cut_str = make_cut_str(sample_ls, eta_range, pT_range, qd0_range)
                ax.text(0.025, 0.78, cut_str, horizontalalignment='left', verticalalignment='center', transform=ax.transAxes,
                            bbox=dict(boxstyle='square', facecolor='white', alpha=0.9))
                
                t_end = time.perf_counter()
                total_entries += stats[0]
            # End pT loop.
            
            plt.tight_layout()
            pdf.savefig()
            plt.close("all")
            print(f"Page {page}/{n_pages} made. Time taken: {t_end - t_start:.2f} s")

        # End eta loop.
        print(f"PDF made at:\n  {outfile_path}")

    total_muons_original = 0
    for smpl in sample_ls:
        total_muons_original += smpl.vdf_prepped.count()
    print(f"Total muons in all samples: {total_muons_original}")
    print(f"Total muons found after cuts: {total_entries}")