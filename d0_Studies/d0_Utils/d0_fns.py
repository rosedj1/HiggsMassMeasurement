import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from matplotlib.backends.backend_pdf import PdfPages

#from d0_Utils.d0_cls import HistInfo    # Not sure why I can't import this...

def get_subset_mask(x_vals, x_min, x_max):
    """
    Return a mask of an array such that: x_min <= element <= x_max
    Useful for selecting fit ranges along the 'x-axis'.
    
    Apply the mask: x_vals[mask]
    
    Parameters
    ----------
    x_vals : array
        The original array from which the subset will be made.
    x_min : float
        The minimum value in x_vals to be kept.
    x_max : float
        The maximum value in x_vals to be kept.
        
    Returns
    -------
    mask : bool array
        An array of booleans. If an element is 'True' then: x_min <= element <= x_max.
        The mask will be the same length as x_vals. 
    """
    mask = (x_min <= x_vals) & (x_vals <= x_max)
    
    return mask

def combine_cut_list(cut_list):
    """Used for ROOT-style cuts while doing tree.Draw()"""
    total_cuts = ''
    for cut in cut_list:
        total_cuts += cut + ' && '
    cut_str = total_cuts.rstrip(' && ')
    return cut_str

def calc_num_bins(bin_min, bin_max, bin_width):
    """
    Calculates the number of bins (int) given: bin_min, bin_max, bin_width

    Parameters
    ----------
    bin_min, bin_max, bin_width
    """
    return int(round( (bin_max-bin_min)/bin_width ))

def calc_x_err_bins(x_val_center_list):
    """
    Returns lists of the x-errors, which may be symmetrical or asymmetrical. 
    """
    high_err_list = []
    for x_center in range(len(x_val_center_list)-1):
        high_err = float( x_val_center_list[x_center+1] - x_val_center_list[x_center] ) / 2
        high_err_list.append(high_err)
        
    low_err_list = [high_err_list[0]] + high_err_list
    high_err_list = high_err_list + [high_err_list[-1]]
    
    return low_err_list, high_err_list

def calc_ymin_for_legend(n_graphs, text_height=0.042):
    """
    ROOT specific! Allows for dynamic expansion of bottom axis of TLegend. 
    There is probably an in-built function, but I want to write my own!
    
    Parameters
    ----------
    n_obj : int
        Number of objects (like graphs, histos) that will show up in your legend.
    text_height : float between 0 and 1
        Essentially sets the font. Was 0.033333 but you'll have to play with it. 
    """
    delta_y = text_height*n_graphs
    return 0.9 - delta_y # When you do TCanvas(), it puts the top axis at y = 0.9.

def make_deltapT2_pdf(df, year, sample, 
                     bspv, 
                     d0_range_ls, 
                     deltapT_range_ls, 
                     pT_cut_ls, 
                     eta_cut_ls, 
                     outfile_path,
                     wrt="pT_gen",
                     save_plots=False, overwrite=False):

    """
    Takes in a DataFrame (df) corresponding to an entire '.root' file and generates a '.pdf' with 
    many deltapT/pT^2 histograms, one histo for each d0 bin. 
    Each histo has the same cuts based on:
        bspv : 'BS' or 'PV'
        deltapT_range_ls : [deltapT_min, deltapT_max, deltapT_bin_width]
        pT_cut_ls : [pT_min, pT_max]
        eta_cut_ls : [eta_min, eta_max]
        
    The df you pass in has a: 
        year : 
        sample : 'MC' or 'Data'
    
    #--- Details ---#
    Go through each event and see if mu1 passes all selection criteria.
    Demand that d0BS1*charge(mu1) fall within specified d0' bin
    Make a DataFrame with all events that pass these selections.
    Do the same thing for mu2 and make another DataFrame.
    
    Returns
    -------
    A list of HistInfo objects. 
    Each HistInfo object contains all the important information from each deltapT histogram. 
    
    If save_plots=True, then will make a '.pdf' showing each deltapT histogram. 
    """

    textsize_legend = 8
    textsize_axislabels = 10
    textsize_ticklabels = 8
    textsize_title = 10
    
    d0_min = d0_range_ls[0]
    d0_max = d0_range_ls[1]
    d0_bin_width = d0_range_ls[2]

    if d0_bin_width < 0.005:
        err_msg = (f"WARNING: d0_bin_width ({d0_bin_width}) is too small (d0_bin_width < 0.005). Stopping now.\n"
                   "To get around this, change the value in the function: make_deltapT2_pdf()")
        raise ValueError(err_msg)    
    deltapT_min = deltapT_range_ls[0]
    deltapT_max = deltapT_range_ls[1]    
    deltapT_bin_width = deltapT_range_ls[2]    
    
    pT_min = pT_cut_ls[0]
    pT_max = pT_cut_ls[1]
    
    eta_min = eta_cut_ls[0]
    eta_max = eta_cut_ls[1]
    
    massZ_min = 70
    massZ_max = 110
    
    d0_bin_arr = np.arange(d0_min, d0_max+0.5*d0_bin_width, d0_bin_width)  # Includes all bin edges: very first to very last.
    d0_bin_arr_shifted = d0_bin_arr[0:-1] + 0.5*d0_bin_width  # Make sure points are plotted in middle of bin window.
    # d0_n_bins = calc_num_bins(d0_min, d0_max, d0_bin_width)
    deltapT_arr = np.arange(deltapT_min, deltapT_max+0.5*deltapT_bin_width, deltapT_bin_width)
    
    pdf_name  = f"{year}_{sample}_{bspv}__"
    title_str_pT_min = f"0{pT_min}" if pT_min < 10 else f"{pT_min}"  # For plot-ordering purposes.
    # Example of Python magic: auto-concatenation of strings.
    pdf_name += (f"{title_str_pT_min}_pT_{pT_max}"
                 f"__{eta_min}_eta_{eta_max}"
                 f"__{d0_min:.3f}_to_{d0_max:.3f}_increm_{d0_bin_width:.3f}"
                 f"__wrt_{wrt}")
    pdf_name = make_str_title_friendly(pdf_name) + ".pdf"
    
    outfile = os.path.join(outfile_path, pdf_name)
    if not os.path.exists(outfile_path):
        os.makedirs(outfile_path)
    if os.path.exists(outfile) and not (overwrite):
        print(f"Skipping {outfile} since it already exists.\nTo write over the file then set overwrite=True.\n")
        return
        
    status = f"Running over: {year} {sample} {bspv}, pT_range={pT_cut_ls}, eta_range={eta_cut_ls}, wrt {wrt}"
    print(status)

    with PdfPages(outfile) as pdf:
        charge1_ser = df['Id1'].replace(13,-1).replace(-13,1)
        charge2_ser = df['Id2'].replace(13,-1).replace(-13,1)

        d0BS1xcharge_ser = charge1_ser * df['d0BS1']
        d0BS2xcharge_ser = charge2_ser * df['d0BS2']
        d0PV1xcharge_ser = charge1_ser * df['d0PV1']
        d0PV2xcharge_ser = charge2_ser * df['d0PV2']  
        
        graph_info_ls = []
        
        for k in range(len(d0_bin_arr)-1):
            this_d0_bin = d0_bin_arr[k]
            next_d0_bin = d0_bin_arr[k+1]
            if bspv in 'BS':
                mask_d0BS1xcharge = (this_d0_bin < d0BS1xcharge_ser) & (d0BS1xcharge_ser < next_d0_bin)
                mask_d0BS2xcharge = (this_d0_bin < d0BS2xcharge_ser) & (d0BS2xcharge_ser < next_d0_bin)
            elif bspv in 'PV':
                mask_d0PV1xcharge = (this_d0_bin < d0PV1xcharge_ser) & (d0PV1xcharge_ser < next_d0_bin)
                mask_d0PV2xcharge = (this_d0_bin < d0PV2xcharge_ser) & (d0PV2xcharge_ser < next_d0_bin)

            #--- Create masks ---#
            mask_massZ = (massZ_min < df['massZ']) & (df['massZ'] < massZ_max)
            # new_df = df[massZ_min < df['massZ']]    # Turns out to not be efficient to trim DataFrame as you go...
                        
            mask_pT1 = (pT_min < df['pT1']) & (df['pT1'] < pT_max)
            mask_pT2 = (pT_min < df['pT2']) & (df['pT2'] < pT_max)

            mask_eta1 = (eta_min < df['eta1']) & (df['eta1'] < eta_max)
            mask_eta2 = (eta_min < df['eta2']) & (df['eta2'] < eta_max)
            
            mask_skip_empty_GENmass2l = 0 < df['GENmass2l']
            
            # Apply all masks.
            if bspv in 'BS':
                original_masked_df1 = df[(mask_d0BS1xcharge & mask_skip_empty_GENmass2l & mask_pT1 & mask_eta1 & mask_massZ)]
                original_masked_df2 = df[(mask_d0BS2xcharge & mask_skip_empty_GENmass2l & mask_pT2 & mask_eta2 & mask_massZ)]
            elif bspv in 'PV':
                original_masked_df1 = df[(mask_d0PV1xcharge & mask_skip_empty_GENmass2l & mask_pT1 & mask_eta1 & mask_massZ)]
                original_masked_df2 = df[(mask_d0PV2xcharge & mask_skip_empty_GENmass2l & mask_pT2 & mask_eta2 & mask_massZ)]

            if wrt in "pT_gen":
                deltapT1_over_pTsqrd_ser = (original_masked_df1.pT1 - original_masked_df1.genLep_pt1) / original_masked_df1.genLep_pt1**2
                deltapT2_over_pTsqrd_ser = (original_masked_df2.pT2 - original_masked_df2.genLep_pt2) / original_masked_df2.genLep_pt2**2
            elif wrt in "pT_reco":
                deltapT1_over_pTsqrd_ser = (original_masked_df1.pT1 - original_masked_df1.genLep_pt1) / original_masked_df1.pT1**2
                deltapT2_over_pTsqrd_ser = (original_masked_df2.pT2 - original_masked_df2.genLep_pt2) / original_masked_df2.pT2**2

            # Muon1 has been checked for each event and one series was made.
            # Then muon2 was checked and a series was made. Combine both series.
            combined_deltapT_ser = deltapT1_over_pTsqrd_ser.append( deltapT2_over_pTsqrd_ser )
            
            # Example of Python magic. Auto-concatenation of strings.
            cuts = (f"{this_d0_bin:.3f}" + r"$<d_{0}^{%s}*\mathrm{charge}(\mu)<$"%bspv + f"{next_d0_bin:.3f}, "
                    "\n"
                    f"{pT_min}<" + r"$p_{T}^{RECO}$" + f"<{pT_max} GeV, "
                    f"{eta_min}<" + r"$\left| \eta \right|$<" + f"{eta_max}, "
                    f"{massZ_min} GeV<" + r"$m_{\mu\mu}$")

            # Statistics.
            n_entries = len(combined_deltapT_ser)
            mean = combined_deltapT_ser.mean()
            mean_err = abs(mean) / np.sqrt(n_entries)  # Standard error of the mean.
            stdev = combined_deltapT_ser.std()
            stdev_err = stdev / np.sqrt(2*n_entries)

            # Store info of this plot in a HistInfo object. 
            this_hist = HistInfo()            
            this_hist.year = year
            this_hist.sample = sample
            this_hist.bspv = bspv
            this_hist.d0_bin_window = [this_d0_bin, next_d0_bin]
            this_hist.x_axis_bounds_list = deltapT_range_ls[0:2] # 2-element list
            this_hist.eta_range = eta_cut_ls
            this_hist.pT_range = pT_cut_ls
            this_hist.massZ_cut = massZ_min
            this_hist.n_entries = n_entries
            this_hist.hist_mean = mean
            this_hist.hist_mean_err = mean_err
            this_hist.hist_stdev = stdev  # Spread of the data.
            this_hist.hist_stdev_err = stdev_err
            
            graph_info_ls.append(this_hist)
            
            # Plotting
            fig, ax = plt.subplots()
            ax.tick_params(axis='both', which='major', direction='in', length=9,   top=True, right=True)
            ax.tick_params(axis='both', which='minor', direction='in', length=4.5, top=True, right=True)
            if wrt in "pT_gen":
                ax.set_xlabel(r"$(p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{GEN})^2$ [GeV]$^{-1}$", fontsize=textsize_axislabels)
            if wrt in "pT_reco":
                ax.set_xlabel(r"$(p_{T}^{RECO} - p_{T}^{GEN}) / (p_{T}^{RECO})^2$ [GeV]$^{-1}$", fontsize=textsize_axislabels)
            ax.set_ylabel(f"Events / [{deltapT_bin_width:.1E} GeV" + r"$^{-1}$]", fontsize=textsize_axislabels)
            ax.set_title(cuts, fontsize=textsize_title)
            ax.set_xlim([deltapT_min, deltapT_max])
            
            # Add a x10^ scale factor to the x-axis.
            formatter = ticker.ScalarFormatter(useMathText=True)
            formatter.set_powerlimits((-1,1))
            ax.xaxis.set_major_formatter(formatter)
            
            ax.grid(which='major',color='k', ls=':')
            
            plt.minorticks_on()

            plt.rc('axes', titlesize=textsize_title)     # fontsize of the axes title
            plt.rc('legend', fontsize=textsize_legend)    # legend fontsize
            
            
            plt.hist(combined_deltapT_ser, bins=deltapT_arr, histtype='step', color='blue')
            # deltapT1_over_pTGENsqrd_BS_ser.plot.hist(bins=100)

            # Stats box.
            textstr = '\n'.join((
                f'{year} {sample}',
                f'Entries = {n_entries}',
                f'Mean = {mean:.2E}' + r' $\pm$ ' + f'{mean_err:.2E}',
                f'Std Dev = {stdev:.2E}' + r' $\pm$ ' + f'{stdev_err:.2E}'))
            props = dict(boxstyle='square', facecolor='whitesmoke')
            ax.text(0.9915, 0.986, textstr, transform=ax.transAxes, fontsize=textsize_legend,
                    verticalalignment='top', 
                    horizontalalignment='right',
                    bbox=props)

            if (save_plots): 
                pdf.savefig()  # saves the current figure into a pdf page

            plt.close()

        print(f"PDF saved at:\n{outfile}\n")

        return graph_info_ls

def make_graph(*graph_tuple, binning_type='eta', verbose=False, save_plots=True, outpath='/Users/Jake/Desktop/'):
    """
    Make one deltapT vs. d0*charge plot, with as many graphs displayed as len(graph_tuple).
    graph_tuple contains many lists (graph_lsX) of histograms (HistInfo).
    
    graph_tuple = (
        graph_ls1,  # Each hist in this graph_ls all have identical cuts, but are made over different d0 bins.
        graph_ls2,  # This graph_ls2 may have slightly different cuts than graph_ls1.
        ...,
    )
    
    #--- Example ---#
    graph_ls1 = [HistInfo1, HistInfo2, ...]  # This list will generate 1 deltapT_vs_d0charge graph. 
    graph_ls2 = [HistInfo1, HistInfo2, ...]  # This list will generate another graph, but may have different cuts. 

    If binning_type is 'eta', then each graph_ls in graph_tuple should have the same pT cuts.
        This way, we can view different eta ranges while making sure pT cuts are all the same. 
    If binning_type is 'pT', then each graph_ls in graph_tuple should have the same eta cuts. 
    
    Returns
    -------
    graph_ls contains all the correct graph info to make a plot. 
    """    
    example_graph_ls = graph_tuple[0]
    example_graph = graph_tuple[0][0] # First graph in first graph list.

    # Put all HistInfo objects into a single list.
    entire_HistInfo_list = []
    for gr_ls in graph_tuple:
        entire_HistInfo_list += gr_ls
        
    #--- Check that things make sense: ---#
    # Example: If binning in eta, make sure each HistInfo object has identical pT_cuts as every other.
    wrong_eta_binning = (binning_type in 'eta') and len(set([(hist.pT_range[0], hist.pT_range[1]) for hist in entire_HistInfo_list])) != 1
    # Do same thing for binning in pT.
    wrong_pT_binning = (binning_type in 'pT') and len(set([(hist.eta_range[0], hist.eta_range[1]) for hist in entire_HistInfo_list])) != 1
    if (wrong_eta_binning or wrong_pT_binning):
        err_msg = f"\n\nBinning type ({binning_type}) specified, "
        err_msg += f"but not all graphs share same {binning_type}_range. Stopping now."
        raise RuntimeError(err_msg)
        
    text_size_legend = 6
    text_size_ax_labels = 12
    text_size_tick_labels = 8
    text_size_title = 12

    fig, ax = plt.subplots()

    al=1  # alpha=0 is transparent
    elw=1  # error bar line width
    ecolor='k'
    ms=4  # marker size
    mec='k'  # marker edge color
    cs=1  # cap size
    mew=0.7  # marker edge width

    d0_charge_latex_str = r'd_{0}^{%s} * \mathrm{charge}(\mu)' % example_graph.bspv
    delta_pT_latex_str = r"\Delta p_{T}/(p_{T}^{GEN})^{2}"
    
    # Having trouble aligning axis label with edges of axis...
    ax.set_xlabel(r'$%s$ [cm]' % d0_charge_latex_str, )
    ax.set_ylabel(r'mean of $%s$ [GeV]$^{-1}$' % delta_pT_latex_str, )

    title_str  = r'$%s$ vs. $%s$' % (delta_pT_latex_str, d0_charge_latex_str)
    title_str += f' for {example_graph.year} {example_graph.sample}' + "\n"
    # If binning in eta, then each graph should the same pT cuts. 
    if binning_type in 'eta':
        title_str += f'{example_graph.pT_range[0]}<' + r'$p_T$' + f'<{example_graph.pT_range[1]} GeV'
    elif binning_type in 'pT':
        title_str += f'{example_graph.eta_range[0]}<' + r'$\left| \eta \right|$' + f'<{example_graph.eta_range[1]}'
    ax.set_title(title_str)

    x_min = min([gr.d0_bin_window[0] for gr in example_graph_ls])
    x_max = max([gr.d0_bin_window[1] for gr in example_graph_ls])
    bin_edge_left = example_graph.d0_bin_window[0]
    bin_edge_right = example_graph.d0_bin_window[1]
    x_bin_width = abs(bin_edge_right - bin_edge_left)

    # ax.set_xlim([x_min*x_axis_scale, x_max*x_axis_scale])
    ax.set_xlim([-0.012, 0.012])
    ax.set_ylim([-0.005, 0.005])
    ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0), useMathText=True)
    
    # # We change the fontsize of minor ticks label (like the numbers on the axis)
    ax.tick_params(axis='both', which='major', direction='in', length=9,   top=True, right=True)
    ax.tick_params(axis='both', which='minor', direction='in', length=4.5, top=True, right=True)
    plt.minorticks_on()

    # Gridlines.
    ax.grid(which='major',color='k', ls=':')

    x_bin_arr = np.arange(x_min, x_max+0.5*x_bin_width, x_bin_width)  # Includes all bin edges: very first to very last.
    x_n_bins = calc_num_bins(x_min, x_max, x_bin_width)

    # Make sure points are plotted in middle of bin window.
    # Don't include right-most bin edge.
    x_bin_arr_shifted = x_bin_arr[0:-1] + 0.5*x_bin_width  

    low_x_err, high_x_err = calc_x_err_bins(x_bin_arr_shifted)

    # Calculate %diff for ratio plot.
    # acc_ZD_arr = np.array(acc_ZD, dtype=float)
    # acc_ALP_arr = np.array(acc_ALP, dtype=float)
    # acc_perc_diff_arr = (acc_ALP_arr - acc_ZD_arr) / acc_ZD_arr

    # Calculate errors on %diff.
    # acc_perc_diff_err_arr = getUncertOfFractionBinomial(acc_ALP_arr, acc_ZD_arr)

    # Plot the data.
    for count, gr_ls in enumerate(graph_tuple, 1):
        example_gr = gr_ls[0]
        y_vals = [gr.hist_mean for gr in gr_ls]
        y_vals_err = [gr.hist_mean_err for gr in gr_ls]
        if (verbose): 
            print(f'y_vals for graph_ls {count}:\n{y_vals}\n')
        if binning_type in 'pT':
            legend_text = f'{example_gr.pT_range[0]} <' + r' $p_T$ ' + f'< {example_gr.pT_range[1]}'
        elif binning_type in 'eta':
            legend_text = f'{example_gr.eta_range[0]} <' + r' $\left| \eta \right|$ ' + f'< {example_gr.eta_range[1]}'
        ax.errorbar(x_bin_arr_shifted, y_vals, xerr=[low_x_err, high_x_err], yerr=y_vals_err, fmt='s',
                label=legend_text,
                color=color_dict[count], elinewidth=elw, ms=ms, markeredgecolor=mec, capsize=cs, mew=mew, ecolor=ecolor)
    ax.legend(loc='upper right', framealpha=1, fontsize=text_size_legend)
    
    if (save_plots):
        plotname  = (f'deltapT_vs_d0{example_graph.bspv}timescharge_{example_graph.sample}{example_graph.year}_'
                     f'{binning_type}binning__')
        if binning_type in 'pT':
            plotname += f'{example_graph.eta_range[0]}_eta_{example_graph.eta_range[1]}'
        elif binning_type in 'eta':
            plotname += f'{example_graph.pT_range[0]}_pT_{example_graph.pT_range[1]}'
        plotname = make_str_title_friendly(plotname)
        
        fullpath = os.path.join(outpath, plotname)
        if not os.path.exists(outpath):
            os.makedirs(outpath)
        plt.savefig(fullpath + '.pdf')
        plt.savefig(fullpath + '.png')
        print("Plot saved at:\n", fullpath)
    
    plt.show()

# def make_kinem_comparison_plot(df, year, sample, 
# 
#     
#     pdf_name  = f"{year}_{sample}_"
#     title_str_pT_min = f"0{pT_min}" if pT_min < 10 else f"{pT_min}"  # For plot-ordering purposes.
#     # Example of Python magic: auto-concatenation of strings.
#     pdf_name += (f"{kinem_gen}_vs_{kinem_rec}__{bspv}"
#                  f"__{title_str_pT_min}_pT_{pT_max}"
#                  f"__{eta_min_abs}_abs_eta_{eta_max_abs}"
#                  f"__{d0_min:.3f}_to_{d0_max:.3f}_increm_{d0_bin_width:.3f}")
#                 #  f"__wrt_{wrt}")
#     pdf_name = make_str_title_friendly(pdf_name) + ".pdf"
#     
#     outfile = os.path.join(outfile_path, pdf_name)
#     if not os.path.exists(outfile_path):
#         os.makedirs(outfile_path)
#     if os.path.exists(outfile) and not (overwrite):
#         print(f"Skipping {outfile} since it already exists.\nTo write over the file then set overwrite=True.\n")
#         return
#         
#     status = f"Running over: {year} {sample} {bspv}, pT_range={pT_cut_ls}, eta_range={eta_cut_ls}, wrt {wrt}"
#     print(status)
# 
# 
#         for k in range(len(d0_bin_arr)-1):
#             this_d0_bin = d0_bin_arr[k]
#             next_d0_bin = d0_bin_arr[k+1]
#             if bspv in 'BS':
#                 mask_d0BSxcharge = (this_d0_bin < d0BSxcharge_ser) & (d0BSxcharge_ser < next_d0_bin)
#                 # mask_d0BS2xcharge = (this_d0_bin < d0BS2xcharge_ser) & (d0BS2xcharge_ser < next_d0_bin)
#             elif bspv in 'PV':
#                 mask_d0PVxcharge = (this_d0_bin < d0PVxcharge_ser) & (d0PVxcharge_ser < next_d0_bin)
#                 # mask_d0PV2xcharge = (this_d0_bin < d0PV2xcharge_ser) & (d0PV2xcharge_ser < next_d0_bin)
# 
# 
#             
#             # Example of Python magic. Auto-concatenation of strings.
#             cuts = f"{this_d0_bin:.3f}" + r"$<d_{0}^{%s}*\mathrm{charge}(\mu)<$"%bspv + f"{next_d0_bin:.3f}, " + "\n"
#             cuts += f"{pT_min}<" + r"$p_{T}$" + f"<{pT_max} GeV, "
#             cuts += f"{eta_min_abs}<" + r"$\left| \eta \right|$<" + f"{eta_max_abs}, "
#             cuts += f"{massZ_min} GeV<" + r"$m_{\mu\mu}$"
#                 
#  
#             # Store info of this plot in HistInfo objects. 
#             hist_gen = HistInfo()
#             hist_rec = HistInfo()            
#             hist_gen.year = year
#             hist_rec.year = year
#             hist_gen.sample = sample 
#             hist_rec.sample = sample
#             hist_gen.bspv = bspv 
#             hist_rec.bspv = bspv
#             hist_gen.d0_bin_window = [this_d0_bin, next_d0_bin] 
#             hist_rec.d0_bin_window = [this_d0_bin, next_d0_bin]
#             hist_gen.x_axis_bounds_list = x_range_ls[0:2] # 2-element list 
#             hist_rec.x_axis_bounds_list = x_range_ls[0:2] # 2-element list
#             hist_gen.eta_range = eta_cut_ls 
#             hist_rec.eta_range = eta_cut_ls
#             hist_gen.pT_range = pT_cut_ls
#             hist_rec.pT_range = pT_cut_ls
#             hist_gen.massZ_cut = massZ_min
#             hist_rec.massZ_cut = massZ_min
#             hist_gen.n_entries = n_entries_gen
#             hist_rec.n_entries = n_entries_rec
#             hist_gen.hist_mean = mean_gen
#             hist_rec.hist_mean = mean_rec
#             hist_gen.hist_mean_err = mean_err_gen
#             hist_rec.hist_mean_err = mean_err_rec
#             hist_gen.hist_stdev = stdev_gen
#             hist_rec.hist_stdev = stdev_rec  # Spread of the data.
#             hist_gen.hist_stdev_err = stdev_err_gen
#             hist_rec.hist_stdev_err = stdev_err_rec
#             
#             graph_info_gen_ls.append(hist_gen)
#             graph_info_rec_ls.append(hist_rec)
 
            
def make_kinem_subplot(lep, ax, data, x_limits, x_bins, x_label, y_label, y_max=-1, log_scale=False):
    """
    MAY BE DEPRECATED.
    Draw a kinematic distribution (e.g. eta1, gen_phi2, etc.) to a subplot in a figure.
    
    Parameters
    ----------
    lep : int
        Either `1` or `2`. Indicates which lepton you are referring to. 
    ax : axis object
        The axes to which plot will be drawn.
    data : array-like
        Data to be histogrammed. 
    x_limits : 2-element list
        A list of the x-min and x-max to be plotted.
    x_bins : array-like 
        Array of bin edges. Should be of length = len(data)+1. <-- I don't think that's true...
    x_label : str
        Label for x-axis.
    y_label : str
        Label for y-axis.
    label_legend : str
        Label for legend.
    y_max : float
        Max on y-axis. If y_max <= 0, then matplotlib will choose y_max.
    log_scale : bool
        If True, sets the y-axis to log scale. 
    """
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    
    ax.set_xlim(x_limits)
    ax.grid(False)
    
    mod_data, n_underflow, n_overflow = add_underoverflow_entries(data, x_limits[0], x_limits[1])
    
    if y_max > 0: ax.set_ylim([0,y_max])
    if (log_scale): ax.set_yscale('log')
        
    stats = get_stats_1Dhist(data)
    label_legend = make_stats_legend_for_1dhist(stats)
    bin_vals, bin_edges, _ = ax.hist(mod_data, bins=x_bins, label=label_legend, histtype='step', color='b')
    ax.legend(loc='upper right', fontsize=7)
    return
  
# def make_1D_dist(ax, data, x_limits, x_bins, x_label, y_label, title, y_max=-1, log_scale=False):
#    """
#    Draw a kinematic distribution (e.g. eta1, gen_phi2, etc.) to an axes object in a figure.
#    
#    Parameters
#    ----------
#    ax : axes object
#        The axes to which plot will be drawn.
#    data : array-like
#        Data to be histogrammed. 
#    x_limits : 2-element list
#        A list of the x-axis range to be plotted: 
#        x_limits=[x-min, x-max]
#    x_bins : array-like 
#        Array of bin edges. Should be of length = len(n_bins)+1.
#    x_label : str
#        Label for x-axis.
#    y_label : str
#        Label for y-axis.
#    title : str
#        Label for title.
#    y_max : float
#        Max on y-axis. If y_max <= 0, then matplotlib will choose y_max.
#    log_scale : bool
#        If True, sets the y-axis to log scale. 
#    """
#    textsize_legend = 9
#    textsize_axislabels = 12
#    textsize_title = 12
#            
#    ax.set_xlabel(x_label, fontsize=textsize_axislabels)
#    ax.set_ylabel(y_label, fontsize=textsize_axislabels)
#    ax.set_title(title, fontsize=textsize_title)
#    
#    ax.set_xlim(x_limits)
#    ax.grid(False)
#    
#    mod_data, n_underflow, n_overflow = add_underoverflow_entries(data, x_limits[0], x_limits[1])
#    
#    if y_max > 0: ax.set_ylim([0,y_max])
#    if (log_scale): ax.set_yscale('log')
#        
#    stats = get_stats_1Dhist(data)
#    label_legend = make_stats_legend_for_1dhist(stats)
#    bin_vals, bin_edges, _ = ax.hist(mod_data, bins=x_bins, label=label_legend, histtype='step', color='b')
#    ax.legend(loc='upper right', framealpha=0.9, fontsize=textsize_legend)
#    
#    return ax, stats
    
def print_header_message(msg):
    n = len(msg)
    octothorpes = (n+12)*'#'  # Who names their variable 'octothorpes'? Honestly?
    buff = 5*'#'
    print(octothorpes)
    print(buff, msg, buff)
    print(octothorpes)
    
def make_binning_array(lim_ls):
    """
    Turn a list of [bin_min, bin_max, bin_width] into an array of length = (bin_max - bin_min)/bin_width,
    whose first element is bin_min, last element is bin_max and has spacing bin_width.
    
    Returns
    -------
    bin_edges : array
        Contains equally-spaced bin edges, where:
            bin_edges[0] = left edge of first bin
            bin_edges[-1] = right edge of last bin
    binw : float
        The width of a bin. All are equal width.
    """
    ls_min = lim_ls[0]
    ls_max = lim_ls[1]
    binw = lim_ls[2]
    
    bin_edges = np.arange(ls_min, ls_max+0.5*binw, binw)
    
    return bin_edges, binw

def shift_binning_array(bin_arr):
    """
    Returns an array of the centers of the bins of bin_arr, excluding the very last bin.
    
    Example: 
        bin_arr = np.array([1, 2, 3, 4])
        bin_arr_shifted = np.array([1.5, 2.5, 3.5])
        
    FIXME: Somewhere in my DarkZ code is this function which can accommodate an asymmetrical bin_arr.
    """
    bin_width = bin_arr[1] - bin_arr[0]
    
    return bin_arr[:-1] + 0.5 * bin_width

def event_counter(start, stop, step, df, branch):
    """
    I THOUGHT I WAS CLEVER WHEN WRITING THIS FUNCTION, BUT THE MUCH MORE EFFICIENT WAY IS:
    
    vals, bin_edges = np.histogram(df[branch], bins=np.arange(start, stop * step/2, step))
    return vals
    
    ORIGINAL:
    Counts the number of events in each bin.
    This function provides a way to make sure that your histogram is binning the correct number of events per bin.
    
    Parameters
    ----------
    start : float 
        The left edge of the very first bin to begin counting.
    stop : float 
        The right edge of the very last bin to stop counting. 
    step : float
        The bin width. 
    df : pandas.DataFrame
        The DataFrame that holds branch.
    branch : str
        The DataFrame column name that holds the data to be binned.
    """
    # vals = np.arange(start, stop, step)
    # for v in vals:
    #     print(f"n_evts with {v:.1f} < eta < {v+step : .1f}:", sum((v < df[branch]) & (df[branch] < v+step) ) )
    vals, bin_edges = np.histogram(df[branch], bins=np.arange(start, stop * step/2, step))
    return vals