"""
# Purpose: 
#   This code produces 1 PDF where on each page a graph is displayed. 
#   A single graph has only 1 eta bin. 
#   Each line on the graph is a different pT bin. 
#   Each point on the line has coordinates: 
#     (x,y) = (avg(q*d0) of that cube, best-fit gaus mean of dpT/pT dist)
#   This code also creates a dict of linear best-fit params
#   and a readable '.txt' file of the same dict.
# SYNTAX: python script.py
# NOTES: 
#   The eta,pT,qd0 etc. data come from a '.pkl' file, specifed by the User.
#   User can specify to scale y-axis by 1/avg(pT) of each cube
#   by setting: scale_by_1divpT = True
# AUTHOR: Jake Rosenzweig
# UPDATED: 2020-06-26
"""
import os
import pickle
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.backends.backend_pdf import PdfPages

# Local imports.
from d0_Utils.d0_cls import KinBin3DOrganizer, GraphLineKinBin3D
from Utils_Python.Utils_Files import make_dirs, make_str_title_friendly, check_overwrite
#---------------------------#
#----- User Parameters -----#
#---------------------------#
fullpath_kinbin_ls_pkl = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pkl_and_csv/kinbins/20200624_fullstats_MC2018JpsiDY_12regwith3000perreg__0p0_eta_2p4__5p0_pT_1000p0_GeV_unbinned_kinbin_ls.pkl"
outdir_plots = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/graphs_dpToverpT/"
outdir_param_dict = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/KinBin_Info/"
pdfname_base = "MC2018_dpToverpTvsqd0_graph_fullstats_paramsrecorded"

plot_style_sheet = "/ufrc/avery/rosedj1/HiggsMassMeasurement/Utils_Python/Plot_Styles_Matplotlib/cmsstyle_plot.mplstyle"

overwrite = False
verbose = True
scale_by_1divpT = False

def prep_area(outdir_plots, outpath, overwrite):
    make_dirs(outdir_plots)
    check_overwrite(outpath, overwrite=overwrite)

def add_params_to_dict(d, param_ls, eta_range, pT_range):
    """
    Take list of best-fit parameters [intercept, slope] and add them to dict (d).
    eta_range and pT_range determine the key.
    [intercept, slope] determine the val.

    Key : val is of the form
        "0.0eta0.2_5.0pT7.0" : {"intercept":1.475E-04, "slope":6.510E-01} 
    """
    eta_min = eta_range[0]
    eta_max = eta_range[1]
    pT_min = pT_range[0]
    pT_max = pT_range[1]
    
    key = f"{eta_min}eta{eta_max}_{pT_min}pT{pT_max}"
    assert key not in d
    d[key] = {"intercept":param_ls[0], "slope":param_ls[1]}

def make_1plot_allgraphs(kb_org, best_params_dict, const_range=[], const_bin="", 
                            ax=None, x_lim=None, y_lim=None, 
                            scale_by_1divpT=False, verbose=True):
    ls_of_kb_ls = kb_org.find_all_KinBin_ls_const_range(const_range=const_range, const_bin=const_bin)
    
    if ax is None:
        fig, ax = plt.subplots()

    for count, kb_ls in enumerate(ls_of_kb_ls, 1):
        # Each kb in kb_ls has the same eta_range and pT_range.
        qd0_avg_ls, dpToverpT_bestfitmean_ls, dpToverpT_bestfitmean_err_ls = kb_org.get_plotting_vals_from_KinBin_ls(kb_ls)

        eta_range = kb_ls[0].eta_range
        pT_range = kb_ls[0].pT_range
        line = GraphLineKinBin3D(x_vals=qd0_avg_ls, y_vals=dpToverpT_bestfitmean_ls, 
                                 x_err_vals=None, y_err_vals=dpToverpT_bestfitmean_err_ls,
                                 eta_range=eta_range, pT_range=pT_range)
        # Doing line.draw_graph will also find and set the best-fit params.
        line.draw_graph(x_label="", y_label="", title="", 
                        kbin_ls=kb_ls, ax=ax, x_lim=x_lim, y_lim=y_lim, count=count,
                        scale_by_1divpT=scale_by_1divpT, verbose=verbose, const_bin=const_bin)
        add_params_to_dict(best_params_dict, line.best_fit_params, eta_range, pT_range)

def make_all_graphs_all_ranges(kb_org, best_params_dict, const_bin="", scale_by_1divpT=False, 
                               x_lim=None, y_lim=None, pdf=None, verbose=verbose):
    if const_bin == "eta":
        ls_of_lists = kb_org.eta_range_ls
    elif const_bin == "pT":
        ls_of_lists = kb_org.pT_range_ls
    else:
        raise ValueError(f"[ERROR] Unknown const_bin str ({const_bin})")
        
    print(f"#----- Making {len(ls_of_lists)} plots -----#")
    print(f" - Constant binning:     {const_bin}")
    print(f" - scaling by 1/avg(pT): {scale_by_1divpT}")
    
    for const_range in ls_of_lists:
        make_1plot_allgraphs(kb_org, best_params_dict=best_params_dict, 
                             const_range=const_range, const_bin=const_bin, 
                             ax=None, x_lim=x_lim, y_lim=y_lim,
                             scale_by_1divpT=scale_by_1divpT, verbose=verbose)
        
        if pdf is not None:
            pdf.savefig()
            plt.close("all")

def get_kinbin_ls_from_pkl(inpkl):
    with open(inpkl, "rb") as f:
        return pickle.load(f)

def set_plot_style():
    plt.style.use(plot_style_sheet)
    rc('lines', linewidth=0.5, markersize=0.1,)
    rc('axes', labelsize=7, titlesize=10)
    rc('legend', fontsize=2)
    # rc('errorbar', capsize=0.0)

def pickle_param_dict(outpath_pkl, param_dict):
    with open(outpath_pkl, 'wb') as output:
        pickle.dump(param_dict, output, protocol=2)
    print(f"[INFO] Best-fit param dict written to pickle:\n{outpath_pkl}")

def make_outdict_names(outdir_param_dict, pdfname_base):
        outpath_file_woextension = os.path.join(outdir_param_dict, f"{pdfname_base}_bestfitparamdict")
        outpath_pkl = outpath_file_woextension + ".pkl"
        outpath_txt = outpath_file_woextension + ".txt"
        return outpath_pkl, outpath_txt
    
def write_param_dict(outpath_txt, outpath_plots, best_params_dict):
    with open(outpath_txt, "a") as f:
        f.write(f"From {outpath_plots}\n\n")
        f.write("qd0_params_dict = {\n")
        for key,val in best_params_dict.items():
            f.write(f"    '{key}' : {val},\n")
        f.write("}")
        print(f"Best-fit param dictionary put into '.txt' file:\n{outpath_txt}")

if __name__ == "__main__":
    best_params_dict = {}
    outpath_plots = os.path.join(outdir_plots, pdfname_base + ".pdf")
    prep_area(outdir_plots, outpath_plots, overwrite)
    set_plot_style()

    # Get data.
    kinbin_ls = get_kinbin_ls_from_pkl(fullpath_kinbin_ls_pkl)
    kb_org = KinBin3DOrganizer(kinbin_ls)

    with PdfPages(outpath_plots) as pdf:
        # Test just 1 graph. 
        # make_1plot_allgraphs(kb_org, const_range=[0.6, 0.8], const_bin="eta", 
        #                     ax=None, x_lim=[-0.0065, 0.0065], y_lim=[-0.08, 0.08], 
        #                     scale_by_1divpT=False, verbose=verbose)
        # pdf.savefig()
        # plt.close("all")

        # Run over all eta ranges.
        if (scale_by_1divpT):
            make_all_graphs_all_ranges(kb_org, best_params_dict=best_params_dict, 
                                        const_bin="eta", scale_by_1divpT=True,
                                        x_lim=[-0.0065, 0.0065], y_lim=[-0.003, 0.003], 
                                        pdf=pdf, verbose=False)
        else:
            make_all_graphs_all_ranges(kb_org, best_params_dict=best_params_dict, 
                                        const_bin="eta", scale_by_1divpT=False,
                                      x_lim=[-0.0065, 0.0065], y_lim=[-0.02, 0.02], 
                                      pdf=pdf, verbose=False)
                        #    x_lim=[-0.0065, 0.0065], y_lim=[-0.08, 0.08], pdf=pdf, verbose=False)
    print(f"PDF made: {outpath_plots}")

    outpath_pkl, outpath_txt = make_outdict_names(outdir_param_dict, pdfname_base)
    pickle_param_dict(outpath_pkl, best_params_dict)
    write_param_dict(outpath_txt, outpath_plots, best_params_dict)

    #     all_pT_ranges = getattr(kb_org, "pT_range_ls")
        
    #     for pT_range in all_pT_ranges:

    # fig, ax = plt.subplots()
        
        
    # def draw_both_plots_in_const_range(kb_org, const_range=[], const_bin="", ax=None, scale_by_1divpT=False, verbose=False):
    #     fig, ax = plt.subplots()
    #     make_1plot_allgraphs(kb_org, const_range=[0.0, 0.2], const_bin="eta", ax=ax, scale_by_1divpT=False, verbose=verbose)
    # #     plt.tight_layout()
        

    #     # Scale by 1/avg(pT).
    #     fig, ax = plt.subplots()
    #     make_1plot_allgraphs(kb_org, const_range=[0.0, 0.2], const_bin="eta", ax=ax, scale_by_1divpT=True, verbose=verbose)
    #     pdf.savefig()

    # draw_both_plots_in_const_range()

    # plt.tight_layout()
    #     pdf.savefig()
    # plt.close("all")