"""
# PURPOSE: 
#  This code reads in values from a pickled dictionary, calculates the 
#  % improvement per [eta, pT] bin, and puts the improvement into 
#  a 2D histogram. Finally it draws the histogram to a PDF.
# NOTES: 
#  The values contained in the pickled dict are the best-fit sigma values
#  from iterated Gaussian fits of kinematic distributions, like qd0 or dpT/pT.
#  The keys in the dict look like: 2.3eta2.4_150.0pT200.0
#  - Be sure to check all the variables in --- User Parameters ---.
#  - If you get blank cells, try modifying your color_lim range. 
# SYNTAX: python this_script.py
# AUTHOR: Jake Rosenzweig
# EDITED: 2021-03-03
"""
import pickle
import os
import numpy as np
import ROOT as rt
# Package imports.
from Utils_Python.Utils_Files import check_overwrite, open_pkl, make_dirs
from Utils_Python.Utils_StatsAndFits import get_status_redchi2_fit
from Utils_ROOT.ROOT_classes import make_TH2F
#--- User Parameters ---#
# inpath_pkl = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/qd0/MC{year}JpsiDY_2D_plot_qd0dist_gausiterfitsigmas_4unbinnedfits_0p0eta2p4_5p0pT1000p0GeV.pkl"
filename = "test15.pdf"
inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/RC_vs_NoRC_itergaussfits_fullstats_pT75then200GeV_extendedxaxis_5iters.pkl"
# outpath_pdf = f"/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/plots/2D_tables_etavspT/MC{year}JpsiDY_2Dplot_dpToverpTimprovement_gausiterfitsigmas_6unbinnedfits_0p0eta2p4_5p0pT1000p0GeV_final.pdf"
outpath_pdf = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/plots/tables2D/test/{filename}"
overwrite = 0
show_plots = 0
make_pdf = 1

min_chi2 = 0.01
max_chi2 = 199

z_label_size = 0.02
text_size = 0.85
n_contour = 300
rt.gStyle.SetPalette(rt.kTemperatureMap) 

unit_factor = 1000.0  # GeV -> MeV.
auto_detect_binedges = True

# eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
# pT_binedge_ls = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]

# Number format: e.g. "6.2f" could be 12.34
# <padding>.<decimals><type>
cell_text_format = ".3g"
#--- Local Functions ---#
def fill_TH2F(h2, x_vals, y_vals, z_vals):
    """Fill a 2-D hist with (x,y) entries (coordinates).
    
    FIXME: use_weight_x and use_weight_y are deprecated.
    Parameters
    ----------
    h2 : ROOT.TH2F
        The 2-D hist to be filled with len(x_coord_ls) entries.
    x_coord_ls : list or array-like
        The x-values of each entry. The first element pairs with the first in y_coord_ls.
    y_coord_ls : list or array-like
        The y-values of each entry. The first element pairs with the first in x_coord_ls.
    """
    assert len(x_vals) == len(y_vals)
    for x, y, z in zip(x_vals, y_vals, z_vals):
        h2.Fill(x, y, z)

def get_binedges_from_kb2d_dct(kb2d_dct):
    """Return a 2-tuple of lists: (eta_bin_edges, pT_bin_edges) which are
    automatically determined from kb2d_dct.
    """
    pT_edges_all = []
    eta_edges_all = [] 
    for kb2d in kb2d_dct.values():
        pT_edges_all.append(kb2d.pT_range[0])
        pT_edges_all.append(kb2d.pT_range[1])
        eta_edges_all.append(kb2d.eta_range[0])
        eta_edges_all.append(kb2d.eta_range[1])
    pT_binedge_ls = sorted(set(pT_edges_all))
    eta_binedge_ls = sorted(set(eta_edges_all))
    return (eta_binedge_ls, pT_binedge_ls)

def get_standarderrorofmean(arr):
    """Return the standard error of the mean."""
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr)
    return np.std(arr) / np.sqrt(len(arr))

def get_mean_and_SEOM(arr):
    """Return a 2-tuple: (mean, standard error of mean)."""
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr)
    mean = np.mean(arr)
    seom = get_standarderrorofmean(arr)
    return (mean, seom)

def get_mean_and_SEOM_fromiterfits(kb2d, which_dct, min_chi2, max_chi2):
    """Return a 2-tuple: (mean, standard error of mean) from iterated Gaussian fits.
    
    Parameters
    ----------
    kb2d : KinBin2D
    which_dct : str
        A dict key to unlock the appropriate fit stats dict.
        Choices:
        - dpT_RC_reco_fit_dct
        - dpTOverpT_RC_reco_fit_dct
        - dpT_RC_gen_fit_dct
        - dpTOverpT_RC_gen_fit_dct
        - dpT_reco_gen_fit_dct
        - dpTOverpT_reco_gen_fit_dct
    min_chi2 : float
        The min chi^2 allowed to be considered a good fit.
    max_chi2 : float
        The max chi^2 allowed to be considered a good fit.
    """
    mean_ls = getattr(kb2d, which_dct)['mean_ls']
    mean_err_ls = getattr(kb2d, which_dct)['mean_err_ls']
    chi2_ls = getattr(kb2d, which_dct)['chi2_ls']
    # Best-fit vals are at the end of lists:
    good_fit = False
    # Start at end of chi^2 list and see if it is within allowed limits.
    for chi2 in chi2_ls[::-1]:
        if get_status_redchi2_fit(chi2, min_chi2, max_chi2):
            good_fit = True
            best_index = chi2_ls.index(chi2)
            break
        else:
            print(f"bad chi^2: {chi2}")
    assert good_fit, "Couldn't find reduced chi^2 within allowed range."

    mean = mean_ls[best_index]
    seom = mean_err_ls[best_index]
    # mean = mean_ls[-1]
    # seom = mean_err_ls[-1]
    return (mean, seom)

def make_hist_and_errhist(internal_name, title=None, 
            n_binsx=5, x_label="", x_units=None, x_min=0, x_max=10, z_min=None,
            n_binsy=5, y_label="", y_units=None, y_min=0, y_max=10, z_max=None,
            z_label_size=None,
            n_contour=100):
    """Return a 2-tuple: (TH2F filled with (x,y,z), TH2F filled with (x,y,z_err)).
    
    NOTE: Make sure """
    h = make_TH2F(internal_name=internal_name, title=title, 
                n_binsx=pT_binedge_ls, x_label=x_label, x_units=x_units,
                n_binsy=eta_binedge_ls, y_label=y_label, y_units=y_units, y_min=y_min, y_max=y_max,
                n_contour=n_contour)
    h_err = h.Clone()
    h_err.SetName(h.GetName().replace("mean", "sterrofmean"))
    return (h, h_err)

def set_bin_vals_and_errs(h2, h2_err):
    """Set the values stored in h2_err as the errors in h2."""
    assert h2.GetNbinsX() == h2_err.GetNbinsX()
    assert h2.GetNbinsY() == h2_err.GetNbinsY()
    for binx in range(1, h2.GetNbinsX() + 1):
        for biny in range(1, h2.GetNbinsY() + 1):
            # Similar cells (e.g. (2,3)) in both hists correspond to each other.
            err = h2_err.GetBinContent(binx, biny)
            h2.SetBinError(binx, biny, err)

if __name__ == "__main__":
    check_overwrite(outpath_pdf, overwrite=overwrite)
    make_dirs(os.path.dirname(outpath_pdf))
    if not show_plots:
        rt.gROOT.SetBatch(True)
    dct = open_pkl(inpath_pkl)

    if auto_detect_binedges:
        eta_binedge_ls, pT_binedge_ls = get_binedges_from_kb2d_dct(dct)
    x_label = "p_{T}^{reco}"
    rt.gStyle.SetOptStat(0)

    # Make hist of mean and standard err of mean for pTs, dpTs, and dpT/pTs.
    #------------------------------------#
    #--- Make 2D hists for plain data ---#
    #------------------------------------# 
    h2_mean_pT_RC, h2_sterrofmean_pT_RC = make_hist_and_errhist(
                                              internal_name="h2_mean_pT_RC", title="<p_{T}^{RC}> (GeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=None, z_max=None, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_pT_reco, h2_sterrofmean_pT_reco = make_hist_and_errhist(
                                              internal_name="h2_mean_pT_reco", title="<p_{T}^{reco}> (GeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=0, z_max=200, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_pT_gen, h2_sterrofmean_pT_gen = make_hist_and_errhist(
                                              internal_name="h2_mean_pT_gen", title="<p_{T}^{gen}> (GeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=0, z_max=200, z_label_size=z_label_size,
                                              n_contour=n_contour)
    # pT differences.
    h2_mean_pTRCminuspTreco, h2_sterrofmean_pTRCminuspTreco = make_hist_and_errhist(
                                              internal_name="h2_mean_pTRCminuspTreco", title="<p_{T}^{RC} - p_{T}^{reco}> (MeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-6, z_max=600, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_pTRCminuspTgen, h2_sterrofmean_pTRCminuspTgen = make_hist_and_errhist(
                                              internal_name="h2_mean_pTRCminuspTgen", title="<p_{T}^{RC} - p_{T}^{gen}> (MeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-6, z_max=600, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_pTrecominuspTgen, h2_sterrofmean_pTrecominuspTgen = make_hist_and_errhist(
                                              internal_name="h2_mean_pTrecominuspTgen", title="<p_{T}^{reco} - p_{T}^{gen}> (MeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-6, z_max=600, z_label_size=z_label_size,
                                              n_contour=n_contour)
    # Relative pT differences.
    h2_mean_relpTRCminuspTreco, h2_sterrofmean_relpTRCminuspTreco = make_hist_and_errhist(
                                              internal_name="h2_mean_relpTRCminuspTreco", title="<(p_{T}^{RC} - p_{T}^{reco})/p_{T}^{reco}> (%)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-2, z_max=20, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_relpTRCminuspTgen, h2_sterrofmean_relpTRCminuspTgen = make_hist_and_errhist(
                                              internal_name="h2_mean_relpTRCminuspTgen", title="<(p_{T}^{RC} - p_{T}^{gen})/p_{T}^{gen}> (%)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-2, z_max=20, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_relpTrecominuspTgen, h2_sterrofmean_relpTrecominuspTgen = make_hist_and_errhist(
                                              internal_name="h2_mean_relpTrecominuspTgen", title="<(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}> (%)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-2, z_max=20, z_label_size=z_label_size,
                                              n_contour=n_contour)
    
    #--------------------------------------#
    #--- Make 2D hists for iterfit data ---#
    #--------------------------------------# 
    h2_mean_pTRCminuspTreco_iterfit, h2_sterrofmean_pTRCminuspTreco_iterfit = make_hist_and_errhist(
                                              internal_name="h2_mean_pTRCminuspTreco_iterfit", title="Best-fit Gaussian #mu#left(p_{T}^{RC} - p_{T}^{reco}#right) (MeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-6, z_max=600, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_pTRCminuspTgen_iterfit, h2_sterrofmean_pTRCminuspTgen_iterfit = make_hist_and_errhist(
                                              internal_name="h2_mean_pTRCminuspTgen_iterfit", title="Best-fit Gaussian #mu#left(p_{T}^{RC} - p_{T}^{gen}#right) (MeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-6, z_max=600, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_pTrecominuspTgen_iterfit, h2_sterrofmean_pTrecominuspTgen_iterfit = make_hist_and_errhist(
                                              internal_name="h2_mean_pTrecominuspTgen_iterfit", title="Best-fit Gaussian #mu#left(p_{T}^{reco} - p_{T}^{gen}#right) (MeV)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-6, z_max=600, z_label_size=z_label_size,
                                              n_contour=n_contour)
    # Relative pT differences.
    h2_mean_relpTRCminuspTreco_iterfit, h2_sterrofmean_relpTRCminuspTreco_iterfit = make_hist_and_errhist(
                                              internal_name="h2_mean_relpTRCminuspTreco_iterfit", title="Best-fit Gaussian #mu#left((p_{T}^{RC} - p_{T}^{reco})/p_{T}^{reco}#right) (%)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-0.5, z_max=0.5, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_relpTRCminuspTgen_iterfit, h2_sterrofmean_relpTRCminuspTgen_iterfit = make_hist_and_errhist(
                                              internal_name="h2_mean_relpTRCminuspTgen_iterfit", title="Best-fit Gaussian #mu#left((p_{T}^{RC} - p_{T}^{gen})/p_{T}^{gen}#right) (%)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-2, z_max=20, z_label_size=z_label_size,
                                              n_contour=n_contour)
    h2_mean_relpTrecominuspTgen_iterfit, h2_sterrofmean_relpTrecominuspTgen_iterfit = make_hist_and_errhist(
                                              internal_name="h2_mean_relpTrecominuspTgen_iterfit", title="Best-fit Gaussian #mu#left((p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}#right) (%)",
                                              n_binsx=pT_binedge_ls, x_label=x_label, x_units="GeV",
                                              n_binsy=eta_binedge_ls, y_label="#left|#eta#right|", y_units=None,
                                              z_min=-2, z_max=20, z_label_size=z_label_size,
                                              n_contour=n_contour)
    
    for kb2d in dct.values():
        # Collect iterated Gaussian fit data.
        avg_bin_val_eta = np.mean(kb2d.eta_range)
        avg_bin_val_pT = np.mean(kb2d.pT_range)

        mean_dpT_RCvsreco_iterfit, sterrofmean_dpT_RCvsreco_iterfit = get_mean_and_SEOM_fromiterfits(kb2d, "dpT_RC_reco_fit_dct", min_chi2, max_chi2)
        mean_dpT_RCvsgen_iterfit, sterrofmean_dpT_RCvsgen_iterfit = get_mean_and_SEOM_fromiterfits(kb2d, "dpT_RC_gen_fit_dct", min_chi2, max_chi2)
        mean_dpT_recovsgen_iterfit, sterrofmean_dpT_recovsgen_iterfit = get_mean_and_SEOM_fromiterfits(kb2d, "dpT_reco_gen_fit_dct", min_chi2, max_chi2)
        mean_reldpT_RCvsreco_iterfit, sterrofmean_reldpT_RCvsreco_iterfit = get_mean_and_SEOM_fromiterfits(kb2d, "dpTOverpT_RC_reco_fit_dct", min_chi2, max_chi2)
        mean_reldpT_RCvsgen_iterfit, sterrofmean_reldpT_RCvsgen_iterfit = get_mean_and_SEOM_fromiterfits(kb2d, "dpTOverpT_RC_gen_fit_dct", min_chi2, max_chi2)
        mean_reldpT_recovsgen_iterfit, sterrofmean_reldpT_recovsgen_iterfit = get_mean_and_SEOM_fromiterfits(kb2d, "dpTOverpT_reco_gen_fit_dct", min_chi2, max_chi2)

        h2_mean_pTRCminuspTreco_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_dpT_RCvsreco_iterfit)  # Does not need unit_factor.
        h2_sterrofmean_pTRCminuspTreco_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_dpT_RCvsreco_iterfit)

        h2_mean_pTRCminuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_dpT_RCvsgen_iterfit * unit_factor)
        h2_sterrofmean_pTRCminuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_dpT_RCvsgen_iterfit * unit_factor)

        h2_mean_pTrecominuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_dpT_recovsgen_iterfit * unit_factor)
        h2_sterrofmean_pTrecominuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_dpT_recovsgen_iterfit * unit_factor)

        h2_mean_relpTRCminuspTreco_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_reldpT_RCvsreco_iterfit)
        h2_sterrofmean_relpTRCminuspTreco_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_reldpT_RCvsreco_iterfit)

        h2_mean_relpTRCminuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_reldpT_RCvsgen_iterfit)
        h2_sterrofmean_relpTRCminuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_reldpT_RCvsgen_iterfit)

        h2_mean_relpTrecominuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_reldpT_recovsgen_iterfit)
        h2_sterrofmean_relpTrecominuspTgen_iterfit.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_reldpT_recovsgen_iterfit)
    # Collect normal data.
    for kb2d in dct.values():
        avg_bin_val_eta = np.mean(kb2d.eta_range)
        avg_bin_val_pT = np.mean(kb2d.pT_range)
        # pT_RC_arr = np.array([])
        # pT_reco_arr = np.array([])
        # pT_gen_arr = np.array([])
        pT_RC_arr = np.array([mu.pT_RC for mu in kb2d.muon_ls])
        pT_reco_arr = np.array([mu.pT for mu in kb2d.muon_ls])
        pT_gen_arr = np.array([mu.gen_pT for mu in kb2d.muon_ls])
        # pT differences:
        dpT_RCvsreco_arr = (pT_RC_arr - pT_reco_arr) * unit_factor  # MeV.
        dpT_RCvsgen_arr = (pT_RC_arr - pT_gen_arr) * unit_factor  # MeV.
        dpT_recovsgen_arr = (pT_reco_arr - pT_gen_arr) * unit_factor  # MeV.
        # Relative pT differences:
        reldpT_RCvsreco_arr = (dpT_RCvsreco_arr / (pT_reco_arr * unit_factor)) * 100.0  # %
        reldpT_RCvsgen_arr = (dpT_RCvsgen_arr / (pT_gen_arr * unit_factor)) * 100.0  # %
        reldpT_recovsgen_arr = (dpT_recovsgen_arr / (pT_gen_arr * unit_factor)) * 100.0  # %
        # Calculate  statistics.
        mean_pT_RC, sterrofmean_pT_RC = get_mean_and_SEOM(pT_RC_arr)
        mean_pT_reco, sterrofmean_pT_reco = get_mean_and_SEOM(pT_reco_arr)
        mean_pT_gen, sterrofmean_pT_gen = get_mean_and_SEOM(pT_gen_arr)
        mean_dpT_RCvsreco, sterrofmean_dpT_RCvsreco = get_mean_and_SEOM(dpT_RCvsreco_arr)
        mean_dpT_RCvsgen, sterrofmean_dpT_RCvsgen = get_mean_and_SEOM(dpT_RCvsgen_arr)
        mean_dpT_recovsgen, sterrofmean_dpT_recovsgen = get_mean_and_SEOM(dpT_recovsgen_arr)
        mean_reldpT_RCvsreco, sterrofmean_reldpT_RCvsreco = get_mean_and_SEOM(reldpT_RCvsreco_arr)
        mean_reldpT_RCvsgen, sterrofmean_reldpT_RCvsgen = get_mean_and_SEOM(reldpT_RCvsgen_arr)
        mean_reldpT_recovsgen, sterrofmean_reldpT_recovsgen = get_mean_and_SEOM(reldpT_recovsgen_arr)
        # Fill all hists.
        h2_mean_pT_RC.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_pT_RC)
        h2_sterrofmean_pT_RC.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_pT_RC)

        h2_mean_pT_reco.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_pT_reco)
        h2_sterrofmean_pT_reco.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_pT_reco)

        h2_mean_pT_gen.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_pT_gen)
        h2_sterrofmean_pT_gen.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_pT_gen)

        h2_mean_pTRCminuspTreco.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_dpT_RCvsreco)
        h2_sterrofmean_pTRCminuspTreco.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_dpT_RCvsreco)

        h2_mean_pTRCminuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_dpT_RCvsgen)
        h2_sterrofmean_pTRCminuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_dpT_RCvsgen)

        h2_mean_pTrecominuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_dpT_recovsgen)
        h2_sterrofmean_pTrecominuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_dpT_recovsgen)

        h2_mean_relpTRCminuspTreco.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_reldpT_RCvsreco)
        h2_sterrofmean_relpTRCminuspTreco.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_reldpT_RCvsreco)

        h2_mean_relpTRCminuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_reldpT_RCvsgen)
        h2_sterrofmean_relpTRCminuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_reldpT_RCvsgen)

        h2_mean_relpTrecominuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, mean_reldpT_recovsgen)
        h2_sterrofmean_relpTrecominuspTgen.Fill(avg_bin_val_pT, avg_bin_val_eta, sterrofmean_reldpT_recovsgen)
    # End loop over kb2ds.
    hist_tup = (h2_mean_pT_RC, h2_mean_pT_reco, h2_mean_pT_gen,
                h2_mean_pTRCminuspTreco, h2_mean_pTRCminuspTgen, h2_mean_pTrecominuspTgen,
                h2_mean_relpTRCminuspTreco, h2_mean_relpTRCminuspTgen, h2_mean_relpTrecominuspTgen)
    hist_err_tup = (h2_sterrofmean_pT_RC, h2_sterrofmean_pT_reco, h2_sterrofmean_pT_gen,
                    h2_sterrofmean_pTRCminuspTreco, h2_sterrofmean_pTRCminuspTgen, h2_sterrofmean_pTrecominuspTgen,
                    h2_sterrofmean_relpTRCminuspTreco, h2_sterrofmean_relpTRCminuspTgen, h2_sterrofmean_relpTrecominuspTgen)
    hist_iterfit_tup = (h2_mean_pTRCminuspTreco_iterfit, h2_mean_pTRCminuspTgen_iterfit, h2_mean_pTrecominuspTgen_iterfit,
                h2_mean_relpTRCminuspTreco_iterfit, h2_mean_relpTRCminuspTgen_iterfit, h2_mean_relpTrecominuspTgen_iterfit)
    hist_err_iterfit_tup = (h2_sterrofmean_pTRCminuspTreco_iterfit, h2_sterrofmean_pTRCminuspTgen_iterfit, h2_sterrofmean_pTrecominuspTgen_iterfit,
                    h2_sterrofmean_relpTRCminuspTreco_iterfit, h2_sterrofmean_relpTRCminuspTgen_iterfit, h2_sterrofmean_relpTrecominuspTgen_iterfit)
    
    z_max_pT = 105
    h2_mean_pT_RC.GetZaxis().SetRangeUser(0, z_max_pT)
    h2_mean_pT_reco.GetZaxis().SetRangeUser(0, z_max_pT)
    h2_mean_pT_gen.GetZaxis().SetRangeUser(0, z_max_pT)

    z_max_dpT = 600
    h2_mean_pTRCminuspTreco.GetZaxis().SetRangeUser(-200, 200)
    h2_mean_pTRCminuspTgen.GetZaxis().SetRangeUser(-1 * z_max_dpT, z_max_dpT)
    h2_mean_pTrecominuspTgen.GetZaxis().SetRangeUser(-1 * z_max_dpT, z_max_dpT)
    h2_mean_pTRCminuspTreco_iterfit.GetZaxis().SetRangeUser(-200, 200)
    h2_mean_pTRCminuspTgen_iterfit.GetZaxis().SetRangeUser(-1 * z_max_dpT, z_max_dpT)
    h2_mean_pTrecominuspTgen_iterfit.GetZaxis().SetRangeUser(-1 * z_max_dpT, z_max_dpT)

    z_max_reldpT = 2.1
    h2_mean_relpTRCminuspTreco.GetZaxis().SetRangeUser(-0.5, 0.5)
    h2_mean_relpTRCminuspTgen.GetZaxis().SetRangeUser(-1 * z_max_reldpT, z_max_reldpT)
    h2_mean_relpTrecominuspTgen.GetZaxis().SetRangeUser(-1 * z_max_reldpT, z_max_reldpT)
    h2_mean_relpTRCminuspTreco_iterfit.GetZaxis().SetRangeUser(-0.5, 0.5)
    h2_mean_relpTRCminuspTgen_iterfit.GetZaxis().SetRangeUser(-1 * z_max_reldpT, z_max_reldpT)
    h2_mean_relpTrecominuspTgen_iterfit.GetZaxis().SetRangeUser(-1 * z_max_reldpT, z_max_reldpT)
    if make_pdf:
        c = rt.TCanvas()
        c.Print(outpath_pdf + "[")
        c.SetLogx(True)
        # Consolidate for loops into a single one.
        for h, h_err in zip(hist_tup, hist_err_tup):
            set_bin_vals_and_errs(h, h_err)
            h.SetContour(n_contour)
            h.SetMarkerSize(text_size)
            # h.SetBarOffset(0.2)
            rt.gStyle.SetPaintTextFormat(cell_text_format)
            h.Draw("colz text e1")
            c.Print(outpath_pdf)
        for h, h_err in zip(hist_iterfit_tup, hist_err_iterfit_tup):
            set_bin_vals_and_errs(h, h_err)
            h.SetContour(n_contour)
            h.SetMarkerSize(text_size)
            # h.SetBarOffset(0.2)
            rt.gStyle.SetPaintTextFormat(cell_text_format)
            h.Draw("colz text e1")
            c.Print(outpath_pdf)
        c.Print(outpath_pdf + "]")