"""Template for Iterative Gaussian Fits on Rochester Corr. vs. Non-RC samples

Requires an input pickled dict of KinBin2Ds.
- You can make the input dict from roch_vs_noroch_kb2dmaker_inclusivehistplotter.py

This is a template file whose values (IN ALL CAPS) get replaced by a SLURM
submission script. The iterative fits are computationally intensive,
and so you should use the following files to control this one:
- d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py
- 

Author: Jake Rosenzweig
Updated: 2021-02-05
"""
from Utils_Python.Utils_Files import open_pkl, save_to_pkl
from Utils_Python.Utils_Files import check_overwrite
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid
from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit
import ROOT
import os
import gc
#----- User Switches -----#
overwrite = 1
verbose = 0
iters = 5
fit_whole_range_first_iter = False  # False gives much more consistent fits.
# eta_ls = equal_entry_bin_edges_eta_mod1_wholenum #[0.0, 0.9, 1.7, 2.4]
# pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum#[5, 25, 50, 100]
inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/MC2018ggH_KB2D_onlymuoninfo_fullstats_pTbinsroundnumbers_75then200GeV.pkl"
# Outpath info is controlled by d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py
eta_range = REPLACE_ETA_LS
# eta_range = [0.0, 0.2]

if __name__ == "__main__":
    eta_str = f"{eta_range[0]}eta{eta_range[1]}".replace(".", "p")
    filename = f"REPLACE_JOB_NAME_{iters}iters_{eta_str}"
    outpath_pkl = os.path.join("REPLACE_OUTDIR_PKL", f"{filename}.pkl")
    outpath_txt = os.path.join("REPLACE_OUTDIR_TXT", f"{filename}.txt")
    outpath_pdf = os.path.join("REPLACE_OUTDIR_PDF", f"{filename}.pdf")
    # f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/output/{filename}.txt"
    check_overwrite(outpath_pkl, overwrite)
    check_overwrite(outpath_txt, overwrite)
    check_overwrite(outpath_pdf, overwrite)
    tdrStyle = setTDRStyle(show_statsbox=True)
    tdrGrid(tdrStyle, gridOn=False)
    ROOT.gROOT.SetBatch(True)
    # org_kb2d = OrganizerKB2D()
    # org_kb2d.read_KB2D_from_pkl(inpath_pkl)
    muon_coll = open_pkl(inpath_pkl)
    new_dct = {}
    for kb2d in muon_coll.KinBin2D_dict.values():
        if kb2d.eta_range == eta_range:
            print(f"...Doing iter Gaus fit on: eta_range={kb2d.eta_range}, pT_range={kb2d.pT_range}")
            # RochCorr vs. reco
            kb2d.dpT_RC_reco_fit_dct, kb2d.dpT_RC_reco_frame = RooFit_iterative_gaus_fit(
                                    data=[(mu.pT_RC - mu.pT) * 1000.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-500,500], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"p_{T}^{RC} - p_{T}^{reco}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="MeV", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose)
            # Scaled by pT.
            kb2d.dpTOverpT_RC_reco_fit_dct, kb2d.dpTOverpT_RC_reco_frame = RooFit_iterative_gaus_fit(
                                    data=[((mu.pT_RC - mu.pT)/mu.pT) * 100.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-2,2], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"(p_{T}^{RC} - p_{T}^{reco})/p_{T}^{reco}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="%", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose)
            # RochCorr vs. gen
            kb2d.dpT_RC_gen_fit_dct, kb2d.dpT_RC_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[(mu.pT_RC - mu.gen_pT) for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-20,20], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"p_{T}^{RC} - p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="GeV", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose)
            # Scaled by pT.
            kb2d.dpTOverpT_RC_gen_fit_dct, kb2d.dpTOverpT_RC_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[((mu.pT_RC - mu.gen_pT)/mu.gen_pT) * 100.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-10,10], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"(p_{T}^{RC} - p_{T}^{gen})/p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="%", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose)
            # reco vs. gen
            kb2d.dpT_reco_gen_fit_dct, kb2d.dpT_reco_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[(mu.pT - mu.gen_pT) for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-20,20], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"p_{T}^{reco} - p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="GeV", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose)
            # Scaled by pT.
            kb2d.dpTOverpT_reco_gen_fit_dct, kb2d.dpTOverpT_reco_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[((mu.pT - mu.gen_pT)/mu.gen_pT) * 100.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-10,10], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="%", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose)
            key = kb2d.get_bin_key()
            new_dct[key] = kb2d
            gc.collect()
    # End loop over KB2Ds.
    def make_pdf_kb2d_frames(dct, attr, outpath_pdf):
        """
        Draw the attribute (an xframe) of each kb2d in dct to a canvas
        and then close the canvas.
        
        attr : str
        """
        c = ROOT.TCanvas()
        c.Print(outpath_pdf + "[")
        for kb2d in dct.values():
            try:
                getattr(kb2d, attr).Draw()
                c.Print(outpath_pdf)
            except:
                print("The following kb2d cannot draw:")
                print(f"  eta={kb2d.eta_range}, pT={kb2d.pT_range}")
                print(f"  Its __dict__ is:\n{kb2d.__dict__}")
        c.Print(outpath_pdf + "]")
    
    # Have all KB2Ds draw each attribute, then close the PDF.
    # Then move on to next attr.
    attr_tup = ("dpT_RC_reco_frame", "dpTOverpT_RC_reco_frame",
                "dpT_RC_gen_frame", "dpTOverpT_RC_gen_frame",
                "dpT_reco_gen_frame", "dpTOverpT_reco_gen_frame")
    for attr in attr_tup:
        suff = attr.rstrip("_frame")
        pdf_path = outpath_pdf.replace(".pdf", f"{suff}.pdf")
        make_pdf_kb2d_frames(new_dct, attr, pdf_path)
    save_to_pkl(new_dct, outpath_pkl)
        