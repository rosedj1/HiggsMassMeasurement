"""Template for Iterative Gaussian Fits on Rochester Corr. vs. Non-RC samples

Requires an input pickled dict of KinBin2Ds.
- You can make the input dict from roch_vs_noroch_kb2dmaker_inclusivehistplotter.py

This is a template file whose values (IN ALL CAPS) get replaced by a SLURM
submission script. The iterative fits are computationally intensive,
and so you should use the following files to control this one:
- d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py
- 

Author: Jake Rosenzweig
Updated: 2021-02-24
"""
from Utils_Python.Utils_Files import open_pkl, save_to_pkl
from Utils_Python.Utils_Files import check_overwrite
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid
from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit
import ROOT
import os
import gc
#----- User Switches -----#
overwrite = 0
verbose = 0
iters = 5
verbose = 1
fit_whole_range_first_iter = False  # False gives much more consistent fits.
use_data_in_xlim = True
# eta_ls = equal_entry_bin_edges_eta_mod1_wholenum #[0.0, 0.9, 1.7, 2.4]
# pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum#[5, 25, 50, 100]
# inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/MC2018ggH_KB2D_onlymuoninfo_fullstats_pTbinsroundnumbers_75then200GeV.pkl"
inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/test/fixitergaussfits_test01_3kbds.pkl"
#--- Output ---#
outdir_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/pickles/test/"
outdir_txt = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/output/test/"
outdir_pdf = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/RochCorr/plots/test/"
filename = "RC_vs_NoRC_itergaussfits_fullstats_fixitergauss_test05_definitelyunbinned"
eta_range = [0.2, 0.4]

if __name__ == "__main__":
    eta_str = f"{eta_range[0]}eta{eta_range[1]}".replace(".", "p")
    filename = f"{filename}_{iters}iters_{eta_str}"
    outpath_pkl = os.path.join(outdir_pkl, f"{filename}.pkl")
    outpath_txt = os.path.join(outdir_txt, f"{filename}.txt")
    outpath_pdf = os.path.join(outdir_pdf, f"{filename}.pdf")
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
    for kb2d in muon_coll.values():
        if kb2d.eta_range == eta_range:
            print(f"...Doing iter Gaus fit on: eta_range={kb2d.eta_range}, pT_range={kb2d.pT_range}")
            # RochCorr vs. reco
            # pT_RC_arr = [mu.pT_RC for mu in kb2d.muon_ls]
            # pT_arr = [mu.pT for mu in kb2d.muon_ls]
            # data = [mu.pT_RC for mu in kb2d.muon_ls]

            # dpT_RC_reco = (pT_RC_arr - pT_arr) * 1000.0
            # dpT_RC_reco[(dpT_RC_reco > x_lim)]

            kb2d.dpT_RC_reco_fit_dct, kb2d.dpT_RC_reco_frame = RooFit_iterative_gaus_fit(
                                    data=[(mu.pT_RC - mu.pT) * 1000.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-500,500], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"p_{T}^{RC} - p_{T}^{reco}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="MeV", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose,
                                    use_data_in_xlim=use_data_in_xlim)
            # Scaled by pT.
            kb2d.dpTOverpT_RC_reco_fit_dct, kb2d.dpTOverpT_RC_reco_frame = RooFit_iterative_gaus_fit(
                                    data=[((mu.pT_RC - mu.pT)/mu.pT) * 100.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=5000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-2,2], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"(p_{T}^{RC} - p_{T}^{reco})/p_{T}^{reco}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="%", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose,
                                    use_data_in_xlim=False)
            # RochCorr vs. gen
            kb2d.dpT_RC_gen_fit_dct, kb2d.dpT_RC_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[(mu.pT_RC - mu.gen_pT) for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-20,20], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"p_{T}^{RC} - p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="GeV", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose,
                                    use_data_in_xlim=use_data_in_xlim)
            # Scaled by pT.
            kb2d.dpTOverpT_RC_gen_fit_dct, kb2d.dpTOverpT_RC_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[((mu.pT_RC - mu.gen_pT)/mu.gen_pT) * 100.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-10,10], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"(p_{T}^{RC} - p_{T}^{gen})/p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="%", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose,
                                    use_data_in_xlim=use_data_in_xlim)
            # reco vs. gen
            kb2d.dpT_reco_gen_fit_dct, kb2d.dpT_reco_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[(mu.pT - mu.gen_pT) for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-20,20], fit_whole_range_first_iter=fit_whole_range_first_iter,
                                    xframe=None, x_label=r"p_{T}^{reco} - p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="GeV", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose,
                                    use_data_in_xlim=use_data_in_xlim)
            # Scaled by pT.
            kb2d.dpTOverpT_reco_gen_fit_dct, kb2d.dpTOverpT_reco_gen_frame = RooFit_iterative_gaus_fit(
                                    data=[((mu.pT - mu.gen_pT)/mu.gen_pT) * 100.0 for mu in kb2d.muon_ls], 
                                    binned_fit=False, switch_to_binned_fit=2000, iters=iters, num_sigmas=2.5,
                                    n_bins=100, x_lim=[-10,10], fit_whole_range_first_iter=True,
                                    xframe=None, x_label=r"(p_{T}^{reco} - p_{T}^{gen})/p_{T}^{gen}", title="%s" % kb2d.make_latex_bin_cut_str(), 
                                    units="%", marker_color=1, force_last_line_color=None, only_draw_last=False, verbose=verbose,
                                    use_data_in_xlim=use_data_in_xlim)
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
        pdf_path = outpath_pdf.replace(".pdf", f"_{suff}.pdf")
        make_pdf_kb2d_frames(new_dct, attr, pdf_path)
    save_to_pkl(new_dct, outpath_pkl)