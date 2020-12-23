import os
import ROOT as r
from array import array
import numpy as np
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid
from Utils_Python.Utils_Files import check_overwrite
from Utils_Python.Utils_StatsAndFits import prop_err_on_dsigoversig
from Utils_ROOT.ROOT_StatsAndFits import RooFit_iterative_gaus_fit

from d0_Studies.Toy_Model.toymodel.scripts.toss_toys import new_pixel_pos, pT_ls
from d0_Studies.Toy_Model.utils.classes import MuonTrack

tstyle = setTDRStyle()
tdrGrid(tstyle)
# r.gROOT.SetBatch(True)

outpath_dir = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/gauss_vs_manual_unc/"
prefix = "gauss_scatter"
suffix = "uniformuncert"
overwrite = False

sigma = 0.01
# Tracker uncertainties in um.
yerr_trkr_ls = [15, 15, 15, 15, 23, 23, 35, 35, 53, 53, 53, 53, 35, 35]
yerr_uniform_ls = [sigma*1E4]*len(new_pixel_pos)  # um.
pT = 5

mutrk = MuonTrack(pT, new_pixel_pos, yerr_trkr_ls, shift_yvals=(True, 0.01))
# Do BS constraint.
mutrk_BS = MuonTrack(pT, new_pixel_pos, yerr_trkr_ls, shift_yvals=(True, 0.01))
mutrk_BS.add_BS(pos=(0,0), bs_unc_y=10)

gr, leg = mutrk.plot_trajectory(leg=None, color=1)
gr_BS, leg_BS = mutrk.plot_trajectory(leg=None, color=1)

gr_BS.SetTitle( gr_BS.GetTitle() + "_withBS" )
c = r.TCanvas("c", "some title", 800, 800)
c.Divide(2,1)  # (2,1) is 2 cols, 1 row
c.cd(1)
gr.Draw("alp")
leg.Draw("same")
c.cd(2)
gr_BS.Draw("alp")
leg_BS.Draw("same")
c.Update()

# fitline, res = mutrk.fit_hits_pol2(gr, color=2)
# gr_bs, leg_bs = mutrk.plot_trajectory(leg=None, color=4)


# def do_itergausfit_

# n_bins_a_hist = 100
# n_bins_aerr_hist = 100

# a_ls = [FILLME]
# aerr_ls = [FILLME]
# d0_ls = [FILLME]
# d0err_ls = [FILLME]
# a_BS_ls = [FILLME]
# aerr_BS_ls = [FILLME]
# d0_BS_ls = [FILLME]
# d0err_BS_ls = [FILLME]

# a_hist_title = r"a_{fit} (p2) from y = ax^2 + bx + c"
# aerr_hist_title = r"#deltaa_{fit} (#deltap2)"
# r"{p_{T}^{reco,BS} - p_{T}^{gen}}{p_{T}^{gen}}"
# r"{p_{T}^{reco} - p_{T}^{gen}}{p_{T}^{gen}}"

# a_max = max(a_ls)
# aerr_max = max(aerr_ls)
# h_afit = r.TH1F("h_afit", a_hist_title, n_bins_a_hist, 0, a_max)
# h_afiterr = r.TH1F("h_afiterr", aerr_hist_title, n_bins_aerr_hist, 0, aerr_max)
# for a, aerr in zip(a_ls, aerr_ls):
#     h_afit.Fill(a)
#     h_afiterr.Fill(aerr)

# fit_stats_dict_a, xframe_a = RooFit_iterative_gaus_fit(a_ls, binned_fit=False, switch_to_binned_fit=100000, iters=2, num_sigmas=2.5, 
#                               n_bins=n_bins_a_hist, x_lim=None, fit_whole_range_first_iter=True,
#                               xframe=None, 
#                               x_label="Independent var", title="", units="",
#                               marker_color=1, force_last_line_color=None, only_draw_last=False, 
#                               verbose=True, view_plot=False)

# fit_stats_dict_aerr, xframe_aerr = RooFit_iterative_gaus_fit(aerr_ls, binned_fit=False, switch_to_binned_fit=100000, iters=1, num_sigmas=2.5, 
#                               n_bins=n_bins_aerr_hist, x_lim=None, fit_whole_range_first_iter=True,
#                               xframe=None, 
#                               x_label="Independent var", title="", units="",
#                               marker_color=1, force_last_line_color=None, only_draw_last=False, 
#                               verbose=True, view_plot=False)

# # Compare best-fit iter. Gaussian sigma of a_fit vals to best-fit mean
# # of a_fiterr dist.
# bf_sig_a_dist = fit_stats_dict_a["std_ls"][-1]
# bf_sigerr_a_dist = fit_stats_dict_a["std_err_ls"][-1]
# mu_aerr_dist = fit_stats_dict_aerr["mean_ls"][-1]
# muerr_aerr_dist = fit_stats_dict_aerr["mean_err_ls"][-1]
# print(f"best-fit sigma(a_fit_dist) +- err:   {bf_sig_a_dist:.6f} +- {bf_sigerr_a_dist:.6f}")
# print(f"best-fit mean(aerr_fit_dist) +- err: {mu_aerr_dist:.6f} +- {muerr_aerr_dist:.6f}")
# compare_val = (bf_sig_a_dist - mu_aerr_dist) / mu_aerr_dist * 100.
# compare_err = prop_err_on_dsigoversig(mu_aerr_dist, bf_sig_a_dist, 
#                                       muerr_aerr_dist, bf_sigerr_a_dist)

# sig_txt = r"#sigma(a)"
# mu_txt = r"#mu(#deltaa)"
# txt = r"#left(#frac{%s - %s}{%s}#right)#times100%" % (sig_txt, mu_txt, mu_txt)
# txt += r" = %.3f#pm%.3f" % (compare_val, compare_err)
# pave = r.TPaveText(0.70, 0.40, 0.90, 0.55, "NDC")
# pave.SetFillColor(0)
# pave.SetBorderSize(0) # Use 0 for no border.
# pave.SetTextAlign(22) # 22 is centered vert and horiz.
# pave.SetTextSize(0.016)
# pave.SetTextColor(1)
# pave.SetFillStyle(1001)  # Solid fill.
# pave.AddText(txt)  # Accommodates LaTeX!

# c = r.TCanvas()
# c.Divide(2,1)  # (2,1) is 2 cols, 1 row
# c.Print(outpath_file + "[")
# c.cd(1)
# xframe_a.Draw()
# c.cd(2)
# xframe_aerr.Draw()
# pave.Draw("same")

# def make_filename(outpath_dir, sigma, hits_ls, yerr_ls=None, prefix=None, suffix=None):
#     """Return an absolute path (str) for final output file.
    
#     NOTE: sigma should be given in cm.
#     """
#     sig_title = str(sigma*1E4)
#     nhits = len(hits_ls)
#     filename = f"allpT_{nhits}hits_{sig_title}umsmear"
#     if prefix is not None:
#         filename = f"{prefix}_{filename}"
#     if (yerr_ls is not None) and (len(set(yerr_ls)) == 1):
#         # Take any element; they're all the same.
#         newend = str(yerr_ls[0])
#         filename = f"{filename}_{newend}"
#     if suffix is not None:
#         filename = f"{filename}_{suffix}"
#     filename = filename.replace('.', 'p') + ".pdf"
#     outpath_file = os.path.join(outpath_dir, filename)
#     return outpath_file

# def do_comparison_GaussErrstoManualErrs(hits_ls, pT, sigma, yerr_ls):
#     """Return a 4-tuple of graphs and fit results comparing y-errs on muon
#     hits as taken from G(mu=0, sigma) vs. manually specified (yerr_ls).
    
#     Parameters
#     ----------
#     hits_ls : list
#         The x-positions (cm) of where the muon hits tracker.
#     pT : float
#         Transverse momentum of muon in GeV.
#     sigma : float
#         The standard deviation of a Gaussian from which y-val errors will be
#         taken. Should be given in cm.
#     yerr_ls : list
#         The "manual" y-val errors to be compared to the Gaussian ones.

#     Returns
#     -------
#     4-tuple (gr, gr_witherrs, fit_func, fit_func_witherrs)
#     """
#     # yerr_smear_ls = [sigma * 1E4] * len(yerr_ls)  # Convert to um.

#     # Prep.
#     assert len(hits_ls) == len(yerr_ls)
#     yerr_ls = np.array(yerr_ls) / 1E4  # cm.

#     # Make the graphs.
#     hitplotorg = HitPlotOrg(pT, hits_ls)
#     hitplotorg_witherrs = HitPlotOrg(pT, hits_ls)

#     gr, leg = hitplotorg.plot_hit_trajectory(leg=None, color=1, smear=True, sigma=sigma)
#     gr_witherrs, leg_witherrs = hitplotorg_witherrs.plot_hit_trajectory(leg=None, color=1, yerr_ls=yerr_ls)

#     for g in [gr, gr_witherrs]:
#         g.SetMaximum(14.5)
#         g.SetMinimum(-0.2)

#     # Fit the graphs.
#     fit_func = hitplotorg.fit_hits_pol2(gr)
#     fit_func_witherrs = hitplotorg_witherrs.fit_hits_pol2(gr_witherrs)

#     return (gr, gr_witherrs, fit_func, fit_func_witherrs, leg, leg_witherrs)

# def make_pdf_all_pTs_GaussvsManuncert(hits_ls, pT_ls, sigma, yerr_ls, outpath_dir, prefix=None, suffix=None):
#     """Make a multi-page PDF with Gauss uncert graphs vs. manual uncert."""
#     outpath_file = make_filename(outpath_dir, sigma=sigma, hits_ls=hits_ls, prefix=prefix, suffix=suffix)
#     check_overwrite(outpath_file, overwrite)
#     c = r.TCanvas()
#     c.Divide(1,2)  # (2,1) is 2 cols, 1 row
#     c.Print(outpath_file + "[")
#     for pT in pT_ls:
#         obj_tup = do_comparison_GaussErrstoManualErrs(hits_ls, pT, sigma, yerr_ls)
#         gr                 = obj_tup[0]
#         gr_witherrs        = obj_tup[1]
#         fit_func           = obj_tup[2]
#         fit_func_witherrs  = obj_tup[3]
#         leg                = obj_tup[4]
#         leg_witherrs       = obj_tup[5]

#         leg.AddEntry(gr, r"Hits in Tracker, shifted by G(#mu=0, #sigma=%.0f #mum)" % (sigma*1E4), "lp")
#         leg.AddEntry(fit_func, r"y = p0 + p1*x + p2*x^{2}", "l")
#         leg_witherrs.AddEntry(gr, r"Hits in Tracker using Tracker uncert.", "lpe")
#         leg_witherrs.AddEntry(fit_func, r"y = p0 + p1*x + p2*x^{2}", "l")

#         # Plot everything.
#         c.cd(1)
#         gr.Draw("ALP")
#         fit_func.Draw("same")
#         leg.Draw("same")
#         c.cd(2)
#         gr_witherrs.Draw("ALP")
#         fit_func_witherrs.Draw("same")
#         leg_witherrs.Draw("same")
#         c.Update()
#         c.Print(outpath_file)
#     c.Print(outpath_file + "]")

# def make_pdf_all_pT_Gaussuncert(hits_ls, pT_ls, sigma, outpath_dir, prefix=None, suffix=None):
#     """Make a multi-page PDF with Gauss uncert graphs."""
#     outpath_file = make_filename(outpath_dir, sigma=sigma, hits_ls=hits_ls, prefix=prefix, suffix=suffix)
#     check_overwrite(outpath_file, overwrite)
#     print(f"Making multi-page PDF for Gaussian smears in file:\n{outpath_file}")
#     # Make plots.
#     c = r.TCanvas()
#     c.Print(outpath_file + "[")
#     for pT in pT_ls:
#         hitplotorg = HitPlotOrg(pT, hits_ls)
#         (gr, leg) = hitplotorg.plot_hit_trajectory(leg=None, color=1, smear=True, sigma=sigma)
#         fit_func = hitplotorg.fit_hits_pol2(gr)
#         gr.SetMaximum(14.5)
#         gr.SetMinimum(-0.2)
#         # gr, gr_witherrs, fit_func, fit_func_witherrs = do_comparison_GaussErrstoManualErrs(new_pixel_pos, pT, sigma, yerr_ls)
#         gr.Draw("ALP")
#         fit_func.Draw("same") 
#         leg.AddEntry(gr, r"Hits in Tracker, shifted by G(#mu=0, #sigma=%.0f #mum)" % sigma, "lp")
#         leg.AddEntry(fit_func, r"y = p0 + p1*x + p2*x^{2}", "l")
#         leg.Draw("same")
#         c.Update()
#         c.Print(outpath_file)
#     c.Print(outpath_file + "]")

# if __name__ == "__main__":
    # make_pdf_all_pTs_GaussvsManuncert(new_pixel_pos, pT_ls, sigma, yerr_uniform_ls, outpath_dir, prefix=prefix, suffix=suffix)
    # make_pdf_all_pT_Gaussuncert(new_pixel_pos, pT_ls, sigma, outpath_dir, prefix=prefix, suffix=suffixe)
