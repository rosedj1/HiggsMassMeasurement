"""Muon Collection pT Correction Factor Dict Updater

If a MyMuonCollection has the "old style" pT corr factor dict,
then it will be updated to be a 3-nested dict with:

fit stats for the dpT/pT vs. q*d0 corrections and
dpT/pT * 1/<pT> vs. q*d0 corrections.
"""
import sys
from Utils_Python.Utils_Files import open_pkl, save_to_pkl
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum #, bin_edges_pT_sevenfifths_to1000GeV_wholenum
from Utils_ROOT.Printer import CanvasPrinter

# printer = CanvasPrinter(show_plots=0)
# printer.make_plots_pretty()

inpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds.pkl"
newpath_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds_new.pkl"

mucoll = open_pkl(inpath_pkl)

mucoll.pT_corr_factor_dict = {
    "dpTOverpT_vs_qd0" : {},
    "dpTOverpTscaled_vs_qd0" : {},  # This is dpT/pT * 1/<pT> vs. qd0.
    "dpTOverpTtimesavgOf1divpT_vs_qd0" : {},
    "dpTOverpTtimesmuOf1divpT_vs_qd0" : {}
}
sys.exit("Finish filling out avgOf1divpT and muOf1divpT parts.")
for kb2d in mucoll.KinBin2D_dict.values():
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct = {
        "dpTOverpT_vs_qd0"       : {"interc_and_err" : None,
                                    "slope_and_err"  : None,
                                    "chi2"           : None,
                                    "NDF"            : None
        },
        "dpTOverpTscaled_vs_qd0" : {"interc_and_err" : None,
                                    "slope_and_err"  : None,
                                    "chi2"           : None,
                                    "NDF"            : None
        },
        # TODO: Finish filling out avgOf1divpT and muOf1divpT parts.
        # "dpTOverpTtimesavgOf1divpT_vs_qd0" : {"interc_and_err" : None,
        #                             "slope_and_err"  : None,
        #                             "chi2"           : None,
        #                             "NDF"            : None
        # },
        # "dpTOverpTtimesmuOf1divpT_vs_qd0" : {"interc_and_err" : None,
        #                             "slope_and_err"  : None,
        #                             "chi2"           : None,
        #                             "NDF"            : None
        # },
    }
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["interc_and_err"] = (kb2d.fit_line.GetParameter(0), kb2d.fit_line.GetParError(0))
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["slope_and_err"]  = (kb2d.fit_line.GetParameter(1), kb2d.fit_line.GetParError(1))
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["chi2"]  = kb2d.fit_line.GetChisquare()
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["NDF"]  = kb2d.fit_line.GetNDF()

    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["interc_and_err"] = (kb2d.fit_line_scaled.GetParameter(0), kb2d.fit_line_scaled.GetParError(0))
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["slope_and_err"]  = (kb2d.fit_line_scaled.GetParameter(1), kb2d.fit_line_scaled.GetParError(1))
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["chi2"]  = kb2d.fit_line_scaled.GetChisquare()
    kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["NDF"]  = kb2d.fit_line_scaled.GetNDF()

    bin_key = kb2d.get_bin_key()
    mucoll.pT_corr_factor_dict["dpTOverpT_vs_qd0"][bin_key] = {
        "intercept"     : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["interc_and_err"][0],
                              "intercept_err" : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["interc_and_err"][1],
                              "slope"         : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["slope_and_err"][0],
                              "slope_err"     : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpT_vs_qd0"]["slope_and_err"][1],
    }
    mucoll.pT_corr_factor_dict["dpTOverpTscaled_vs_qd0"][bin_key] = {
        "intercept"     : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["interc_and_err"][0],
                                    "intercept_err" : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["interc_and_err"][1],
                                    "slope"         : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["slope_and_err"][0],
                                    "slope_err"     : kb2d.dpTOverpT_vs_qd0_fit_stats_dct["dpTOverpTscaled_vs_qd0"]["slope_and_err"][1],
    }

save_to_pkl(mucoll, newpath_pkl, overwrite=overwrite)
#--- Draw all multigraphs ---#
# mucoll.make_all_multigraphs(equal_entry_bin_edges_eta_mod1_wholenum[:-1], scale_by_1divpT=False)
# printer.canv.Print(outpath_pdf + "[")
# mucoll.draw_all_multigraphs(outpath_pdf, printer, scale_by_1divpT=False)
# printer.canv.Print(outpath_pdf + "]")

