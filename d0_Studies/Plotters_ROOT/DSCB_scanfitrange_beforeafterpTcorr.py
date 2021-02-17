import ROOT
# # Local imports.
from d0_Studies.Plotters_ROOT.DSCB_scanner_relmasserr import make_fullpath_file
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid, fixOverlay
from Utils_Python.Utils_Files import check_overwrite
from Utils_ROOT.ROOT_fns import read_cpp_as_txt, load_cpp_code
from Utils_ROOT.ROOT_StatsAndFits import DSCBFitter, DSCBFitScanner, DSCBFitPlotter
from Utils_ROOT.ROOT_classes import make_TLegend, make_TMultiGraph_and_Legend
from Utils_ROOT.ROOT_Plotting import make_ratio_pads

# infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromDYJpsifactors_fullstats_noFSR.root"
infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromAdHocfactors_fullstats_noFSR_zerointerc.root"
# infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2017_m4mu_m4mucorrfromGeoFitfactors_fullstats_noFSR.root"
# infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2016_m4mu_m4mucorrfromGeoFitfactors_fullstats_noFSR.root"
outfile_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/plots/applypTcorrplots/CorrFromMC/tests/"
# infile_path = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/MC2018ggH_passFull_fullstats.root"
# outfile_dir = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/plots/test/DSCBscanoutput/"
outfile_prefix = "MCYEARggH_applypTcorrMETHOD_scanfitrange_zoom_zerointerc_test13"
# outfile_prefix = "MC2017ggH_applypTcorrGeoFit_scanfitrange"

#--- WARNING ---#
# When you change this, you must manually change the method in:
# Utils_ROOT/ROOT_StatsAndFits.py :: DSCBFitter.get_DSCBfit_and_draw()
# Looks something like: ROOT.fit_and_draw()
cpp_DSCB_code_path = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/fit_and_draw_DSCB_withAdHocCorr.C"
so_filepath = "./RooMyPDF_DSCB_C.so"
draw_beforeafter_corr = 1
overwrite = 1
verbose = 1
year = "2018"
method = "AdHoc" # "GeoFit"
zoom = False

m4mu_min, m4mu_max = 105.0, 145.0
relmasserr_min, relmasserr_max = 0.0, 5.0
max_lower_edge = 122.0
min_upper_edge = 128.0
step_GeV = 4
n_bins = 100
# m4mu_min = 121.2
# m4mu_max = 128.4

if __name__ == "__main__":
    print("...Preparing your area.")
    assert year in infile_path
    outfile_pdf = make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, relmasserr_min, relmasserr_max, step_GeV)
    outfile_pdf = outfile_pdf.replace("METHOD", method)
    outfile_pdf = outfile_pdf.replace("YEAR", year)
    outfile_pdf = outfile_pdf.replace(".txt", ".pdf")
    check_overwrite(outfile_pdf, overwrite=overwrite)
    tdrStyle = setTDRStyle(show_statsbox=True)
    tdrGrid(tdrStyle, gridOn=False)
    ROOT.gROOT.SetBatch(True)

    print("...Loading DSCB fitter code...")
    load_cpp_code(so_filepath, cpp_DSCB_code_path)
    print("...Opening root file...")
    f = ROOT.TFile(infile_path)
    tree = f.Get("tree")

    print("Using step_size = {} GeV".format(step_GeV))
    # print("Using bins:\n", fit_range_edges)
    print("...Preparing canvas.")
    canv = ROOT.TCanvas()
    canv.Print(outfile_pdf + "[")

    print("Performing DSCB fits.")
    scanner = DSCBFitScanner(relmasserr_min, relmasserr_max)  # Loop over these.
    scanner.do_DSCB_fits_over_mass4mu_range(m4mu_min, m4mu_max, step_GeV,
                                            max_lower_edge, min_upper_edge,
                                            tree, canv, outfile_pdf, n_bins,
                                            method=method, year=year,
                                            verbose=verbose,
                                            draw_beforeafter_corr=draw_beforeafter_corr,
                                            zoom=zoom)
    plotter = DSCBFitPlotter()
    # make_new_page_after is temporarily deprecated.
    gr_mean_uncorr = plotter.plot_X_vs_GeVfitrange("mean", scanner, canv, outfile_pdf, make_new_page_after=True)
    gr_mean_corr = plotter.plot_X_vs_GeVfitrange("mean", scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=True)
    gr_sigma_uncorr = plotter.plot_X_vs_GeVfitrange("sigma", scanner, canv, outfile_pdf, make_new_page_after=True)
    gr_sigma_corr = plotter.plot_X_vs_GeVfitrange("sigma", scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=True)
    gr_integral = plotter.plot_X_vs_GeVfitrange("integral", scanner, canv, outfile_pdf, make_new_page_after=True)
    gr_integral_corr = plotter.plot_X_vs_GeVfitrange("integral", scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=True)
    gr_sigma_improv = plotter.plot_X_vs_GeVfitrange("sigma_improv", scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=True)
    gr_mean_shift = plotter.plot_X_vs_GeVfitrange("mean_shift", scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=True)

    mg_mean, leg_mean = make_TMultiGraph_and_Legend(gr_ls=[gr_mean_uncorr, gr_mean_corr],
                                        leg_txt_ls=[r"Before p_{T} corr.", r"After p_{T} corr."],
                                        y_min=124.80, y_max=124.96)

    mg_sigma, leg_sigma = make_TMultiGraph_and_Legend(gr_ls=[gr_sigma_uncorr, gr_sigma_corr],
                                        leg_txt_ls=[r"Before p_{T} corr.", r"After p_{T} corr."],
                                        y_min=0.9, y_max=1.25)

    mg_integ, leg_integ = make_TMultiGraph_and_Legend(gr_ls=[gr_integral, gr_integral_corr],
                                        leg_txt_ls=[r"Before p_{T} corr.", r"After p_{T} corr."])

    mg_sigma_improv, leg_sigma_improv = make_TMultiGraph_and_Legend(gr_ls=[gr_sigma_improv],
                                        leg_txt_ls=[year],
                                        y_min=3.5, y_max=8.5)

    # mg_ls = [mg_mean, mg_sigma, mg_sigma_improv, mg_integ]
    # leg_ls = [leg_mean, leg_sigma, leg_sigma_improv, leg_integ]
    mg_ls = [mg_mean, mg_sigma]
    leg_ls = [leg_mean, leg_sigma]
    ratio_ls = [gr_mean_shift, gr_sigma_improv]
    text_size = 0.07
    # Use any graph to get min/max x-vals.
    x_min = min(gr_integral.x_vals) - 1
    x_max = max(gr_integral.x_vals) + 1
    for mg, leg, rat in zip(mg_ls, leg_ls, ratio_ls):
        y_err_min = (min(rat.y_vals) - 0.5 * max(rat.y_vals_err)) * 0.9
        y_err_max = (max(rat.y_vals) + 0.5 * max(rat.y_vals_err)) * 1.1
        # example_gr = list(mg.GetListOfGraphs())[0]
        # y_min = min(example_gr.y_vals)
        # y_max = max(example_gr.y_vals) + 1.1
        ptop, pbot = make_ratio_pads()  # Draws pads to canvas.
        ptop.Draw()
        pbot.Draw()
        ptop.cd()
        mg.Draw("a")
        ROOT.gPad.Modified()
        mg.GetXaxis().SetLimits(x_min, x_max)
        leg.Draw("same")
        pbot.cd()
        rat.GetXaxis().SetLabelSize(text_size)
        rat.GetYaxis().SetLabelSize(text_size)
        rat.GetXaxis().SetTitleSize(text_size)
        rat.GetYaxis().SetTitleSize(text_size)
        rat.GetYaxis().SetTitleOffset(0.6)  # 0.3, Default is 0.005.
        # rat.SetXTitle(x_title)
        # rat.SetYTitle("corr. / uncorr.")
        # rat.SetMinimum(0.8)
        # rat.SetMaximum(1.2)
        rat.GetYaxis().SetNdivisions(207)
        rat.GetXaxis().SetTickLength(0.12)
        rat.GetXaxis().SetLimits(x_min, x_max)
        rat.GetYaxis().SetRangeUser(y_err_min, y_err_max)
        rat.SetTitle("")
        rat.Draw()
        canv.Print(outfile_pdf)
        del ptop, pbot
        canv.Clear()
    for mg, leg in zip([mg_integ], [leg_integ]):
        mg.Draw("a")
        leg.Draw("same")
        canv.Print(outfile_pdf)
    canv.Print(outfile_pdf + "]")
    print("Done.")