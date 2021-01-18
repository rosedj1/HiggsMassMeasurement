from __future__ import print_function
import ROOT
import os
import sys
from array import array
from pprint import pprint
from collections import OrderedDict
import numpy as np
# Local imports.
sys.path.append("/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/")
from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import setTDRStyle, tdrGrid, fixOverlay
from Utils_ROOT.ROOT_classes import make_TH1F, make_TGraphErrors, make_TLegend
# from d0_Studies.d0_Utils.d0_fns import calc_x_err_bins_from_bin_edges

infile_path = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/MC2018ggH_passFull_fullstats.root"
outfile_dir = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/plots/test/DSCBscanoutput/"
outfile_prefix = "MC2018ggH_singlemasserrorbin_saveintegral"
cpp_DSCB_code_path = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/Plotters_ROOT/fit_and_draw_DSCB.C"

overwrite = 0
m4mu_min, m4mu_max = 105.0, 140.0
relmasserr_min, relmasserr_max = 0.0, 0.75
m4muErr_binedge_ls = [0.0, 0.75, 1.25, 5.0]  # This is a percent!
# m4muErr_binedge_ls = [0.0, 0.75]
max_lower_edge = 116.5
min_upper_edge = 128.0
step_GeV = 0.5
n_bins = 100
# m4mu_min = 121.2
# m4mu_max = 128.4

def make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, relmasserr_min, relmasserr_max, step_GeV):
    f = "{}_{}m4mu{}_step{:.1f}GeV_{}relmasserr{}".format(outfile_prefix, m4mu_min, m4mu_max, step_GeV, relmasserr_min, relmasserr_max)
    outfile_name = f.replace(".", "p") + ".txt"
    return os.path.join(outfile_dir, outfile_name)

def check_overwrite(outfile, overwrite=False):
    """
    Raises an error if outfile exists and overwrite == False. 
    """
    if os.path.exists(outfile) and not (overwrite):
        err_msg = """
            Not allowed to overwrite file since it already exists
            \n{}\n
            To write over this file, set overwrite = True.\n
            """.format(outfile)
        raise RuntimeError(err_msg)

def read_cpp_as_txt(cpp_file):
    """Return a string of all the C++ code store in cpp_file."""
    with open(cpp_file, "r") as f:
        line_ls = f.readlines()
        # Each element in the list is proceeded by '\n'.
        clean_line_ls = [x.rstrip('\n') for x in line_ls]
        # Turn into a single string:
        code_str = '\n'.join(clean_line_ls)
        return code_str

def load_cpp_code(dot_so_filepath, cpp_code_path):
    """Load the C++ code stored in cpp_code_path and a corresponding .so file."""
    code_str = read_cpp_as_txt(cpp_DSCB_code_path)
    ROOT.gSystem.Load(dot_so_filepath)
    ROOT.gInterpreter.Declare(code_str)

#***** (Probably) useless functions below *****#
def get_params_and_errs(tree, m4mu_min, m4mu_max):
    """Return mean and mean_err from DSCB on m4mu range.
    
    FIXME: Makes a call to deprecated fn: ROOT.DSCB_unbinnedfit
    """
    print("Performing DSCB unbinned fit.")
    result = ROOT.DSCB_unbinnedfit(t, m4mu_min, m4mu_max)
    result_tup = parse_DSCBfit_result(result)
    return result_tup

def make_fit_dict(fit_range_edges, max_lower_edge, min_upper_edge):
    """Return a dict of fit stats from unbinned DSCB...
    
    FIXME: Makes a call to deprecated fn: ROOT.DSCB_unbinnedfit
    """
    fit_dct = OrderedDict()
    # Use the lowest vales in fit_range_edges and
    # the highest vals in its flipped array.
    for minval, maxval in zip(fit_range_edges, fit_range_edges[::-1]):
        if (minval >= max_lower_edge) or (maxval < min_upper_edge):
            break
        print("\n...Filling dict over range: [{}, {}]".format(minval, maxval))
        # mean, mean_err, sigma, sigma_err = get_params_and_errs(t, minval, maxval)
        key = (minval, maxval)
        fit_dct[key] = {"mean":mean, "sigma":sigma, "mean_err":mean_err, "sigma_err":sigma_err}
    return fit_dct

def dscb_fits_GeVsteps(infile_path, outfile_dir, outfile_prefix, m4mu_min, m4mu_max, step_GeV, max_lower_edge, bins=100, overwrite=False):
    """
    Make a txt file with the fit stats from multiple DSCB fits on m4mu dist.
    Subsequent fits will narrow the fit range by `step_GeV`.
    """
    # outfile_txt = make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, step_GeV)
    # check_overwrite(outfile_txt, overwrite)
    # print("Creating file:\n", outfile_txt)
    # print("Loading C++ DSCB unbinned fitting code...")
    # load_cpp_code()
    # print("Opening file:\n", infile_path)
    # f = ROOT.TFile(infile_path)
    # t = f.Get("passedEvents")

    bins = (m4mu_max - m4mu_min) / float(step_GeV)
    fit_range_edges = np.linspace(m4mu_min, m4mu_max, bins+1)
    print("Using step_size = {} GeV".format(step_GeV))
    print("Using bins:\n", fit_range_edges)
    # lower_bounds = np.linspace(105, 120, 4)
    # upper_bounds = np.linspace(130, 140, 4)[::-1]  # Flip it around to zip properly.
    fit_dct = make_fit_dict(t, fit_range_edges, max_lower_edge)
    pprint(fit_dct)
    write_output_file(outfile_txt, fit_dct)
    # End loop over bounds.
    #print("\nFinal means and errors:")
    #pprint(fit_list)
#***** (Probably) useless functions above *****#
def write_output_file(outfile_txt, fit_dct):
    with open(outfile_txt, "w") as f:
        header = "fit_range_min\tfit_range_max\tmean\tmean_err\tsigma\tsigma_err\n"
        f.write(header)
        for key,val in fit_dct.items():
            tup_of_keys = (key[0], key[1], val["mean"], val["mean_err"], val["sigma"], val["sigma_err"])
            print("...Writing these vals to file:", tup_of_keys)
            line = "{:13.8f}\t{:13.8f}\t{:13.8f}\t{:13.8f}\t{:13.8f}\t{:13.8f}\n".format(*tup_of_keys)
            f.write(line)

def make_and_fill_hists(t):
    """
    Return a list of 3 filled histograms:
    - mass4mu
    - mass4lErr
    - mass4mu/mass4lErr * 100%
    """
    # Make histograms.
    h_m4mu = make_TH1F(r"h_m4mu", title=r"Mass", xlabel=r"m_{4#mu}", n_bins=140, x_min=105, x_max=140, units=r"GeV")
    h_dm4mu = make_TH1F(r"h_dm4mu", title=r"Mass Error", xlabel=r"#deltam_{4#mu}", n_bins=100, x_min=0, x_max=5, units=r"GeV")
    h_reldm4mu = make_TH1F(r"h_reldm4mu", title=r"Relative Mass Error", xlabel=r"#deltam_{4#mu}/m_{4#mu} (%)", n_bins=100, x_min=0, x_max=5.0)
    cuts = "passedFullSelection==1 && finalState==1"
    # Fill histograms.
    t.Draw("mass4l >> h_m4mu", cuts, "goff")
    t.Draw("mass4lErr >> h_dm4mu", cuts, "goff")
    t.Draw("mass4lErr/mass4l*100 >> h_reldm4mu", cuts, "goff")
    return [h_m4mu, h_dm4mu, h_reldm4mu]

def make_pdf_mass4mu_mass4lErr_relmass4lErr(tree):
    """Create a PDF with 3 dists: mass4mu, mass4lErr, their ratio."""
    print("...Making and filling hists.")
    hist_ls = make_and_fill_hists(tree)
    print("...Hists made. Making plots.")

    canv = ROOT.TCanvas()
    canv.Print(outfile + "[")
    for h in hist_ls:
        h.Draw("hist")
        canv.Print(outfile)
    canv.Print(outfile + "]")
    print("Done.")

class DSCBFitter:
    """Organizer for an unbinned single DSCB fit."""

    def __init__(self):
        """Instantiate a DSCB fit object."""
        self.m4mu_min = None
        self.m4mu_max = None
        self.relm4muerr_min = None
        self.relm4muerr_max = None
        self.n_bins = None
        self.fit_result = None
        self.integral = None

    def get_DSCBfit_and_draw(self, tree, canv, outfile_pdf,
        m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max, n_bins):
        """Draw a DSCB fit to a canvas, store the info, and return the fit results.
        
        NOTE: Make sure you have an open TCanvas on which to draw.

        Parameters
        ----------
        tree : TTree
            Contains branches: mass4mu, mass4muErr, passedFullSelection, finalState
        """
        # print("    * Filling dict over m4mu range: [{:.3f}, {:.3f}]".format(m4mu_min, m4mu_max))
        print("  * Performing DSCB fit: ")
        self.info_printer(m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max, n_bins)
        # Perform selections and draw to canvas.
        integral = ROOT.vector('Double_t')()
        result = ROOT.fit_and_draw_DSCB(tree, m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max, integral, canv, outfile_pdf, n_bins)
        self.m4mu_min = m4mu_min
        self.m4mu_max = m4mu_max
        self.relm4muerr_min = relm4muerr_min
        self.relm4muerr_max = relm4muerr_max
        self.n_bins = n_bins
        self.fit_result = result
        self.integral = integral[0]
        return result

    def info_printer(self, m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max, n_bins):
        """Print useful debug info."""
        print("  * Unbinned fit show on plot with n_bins = {}".format(n_bins))
        print("  * mass4mu range =         [{}, {}]".format(m4mu_min, m4mu_max))
        print("  * rel_mass4mu_err range = [{}, {}]".format(relm4muerr_min, relm4muerr_max))

    def parse_DSCBfit_result(self, res):
        """
        Return a 4-tuple of the fit results from a DSCB fit.
        Store the fit result info.
        
        Parameters
        ----------
        res : RooFitResult* (pointer to a RooFitResult)
        """
        try:
            mean = res.floatParsFinal().find("mean").getVal()
            mean_err = res.floatParsFinal().find("mean").getError()
            sigma = res.floatParsFinal().find("sigma").getVal()
            sigma_err = res.floatParsFinal().find("sigma").getError()
        except AttributeError:
            mean = res.floatParsFinal().find("#mu").getVal()
            mean_err = res.floatParsFinal().find("#mu").getError()
            sigma = res.floatParsFinal().find("#sigma").getVal()
            sigma_err = res.floatParsFinal().find("#sigma").getError()
        self.mean = mean
        self.mean_err = mean_err
        self.sigma = sigma
        self.sigma_err = sigma_err
        return (mean, mean_err, sigma, sigma_err)
        
class DSCBFitScanner:
    """
    Class to organize various DSCB fits across different mass4mu ranges.
    NOTE: Fits are all within a single rel_m4mu_err bin.
    """
    def __init__(self, relm4muerr_min, relm4muerr_max):
        self.relm4muerr_min = relm4muerr_min
        self.relm4muerr_max = relm4muerr_max
        self.dscb_ls = []

    def do_DSCB_fits_over_mass4mu_range(self, x_min, x_max, stepsize,
                                              max_lower_edge, min_upper_edge,
                                              tree, canv, outfile_pdf, n_bins):
        """
        Draw multiple DSCB fits over a mass4mu range to a TCanvas.

        Example for explanation of parameters:
            fit_range_edges = [1,2,3,4,5,6,7,8]
            So fit1 has starting range: [1,8]
            So fit2 has starting range: [2,7]
            So fit3 has starting range: [3,6]
        Parameters
        ----------
        x_min : float
            The starting point of your binning range.
        x_max : float
            The ending point of your binning range.
        stepsize : float
            The bin width.
        max_lower_edge : float
            An x value which terminates the fitting procedure.
            Triggers when the (increasing) lower bound of the fit > max_lower_edge.
            In the above example, if max_lower_edge == 2, then fit3 would not happen.
        min_upper_edge : float
            An x value which terminates the fitting procedure.
            Triggers when the (decreasing) upper bound of the fit < min_upper_edge.
            In the above example, if min_upper_edge == 8, then fits2,3 would not happen.
        tree : TTree
            Contains branches: mass4mu, mass4muErr, passedFullSelection, finalState
        canv : ROOT.TCanvas
            The canvas on which to draw the plots.
        outfile_pdf : str
            The full filepath to store the pdf.
        n_bins : int
            The number of bins in which to bin the mass4mu distribution.
        """
        # Make sure n_bins can be represented as an int.
        if not isinstance(n_bins, int):
            assert (n_bins).is_integer() 
        fit_range_edges = self.make_binedges_from_stepsize(x_min, x_max, stepsize)
        print("This scanner made fit_range_edges:\n", fit_range_edges)
        for m4mu_min, m4mu_max in zip(fit_range_edges, fit_range_edges[::-1]):
            if (m4mu_min >= max_lower_edge) or (m4mu_max < min_upper_edge):
                break
            dscb = DSCBFitter()
            res = dscb.get_DSCBfit_and_draw(tree, canv, outfile_pdf,
                          m4mu_min, m4mu_max, self.relm4muerr_min, self.relm4muerr_max, n_bins)
            mean, sigma, mean_err, sigma_err = dscb.parse_DSCBfit_result(res)
            self.dscb_ls.append(dscb)

    def calc_num_bins(self, x_min, x_max, stepsize, allow_rebinning=False):
        """Return the number of bins from x_min to x_max in bin widths of `stepsize`.
        
        NOTE: 
          - If the parameters do not yield an int number of bins,
            then the number of bins will be rounded to nearest int.
          - This fn will not change x_min or x_max, but you are not guaranteed
            to get the stepsize you requested

        Parameters
        ----------
        x_min : float
            The starting point of your binning range.
        x_max : float
            The ending point of your binning range.
        stepsize : float
            The bin width.
        allow_rebinning : bool
            If True, then will round stepsize to the nearest int to make sure
            there are an int number of bins between x_min, x_max.
        """
        bins = (x_max - x_min) / float(stepsize)
        try:
            assert bins.is_integer()
        except AssertionError:
            msg = "  Cannot create int number of bins using:\n"
            msg += "  x_min={}, x_max={}, stepsize={}".format( x_min, x_max, stepsize)
            if allow_rebinning:
                print("...Warning! You will not get the stepsize you requested")
                print(msg)
                bins = round(bins)
            else:
                raise ValueError(msg)
        return bins

    def make_binedges_from_stepsize(self, x_min, x_max, stepsize):
        """Return an array of bin edges: [x_min, x_max, stepsize]."""
        bins = self.calc_num_bins(x_min, x_max, stepsize, allow_rebinning=True)
        assert bins.is_integer()
        return np.linspace(x_min, x_max, bins+1)

            # key = (m4mu_min, m4mu_max)
            # self.m4mu_fit_result_dct[key] = {"mean":mean,
            #                                  "sigma":sigma,
            #                                  "mean_err":mean_err,
            #                                  "sigma_err":sigma_err,
            #                                  "fit_result":res
            #                                  }
        # End loop over m4mu bins.

class DSCBFitPlotter:
    """Class to show how DSCB fit parameters change over different fit ranges."""
    
    def __init__(self):
        pass
        # self.relmass4mu_ls = relmass4mu_ls
    
    def plot_X_vs_GeVfitrange(self, var, scanner, canv, outfile_pdf, make_new_page_after=True):
        """Draw plots of var vs. different GeV fit ranges,
        where var is something like: mean(DSCB) and sigma(DSCB).

        Example:
            Fit 1 range: [105, 140] GeV
            Fit 2 range: [110, 135] GeV

        Parameters
        ----------
        var : str
            Attribute of DSCBFitter() object. Will be plotted as the y-coord.
            Supports: "mean", "sigma", "integral"
        scanner : DSCBFitScanner()
            FIXME: Finish docstring.
        """
        # Get the data from the fits.
        x = [dscb.m4mu_min for dscb in scanner.dscb_ls]
        y_mean = [getattr(dscb, var) for dscb in scanner.dscb_ls]  # Either dscb.mean or .sigma.
        y_mean_err = [getattr(dscb, var+"_err") for dscb in scanner.dscb_ls] if var not in ["integral"] else np.zeros_like(y_mean)

        if var in ["mean"]:
            var_latex = r"#mu"
            y_min = 124.8 # 124.7
            y_max = 124.96
            y_units = "GeV"
        elif var in ["sigma"]:
            var_latex = r"#sigma"
            y_min = 0.85  # min(y_mean) * 0.90
            y_max = 2.4   # max(y_mean) * 1.10
            y_units = "GeV"
        elif var in ["integral"]:
            var_latex = r"integral"
            y_min = min(y_mean) * 0.90
            y_max = max(y_mean) * 1.10
            y_units = ""
        else:
            raise ValueError

        relmin = scanner.relm4muerr_min
        relmax = scanner.relm4muerr_max
        # Prepare the plot decor.
        title_template = r"Variation of DSCB fit  %s over different fit ranges" % var_latex
        # title_template = r"#splitline{Variation of DSCB fit  ? over different fit ranges}"
        # title_template += r"{%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%}" % (relmin, relmax)
        x_label = r"lower bound of fit range"
        y_label = r"%s(DSCB)" % var_latex
        x_min, x_max = 102, 118
        x_units = "GeV"
        line_color=1
        line_width=2
        marker_color=1
        marker_style=20
        marker_size=1
        # Make the plot and draw it to the canvas.
        gr = make_TGraphErrors(x, y_mean, x_err=None, y_err=y_mean_err,
                        use_binwidth_xerrs=False,
                        title=title_template,
                        x_label=x_label, x_min=x_min, x_max=x_max, x_units=x_units,
                        y_label=y_label, y_min=y_min, y_max=y_max, y_units=y_units,
                        line_color=line_color, line_width=line_width,
                        marker_color=marker_color, marker_style=marker_style, marker_size=marker_size,
                        )
        draw_options = "AP" if make_new_page_after else "AP same"
        gr.Draw(draw_options)

        leg = make_TLegend(x_dim=(0,1), y_dim=(0,1),
                           screenshot_dim=(878,872), buffer_dim=(176,438,610,131))
        # leg = ROOT.TLegend(0.20, 0.70, 0.50, 0.85)

        leg_text = r"%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%" % (relmin, relmax)
        leg.SetTextSize(0.02)
        leg.AddEntry(gr, leg_text, "lpfe")
        # leg.SetLineWidth(3)
        leg.SetBorderSize(1)
        leg.Draw("same")

        if make_new_page_after:
            canv.Print(outfile_pdf)

    # def plot_mean_vs_GeVfitrange(self, scanner, canv, outfile_pdf, draw_same_canvas=False):
    #     """Draw plots of mean(DSCB) and sigma(DSCB) across various GeV fit ranges.

    #     Example:
    #         Fit 1 range: [105, 140] GeV
    #         Fit 2 range: [110, 135] GeV

    #     Parameters
    #     ----------
    #     scanner : DSCBFitScanner()

    #     """
    #     relmin = scanner.relm4muerr_min
    #     relmax = scanner.relm4muerr_max
    #     internal_name = "gr_{}relm4muerr{}".format(relmin, relmax)
    #     # Get the data from the fits.
    #     x = [dscb.m4mu_min for dscb in scanner.dscb_ls]
    #     y_mean = [dscb.mean for dscb in scanner.dscb_ls]
    #     y_mean_err = [dscb.mean_err for dscb in scanner.dscb_ls]
    #     y_sigma = [dscb.sigma for dscb in scanner.dscb_ls]
    #     y_sigma_err = [dscb.sigma_err for dscb in scanner.dscb_ls]
    #     # y_max_plot = 1.1 * max(y)
    #     # Prepare the plot decor.
    #     title_template = r"Variation of DSCB fit  ? over different fit ranges"
    #     # title_template = r"#splitline{Variation of DSCB fit  ? over different fit ranges}"
    #     # title_template += r"{%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%}" % (relmin, relmax)
    #     leg_text = r"%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%" % (relmin, relmax)
    #     x_label = r"lower bound of fit range"
    #     y_label = r"?(DSCB)"
    #     x_min, x_max = 100, 125
    #     x_units, y_units = "GeV", "GeV"
    #     line_color=1
    #     line_width=2
    #     marker_color=1
    #     marker_style=20
    #     marker_size=1

    #     leg = ROOT.TLegend()

    #     grerr_mean = make_TGraphErrors(internal_name, x, y_mean, x_err=None, y_err=y_mean_err,
    #                     use_binwidth_xerrs=False,
    #                     title=title_template.replace("?", r"#mu"),
    #                     x_label=x_label, x_min=x_min, x_max=x_max, x_units=x_units,
    #                     y_label=y_label.replace("?", r"#mu"), y_min=124.7, y_max=125, y_units=y_units,
    #                     line_color=line_color, line_width=line_width,
    #                     marker_color=marker_color, marker_style=marker_style, marker_size=marker_size,
    #                     )
    #     grerr_mean.Draw("AP")
    #     canv.Print(outfile_pdf)

    #     grerr_sigma = make_TGraphErrors(internal_name, x, y_sigma, x_err=None, y_err=y_sigma_err,
    #                     use_binwidth_xerrs=False,
    #                     title=title_template.replace("?", r"#sigma"),
    #                     x_label=x_label, x_min=x_min, x_max=x_max, x_units=x_units,
    #                     y_label=y_label.replace("?", r"#sigma"), y_min=0.7, y_max=1.2, y_units=y_units,
    #                     line_color=line_color, line_width=line_width,
    #                     marker_color=marker_color, marker_style=marker_style, marker_size=marker_size,
    #                     )
 
    #     draw_options = "AP same" if draw_same_canvas else "AP"
    #     grerr_sigma.Draw(draw_options)


    #     canv.Print(outfile_pdf)

if __name__ == "__main__":
    print("...Preparing your area.")
    outfile_pdf = make_fullpath_file(outfile_dir, outfile_prefix, m4mu_min, m4mu_max, relmasserr_min, relmasserr_max, step_GeV)
    outfile_pdf = outfile_pdf.replace(".txt", ".pdf")
    check_overwrite(outfile_pdf, overwrite=overwrite)
    tdrStyle = setTDRStyle(show_statsbox=True)
    tdrGrid(tdrStyle, gridOn=False)
    ROOT.gROOT.SetBatch(True)

    print("...Loading DSCB fitter code...")
    load_cpp_code("./RooMyPDF_DSCB_C.so", cpp_DSCB_code_path)
    print("...Opening root file...")
    f = ROOT.TFile(infile_path)
    tree = f.Get("passedEvents")

    print("Using step_size = {} GeV".format(step_GeV))
    # print("Using bins:\n", fit_range_edges)
    print("...Preparing canvas.")
    canv = ROOT.TCanvas()
    canv.Print(outfile_pdf + "[")

    print("Performing DSCB fits.")
    scanner = DSCBFitScanner(relmasserr_min, relmasserr_max)  # Loop over these.
    scanner.do_DSCB_fits_over_mass4mu_range(m4mu_min, m4mu_max, step_GeV,
                                                max_lower_edge, min_upper_edge,
                                                tree, canv, outfile_pdf, n_bins)
    plotter = DSCBFitPlotter()
    plotter.plot_X_vs_GeVfitrange("mean", scanner, canv, outfile_pdf, make_new_page_after=True)
    plotter.plot_X_vs_GeVfitrange("sigma", scanner, canv, outfile_pdf, make_new_page_after=True)
    plotter.plot_X_vs_GeVfitrange("integral", scanner, canv, outfile_pdf, make_new_page_after=True)
    # plotter.plot_X_vs_GeVfitrange("sigma", scanner, canv, outfile_pdf, make_new_page_after=True)
    canv.Print(outfile_pdf + "]")
    print("Done.")