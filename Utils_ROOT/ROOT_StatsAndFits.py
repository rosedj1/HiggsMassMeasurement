import ROOT as rt
import numpy as np
from array import array
from pprint import pprint

from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import fixOverlay
from Utils_Python.Utils_StatsAndFits import prop_err_on_dsigoversig
from Utils_Python.printing import print_header_message
from Utils_ROOT.ROOT_Plotting import Root_Hist_GetLastBinRightEdge, make_new_xframe
from Utils_ROOT.ROOT_fns import skip_black_yellow_fit_line_colors
from Utils_ROOT.ROOT_classes import make_TGraphErrors, make_TLegend
from d0_Studies.d0_Utils.d0_dicts import color_dict_RooFit

def RooFit_gaus_fit(data, binned_fit=True, fit_range=None, xframe=None, 
                    count=1, 
                    x_label="Independent Var", units="", x_lim=None,
                    verbose=False,
                    n_bins=100,
                    line_color=4, marker_color=1, 
                    force_line_color=False, view_plot=False):
    """
    Perform a binned or unbinned Gaussian fit to some data using RooFit.
    Return a 2-tuple: 
        (best_fit_stats_ls, xframe_with_gauss)
        where: best_fit_stats_ls: [mean, mean_err, sigma, sigma_err].

    TODO:
    [ ] Fix the tick marks which show up on top of stats box
        in top-right corner. Stats box should be topmost.
    [ ] A range error always pops up, even though all ranges seem to
        be taking effect.
        [#0] ERROR:Plotting -- Range 'fit_nll_gauss_roodata' not
             defined for variable 'x_var'. Ignoring ...

    Parameters
    ----------
    data : list (or array-like) or ROOT.TH1
        The data to be fit with a Gaussian function.
        NOTE:
        If the data is of type `list` or `numpy.ndarray`,
        then either a binned OR *unbinned fit* can be performed,
        depending on the value of `binned_fit`.
        However, if the data is of type `ROOT.TH1` then only 
        a *binned fit* is performed.
    binned_fit : bool, optional
        If True, force a binned fit.
        Otherwise: 
            if type(data) is `ROOT.TH1F`: binned fit.
            if type(data) is `list` or `numpy.ndarray`: unbinned fit.
    fit_range : 2-elem list, optional
        [x_min, x_max]
        Range along x-axis over which to do the fit.
        If None, fit range will be [min(data), max(data)].
    xframe : ROOT.RooRealVar.frame, optional
        Essentially a canvas onto which the data and Gaussian fit will be drawn.
        If a frame is not provided, a frame will be created.
    count : int, optional
        An iterator to help with setting the legend of best-fit params.
        Used to position the fit stats labels on the plot. 
        Also used to color the fit lines. 
        Example: count=1 is kBlack, count=2 is kRed
    x_label : str, optional
        Label of the x-axis. Can accommodate ROOT LaTeX by using '#'.
    units : str, optional
        Units of the x-axis variable.
    x_lim : 2-elem list, optional
        The min and max x-axis values to show on the plot: [x_min, x_max]
    verbose : bool
        If True, print juicy debug info.
    n_bins : int, optional
        The number of bins to use for a binned fit on an array.
        Used when type(data) is either `numpy.ndarray` or `list`
        AND the User specified to do a binned fit.
    line_color : ROOT.TColor (int), optional
    marker_color : ROOT.TColor (int), optional
    force_line_color : bool, optional
        If False, then line color is be determined by `line_color`.
        When count=5, the line color is changed to kOrange (kYellow is bad).
        However, if `force_line_color`, then line_color is necessarily
        determined by `line_color`
    view_plot : bool
        If True, then the canvas will be drawn to the screen.
        Used for testing purposes.

    Returns
    -------
    A 2-tuple: (fit_stats_ls, xframe)
        fit_stats_ls : 4-element list
            The best-fit parameters from the fit: [mean, mean_err, sigma, sigma_err]
        xframe : ROOT.RooFit.RooRealVar.frame
            Essentially a canvas with the data and Gaussian fit drawn on.
            Still need to do xframe.Draw() onto a TCanvas.
    """
    rt.RooMsgService.instance().setStreamStatus(1,False)
    # Investigate data.
    if isinstance(data, rt.TH1):
        msg = (
            f"[ERROR] Your data are in a histogram, but you are requesting an unbinned fit.\n"
            f"  Perhaps you should set `binned_fit = True`."
            )
        assert binned_fit, msg
        x_min = data.GetBinLowEdge(0)
        x_max = Root_Hist_GetLastBinRightEdge(data)
        x_avg = data.GetMean()
        x_std = data.GetStdDev()
        n_entries = data.GetEntries()
    elif (isinstance(data, list) or isinstance(data, np.ndarray)):
        data = np.array(data)
        x_min = data.min()
        x_max = data.max()
        x_avg = data.mean()
        x_std = data.std()
        n_entries = len(data)
    else:
        raise TypeError(f"Data (type: {type(data)}) must of type `list`, `numpy.ndarray`, or `ROOT.TH1`.")

    if binned_fit:
        print(f"[INFO] Performing a BINNED Gaussian fit.")
    else:
        print(f"[INFO] Performing an UNBINNED Gaussian fit.")
    # Make fit variables.
    x_var = rt.RooRealVar("x_var", x_label, x_min, x_max, units)  # (name, title, min, max, units)
    # x_var.setRange("test_range", x_min, x_max)
    mean = rt.RooRealVar("mean","Mean of Gaussian", x_avg, x_min, x_max)
    sigma = rt.RooRealVar("sigma","Width of Gaussian", abs(x_std), 0, abs(x_std)*10.0)
    gauss = rt.RooGaussian("gauss","The Gaussian Itself", x_var, mean, sigma)
    
    # Prepare appropriate RooFit data container.
    if isinstance(data, rt.TH1):
        # Must do binned fit. Make RooDataHist.
        tmp_hist = data
        # x_var.setBins(n_bins)
        ls = rt.RooArgList(x_var)
        roodata = rt.RooDataHist("roodata", "RooDataHist", ls, data)
    else:
        # Using array of data. Not sure yet if binned or unbinned fit.
        # Make a hist with data array for stats purposes.
        tmp_hist = rt.TH1F()
        tmp_hist.SetBins(n_bins, 
                         x_min if x_lim is None else x_lim[0], 
                         x_max if x_lim is None else x_lim[1]
                         )
        tmp_hist.StatOverflows(True)  # Count under/overflow bins in stats.
        tmp_hist.Sumw2()
        for val in data:
            tmp_hist.Fill(val)
        # Check type of binning.
        if binned_fit:
            # x_var.setBins(n_bins)
            ls = rt.RooArgList(x_var)
            roodata = rt.RooDataHist("roodata", "RooDataHist", ls, tmp_hist)
        else:
            # Do unbinned fit.
            ptr = array('f', [0.])
            tree = rt.TTree("tree", "tree")
            tree.Branch("x_var", ptr, "x_var/F")
            for val in data:
                ptr[0] = val
                tree.Fill()
            roodata = rt.RooDataSet("roodata","dataset from tree", tree, rt.RooArgSet(x_var))

    # Modify fit range, if needed.
    if fit_range is None:
        fit_x_min = x_min
        fit_x_max = x_max
    else:
        fit_x_min = fit_range[0]
        fit_x_max = fit_range[1]
        
    print(f"JAKE: find the fit range error1")
    # Find and fill the mean and sigma variables.
    if verbose:
        # result = gauss.fitTo(roodata, rt.RooFit.Range(fit_x_min, fit_x_max))
        result = gauss.fitTo(roodata, rt.RooFit.Range(fit_x_min, fit_x_max), rt.RooFit.PrintLevel(3),
        # rt.RooFit.NumCPU(4, 3) # (num_cpus, strategy)  # This just hangs...
        # rt.RooFit.Timer(True),
        # rt.RooFit.BatchMode(True), # Compute batch of likelihood values. Computations are faster.
        # rt.RooFit.Offset(True)
        )
        print(f"JAKE: find the fit range error2")
    else:
        result = gauss.fitTo(roodata, rt.RooFit.Range(fit_x_min, fit_x_max), rt.RooFit.PrintLevel(-1),
        # rt.RooFit.NumCPU(4, 3) # (num_cpus, strategy)  # This just hangs...
        # rt.RooFit.Timer(True),
        # rt.RooFit.BatchMode(True) # Compute batch of likelihood values. Computations are faster.
        )
        print(f"JAKE: find the fit range error3")
    print(f"JAKE: find the fit range error4")
    # Draw data and then the fit.
    # NOTE: After a MONSTROUS amount of testing,
    # I finally managed to figure out how to add all
    # objects to xframe and keep them persistent.
    # Be careful, if you must modify the below!
    if not view_plot: 
        rt.gROOT.SetBatch(True)
    # c = rt.TCanvas()
    # c.Draw("same")
    # If frame does not exist, make it.
    # x_var.setBins(n_bins)
    if xframe is None:
        xframe = x_var.frame(rt.RooFit.Range(x_min, x_max), rt.RooFit.Title(x_label))
    #---------------------#
    # FIXME: Needs testing:
    # if x_lim is not None:
    #     xframe.setRange("x_window", *x_lim)
    #---------------------#
    # tmp_hist.SetLineColor(rt.kBlue)
    # tmp_hist.SetStats(True)
    tmp_hist.SetLineWidth(0)  # Don't show hist.
    tmp_hist.Draw("sames")

    # st = tmp_hist.FindObject("stats")
    # xframe.addObject(st)
    # xframe.Draw("same")
    # xframe.SetStats(True)  # In RooFitStats box 

    # Put the Gaussian fit parameters on the plot.
    if count == 5 and not force_line_color:
        # Avoid yellow because it is difficult to read.
        line_color = rt.kOrange

    # gauss.paramOn(xframe, 
    #     rt.RooFit.Layout(0.13, 0.40, text_y_min)  # (x_min, x_max, y_max) as fraction of canvas width.
    #     )
    # xframe.getAttText().SetTextSize(0.020)
    # xframe.getAttText().SetTextColor(line_color)

    # roodata.plotOn(xframe, rt.RooLinkedList())

    # tmp_hist.Draw("same")
    xframe.addObject(tmp_hist, "sames", True)  # Draw on same canvas. True means hist will not be drawn.
    # st.Draw("same")
    # rt.gPad.Update()
    # rt.gPad.GetPrimitive("stats")
    roodata.plotOn(xframe, rt.RooFit.MarkerColor(marker_color))
    # roodata.statOn(xframe,  
    #                 # rt.RooFit.Label("Unbinned data:"),
    #                 rt.RooFit.Layout(0.68, 0.95, 0.95), 
    #                 rt.RooFit.Format("NMRE")
    #                 )
    # xframe.getAttText().SetTextSize(0.020)
    # xframe.getAttText().SetTextColor(1)

    # RooFit insists on plotting the very first Gaussian along the x-axis.
    # Therefore, draw it with linewidth=0 and then redraw a good line.
    gauss.plotOn(xframe, rt.RooFit.LineColor(line_color), rt.RooFit.Range(fit_x_min, fit_x_max), rt.RooFit.LineWidth(0))
    gauss.plotOn(xframe, rt.RooFit.LineColor(line_color), rt.RooFit.Range(fit_x_min, fit_x_max), rt.RooFit.LineWidth(2))
    # if count == 1:
    #     gauss.plotOn(xframe,   rt.RooFit.LineColor(line_color), rt.RooFit.Range(fit_x_min, fit_x_max), rt.RooFit.LineWidth(0))

    # Put fit stats on plot.
    text_y_min = 0.86 - 0.11*(count-1)
    pave = rt.TPaveText(0.17, text_y_min-0.11, 0.40, text_y_min, "NDC")
    pave.SetFillColor(0)
    pave.SetBorderSize(0) # Use 0 for no border.
    pave.SetTextAlign(12) # 22 is centered vert and horiz.
    pave.SetTextSize(0.016)
    pave.SetTextColor(line_color)
    pave.SetFillStyle(1001)  # Solid fill.
    pave.AddText(f"Fit {count}:")  # Accommodates LaTeX!
    pave.AddText(f"  #mu = {mean.getVal():.4g} #pm {mean.getError():.4g}")
    pave.AddText(f"  #sigma = {sigma.getVal():.4g} #pm {sigma.getError():.4g}")
    pave.AddText("  #chi^{2}/ndf = %.4f" % xframe.chiSquare())
    pave.Draw("same")

    # latex = rt.TLatex()
    # latex.SetNDC()
    # latex.SetTextSize(0.016)
    # latex.SetTextColor(line_color)

    # Ensures a stats box exists so that the hist
    # can access it down below.
    rt.gStyle.SetOptStat("iouRMe")

    # roodata.statOn(xframe)
    # tmp_hist.SetStats(True)

    # Add all objects to xframe and THEN draw it.
    # latex.DrawLatex(0.13, text_y_min, f"Fit {count}:")
    # print(f"length = {len(rt.gPad.GetListOfPrimitives())}")
    # print(f"last one: {rt.gPad.GetListOfPrimitives()[-1]}")
    # latex.DrawLatex(0.13, text_y_min-0.02, f"  #mu = {mean.getVal():.4g} #pm {mean.getError():.4g}")
    # latex.DrawLatex(0.13, text_y_min-0.04, f"  #sigma = {sigma.getVal():.4g} #pm {sigma.getError():.4g}")
    # latex.DrawLatex(0.13, text_y_min-0.06,  "  #chi^{2}/ndf = %.4g" % xframe.chiSquare())
    xframe.addObject(pave.Clone())  # Not sure if .Clone() is necessary, but it works!
    xframe.addObject(tmp_hist.FindObject("stats").Clone())
    xframe.Draw("same")
    xframe.findObject("stats").Draw()
    xframe.findObject("stats").SetX2NDC(0.95)
    xframe.findObject("stats").SetY2NDC(0.90)
    # rt.gPad.RedrawAxis()
    # rt.gPad.Modified()
    # rt.gPad.Update()

    # Make a legend.
    #     leg_text  = "#splitline{#mu = %.5f #pm %.5f}" % (mean.getVal(),  mean.getError())
    #     leg_text += "{#sigma = %.5f #pm %.5f}" % (sigma.getVal(),  sigma.getError())
    #     leg = rt.TLegend(0.10, 0.75, 0.35, 0.9)
    #     leg = rt.TLegend(0.03, 0.80, 0.20, 0.9)
    #     leg = rt.TLegend()
    #     leg.AddEntry("hist", leg_text, "pel")
    #     leg.Draw("same")

    #     leg_text  = "#splitline{#mu = %.3f #pm %.3f}" % (mean.getVal(),  mean.getError())
    #     leg_text += "{#sigma = %.3f #pm %.3f}" % (sigma.getVal(),  sigma.getError())

    # fixOverlay()

    # try:
    #     del x_var #tmp_hist
    # except NameError:
    #     pass
    fit_stats_ls = [mean.getVal(), mean.getError(), sigma.getVal(), sigma.getError(), xframe.chiSquare()]
    return (fit_stats_ls, xframe)

def RooFit_iterative_gaus_fit(data, binned_fit=False, switch_to_binned_fit=2000, iters=1, num_sigmas=2.5, 
                              n_bins=100, x_lim=None, fit_whole_range_first_iter=False,
                              xframe=None, x_label="Independent var", title="", units="", marker_color=1,
                              force_last_line_color=None, only_draw_last=False, verbose=False, view_plot=False,
                              use_data_in_xlim=False, use_smart_window=False):
    """Return a 2-tuple: 
        (dict of the fit statistics from iterative Gaussian fits on given data, 
         xframe on which the fits can be drawn)
    
    FIXME:
        - Make option to show hist stats on xframe.
    NOTE: 
        - The fit can be *binned* or *unbinned* depending on the bool `binned_fit`.
          A histogram can only do binned fits, but an array can be binned or unbinned.
        - The fit can also switch from unbinned to binned fit, if you need to save time for instance.
          if len(data) > `switch_to_binned_fit`, then a binned fit will be performed.
        - The fit range = [mean - num_sigmas*std, mean + num_sigmas*std], 
            and changes for each subsequent fit, since the mean and std change with each fit.
        for each iterative fit.

    Procedure (suppose num_sigmas=2):
        Fit 1 fit range = [data.mean - 2*data.std,   data.mean + 2*data.std]
        Fit 2 fit range = [Fit1.mu   - 2*Fit1.sigma, Fit1.mu   + 2*Fit1.sigma]
        etc.

    Parameters
    ----------
    data : ROOT.TH1 or list or numpy.ndarray
        Histogram or array-like to be fit with a Gaussian curve.
    binned_fit : bool, optional
        If True, force a binned fit. otherwise an unbinned fit will be performed.
    xframe : ROOT.RooRealVar.frame, optional
        Essentially a canvas onto which the data and Gaussian fits will be drawn.
        If a frame is not provided, a frame will be created.
    iters : int, optional
        The number of fit iterations to perform.
    num_sigmas : int, optional
        The number of sigmas to go away from the mean in each direction. 
        Used in determining the fit range. 
    n_bins : int, optional
        The number of bins along the x-axis.
    x_lim : 2-elem list, optional
        The min and max x-axis values to show on the plot: [x_min, x_max]
        If None, then min(data) and max(data) are used as plot limits.
    fit_whole_range_first_iter : bool, optional
        If True, the entire x-axis range will be fit over for first iteration.
        Otherwise, first fit range is: mean(data) +- num_sigmas * rms(data)
        NOTE:
        - When the std(distribution) >> sigma(core Gaussian),
        then mean(dist) +- N*std(dist) is HUGE.
        This has a tendency to give inconsistent fits.
        You should set this param to True.
        Otherwise, you will probably get consistent fits when this is False.
    x_label : str, optional
        Label of the x-axis. Can accommodate ROOT LaTeX by using '#'.
    units : str, optional
        Units of the x-axis.
    title : str, optional
    only_draw_last : bool, optional
        If True, only draw the very last fit on the frame. 
    switch_to_binned_fit : int, optional
        The max number of entries in an array on which an UNBINNED fit should
        be performed. If n_entries in array > switch_to_binned_fit, then a
        binned fit will be done.
    line_color : ROOT.TColor (int), optional
    marker_color : ROOT.TColor (int), optional
    force_last_line_color : ROOT.TColor (int), optional
        Last fit will have this line color.
        Useful for controlling the line color when doing `only_draw_last`.
    verbose : bool
        If True, print juicy debug info.
    view_plot : bool
        If True, then the canvas will be drawn to the screen
        and the terminal will hang, waiting for user input.
        Used for testing purposes.
    use_data_in_xlim : bool
        If True, only values between x_lim[0] and x_lim[1] are used to begin
        the fits.
    use_smart_window : bool
        If True, then set left edge of window at the first x val whose bin
        height > 5% of the max bin height. Right edge is set to last such bin.
        Still works with `use_data_in_xlim`, which just trims data range.
        Overrides any values passed to x_lim.
        Overrides `fit_whole_range_first_iter`.

    Returns
    -------
    A 2-tuple: (fit_stats_dict, xframe)
        fit_stats_dict : dict
            A dictionary of the iterated fit statistics. 
            Key : str
            Val : list of each iteration. The last element is the last fit. 
            {
                "mean_ls"      : [mean_fit1, mean_fit2, ...],
                "mean_err_ls"  : [mean_err_fit1, mean_err_fit2, ...],
                "std_ls"       : [sigma_fit1, sigma_fit2, ...],
                "std_err_ls"   : [sigma_err_fit1, sigma_err_fit2, ...],
                "fit_range_ls" : [(x_min_fit1, x_max_fit1), (x_min_fit2, x_max_fit2), ...]
            }
        xframe : ROOT.RooFit.RooRealVar.frame
            Essentially a canvas with the data and Gaussian fit drawn on.
            Still need to do xframe.Draw() onto a TCanvas.
    """
    if len(data) == 0:
        fit_stats_dict = {
            "mean_ls"      : [0] * iters,
            "mean_err_ls"  : [0] * iters,
            "std_ls"       : [0] * iters,
            "std_err_ls"   : [0] * iters,
            "chi2_ls"      : [0] * iters,
            "fit_range_ls" : [0] * iters,
        }
        if xframe is None:
            xframe = make_new_xframe(0, 1, x_lim, x_label, units, n_bins, title)[3]
        return (fit_stats_dict, xframe)
    if isinstance(data, rt.TH1):
        data_mean = data.GetMean()
        data_rms = data.GetRMS()
        data_x_min = data.GetBinLowEdge(1)
        data_x_max = Root_Hist_GetLastBinRightEdge(data)
    elif isinstance(data, np.ndarray) or isinstance(data, list):
        data = np.array(data)
        data_mean = data.mean()
        data_rms = data.std()
        if use_data_in_xlim:
            print(f"Truncating data to fit between: [{x_lim[0]}, {x_lim[1]}]")
            data = data[(x_lim[0] < data) & (data < x_lim[1])]
        data_x_min = data.min()
        data_x_max = data.max()
    else:
        msg = f"The type of data ({type(data)}) must be `list`, `numpy.ndarray`, or `ROOT.TH1`"
        raise TypeError(msg)

    if use_smart_window:
        print("...Using smart window for iter fit.")
        if x_lim is not None:
            print_header_message(
                f"[WARNING] Overriding specified x_lim: {x_lim}."
                )
        # Use Suzanne's clever trick of setting x_lim based on BIN HEIGHT!
        binedges = np.linspace(data_x_min, data_x_max, n_bins + 1)
        print(f"binedges:\n{binedges}")
        binentries, _ = np.histogram(data, bins=binedges)
        print(f"binentries:\n{binentries}")
        # Set x window to start at first bin whose height > 5% of max.
        mask = binentries > (0.05 * max(binentries))
        print(f"mask[:50]:\n{mask[:50]}")
        # Mask is shorter than binedges. Look at left/right cases separately.
        _xmin = min(binedges[:-1][mask]) # Left edge of first good bin.
        _xmax = max(binedges[1:][mask])  # Right edge of last good bin.
        x_lim = (_xmin, _xmax)
        print(f"new x_lim: {x_lim}")

    if xframe is None:
        x_min, x_max, x_roofit, xframe = make_new_xframe(data_x_min, data_x_max, x_lim,
                                                         x_label, units, n_bins, title)

    # Prepare final lists of fit stats.
    mean_ls = []
    mean_err_ls = []
    std_ls = []
    std_err_ls = []
    chi2_ls = []
    fit_range_ls = []
    
    if verbose:
        print(f"...Performing {iters} iterated Gaussian fits...")# Iterated Gaussian fits.
    count = 0
    while count < iters:
        count += 1 
        if count == 1:
            # Determine starting fit range:
            if use_smart_window:
                # Use full smart window.
                x_min_fit = x_lim[0]
                x_max_fit = x_lim[1]
            elif fit_whole_range_first_iter:
                # Use entire data range.
                x_min_fit = data_x_min
                x_max_fit = data_x_max
            else:
                # Use a range based on hist stats.
                x_min_fit = data_mean - num_sigmas * abs(data_rms)
                x_max_fit = data_mean + num_sigmas * abs(data_rms)
        else:
            # Make a fit range based on previous fit result.
            x_min_fit = mean_ls[-1] - num_sigmas * abs(std_ls[-1])
            x_max_fit = mean_ls[-1] + num_sigmas * abs(std_ls[-1])
        # If we pass data limits, then use min/max.
        if x_min_fit < data_x_min:
            x_min_fit = data_x_min
        if x_max_fit > data_x_max:
            x_max_fit = data_x_max
        fit_range = (x_min_fit, x_max_fit)
        # Set color of fit line.
        color = count

        if count == iters:
            # Last fit. Maybe change final line color. Maybe make a new frame.
            if force_last_line_color is not None:
                color = force_last_line_color
            if only_draw_last:
                # Only show last iterated fit and stats.
                print(f"Drawing only the last Iter. Gauss. Fit (#{count}) to xframe.")
                # Make a new xframe to clear the previous one(s). 
                x_min, x_max, x_roofit, xframe = make_new_xframe(data_x_min, data_x_max, x_lim,
                                                            x_label, units, n_bins, title)
        # Determine whether to use binned or unbinned fit.
        do_binned = binned_fit
        if not binned_fit:
            # Trying to do unbinned fit. 
            # See if it will take too long.
            if len(data) > switch_to_binned_fit:
                if verbose:
                    msg = (f"[INFO] Switching to BINNED fit since len(data) {len(data)} > the limit set {switch_to_binned_fit}:\n")
                    print(msg)
                do_binned = True


        fit_stats_ls, xframe = RooFit_gaus_fit(data, binned_fit=do_binned, fit_range=fit_range, xframe=xframe, 
                                               count=count,
                                               x_label=x_label, units="", x_lim=x_lim,
                                               verbose=verbose,
                                               n_bins=n_bins,
                                               line_color=color,
                                               marker_color=marker_color,
                                               force_line_color=only_draw_last,
                                               view_plot=view_plot)

        mean_ls.append(fit_stats_ls[0])
        mean_err_ls.append(fit_stats_ls[1])
        std_ls.append(fit_stats_ls[2])
        std_err_ls.append(fit_stats_ls[3])
        chi2_ls.append(fit_stats_ls[4])
        fit_range_ls.append(fit_range)
    # End while loop.
    if view_plot:
        input("Press enter to continue.")

    # Iterated fits are done. Collect stats.
    fit_stats_dict = {
        "mean_ls"      : mean_ls,
        "mean_err_ls"  : mean_err_ls,
        "std_ls"       : std_ls,
        "std_err_ls"   : std_err_ls,
        "chi2_ls"      : chi2_ls,
        "fit_range_ls" : fit_range_ls,
    }
    if verbose:
        print("Final iterated Gaussian fit statistics:")
        pprint(fit_stats_dict)
        print()
    return (fit_stats_dict, xframe)

def get_BWxCBplusExp_fit_stats(mean, sigma, alpha, n, tau, fsig):
    """
    Return a dict which stores the best-fit parameters and errors.
    
    Parameters
    ----------
    mean : float
        Mean of Gaussian core of Crystal Ball.
    sigma : float 
        Sigma of Gaussian core of Crystal Ball.
    alpha : float
        Describes where the Gaussian-to-left-tail-power-law switch takes place. 
        Gives the number of standard deviations when the switch happens.
    n : float
        The power of the power-law function.
        The greater n is, the more of a tail there will be.
    tau : float 
        The exponential decay rate. 
    fsig : float
        The normalizing coefficient in front of the BWxCB:
        fsig * BWxCB + [%] * bkg
    """
    fit_stats_dict = {
        "mean" : mean.getVal(), 
        "mean_err" : mean.getError(),
        "sigma" : sigma.getVal(),
        "sigma_err" : sigma.getError(),
        "alpha" : alpha.getVal(), 
        "alpha_err" : alpha.getError(),
        "n" : n.getVal(), 
        "n_err" : n.getError(),
        "tau" : tau.getVal(), 
        "tau_err" : tau.getError(),
        "fsig" : fsig.getVal(), 
        "fsig_err" : fsig.getError(),
    }
    return fit_stats_dict

def RooFit_CBxBWplusExp_fit_binned(hist, x_lim, fit_range=None, show_params=True, params_box=[0.2, 0.4, 0.85], 
                                   linecolor=rt.kBlue, markercolor=rt.kBlue):
    """
    Fit a histogram with a Breit-Wigner function convoluted with a Crystal Ball function
    while adding an exponential background function.
    
    NOTE: The BW parameters are fixed at the PDG values: 
        BW_MEAN_PDG = 91.19
        BW_SIGMA_PDG = 2.44

    Parameters
    ----------
    hist : ROOT.TH1
        The histogram already filled with data.
    x_lim : 2-element list
        The x-axis range to view the hist. [x_fit_min, x_fit_max]
    fit_range : 2-element list, optional
        [x_fit_min, x_fit_max]
    show_params : bool, optional
        If True, show the parameters in a box on the x_frame. 
    params_box : 3-element list, optional
        The coordinates of the fit parameters box.
        [x_left_side, x_right_side, y_top]  # as a fraction of canvas width.
    linecolor :
        
    Returns
    -------
    x_frame : ROOT.RooRealVar.frame
        The frame object which holds the plots. 
        Can be drawn to a TCanvas.
    """
    #     rt.RooMsgService.instance().setStreamStatus(1,False)

    BW_MEAN_PDG = 91.19
    BW_SIGMA_PDG = 2.44
    
    x_min = x_lim[0]
    x_max = x_lim[1]
    
    x = rt.RooRealVar("x", "Mass (GeV/c^{2})", x_min, x_max)
    h_Zboson = rt.RooDataHist("h_Zboson", "h_Zboson", rt.RooArgList(x), hist)
    
    # Prepare the fit model.
    ## BW
    BW_mean_DY = rt.RooRealVar("BW_mean_DY", "BW_mean_DY", BW_MEAN_PDG)
    BW_sigma_DY = rt.RooRealVar("BW_sigma_DY", "BW_sigma_DY", BW_SIGMA_PDG)
    BW_DY = rt.RooBreitWigner("BW_DY", "BW_DY", x, BW_mean_DY, BW_sigma_DY)
    
    ## CB
    Mean = rt.RooRealVar("Mean", "Mean", 0.0, -5.0, 5.0)
    Sigma = rt.RooRealVar("Sigma", "Sigma", 1, 0.0, 10.0)
    CB_alpha = rt.RooRealVar("CB_alpha", "CB_alpha", 1, 0., 10)
    CB_exp = rt.RooRealVar("CB_exp", "CB_exp", 5, 0., 30)
    CB = rt.RooCBShape("CB", "CB", x, Mean, Sigma, CB_alpha, CB_exp)
    
    ## expo
    tau = rt.RooRealVar("tau", "tau", 0, -1, 1)
    bkg = rt.RooExponential("bkg","bkg", x, tau)
    fsig = rt.RooRealVar("fsig","signal fraction", 0.7, 0.5, 1.2)
    
    # Combine models.
    BWxCB = rt.RooFFTConvPdf("BWxCB","BWxCB", x, BW_DY, CB)
    Final_DY = rt.RooAddPdf("Final_DY","Final_DY", rt.RooArgList(BWxCB, bkg), rt.RooArgList(fsig))

    # Plot the data and the fit.
    xframe = x.frame(rt.RooFit.Title("BW x CB + exp"))
    h_Zboson.plotOn(xframe, rt.RooFit.MarkerColor(markercolor))
    # h_Zboson_corr.plotOn(xframe, rt.RooFit.LineColor(rt.kGreen+2))
    
    # Draw the BWxCB.
    if fit_range is None:
        Final_DY.fitTo(h_Zboson)
    else:
        Final_DY.fitTo(h_Zboson, rt.RooFit.Range(fit_range[0], fit_range[1]))
    Final_DY.plotOn(xframe, rt.RooFit.LineColor(linecolor),
                            rt.RooFit.MarkerColor(markercolor),
                            rt.RooFit.MarkerSize(0.05))
    # Draw the bkg.
    Final_DY.plotOn(xframe, rt.RooFit.Components("bkg"), 
                            rt.RooFit.LineColor(rt.kBlue), 
                            rt.RooFit.LineStyle(rt.kDashed))
    if (show_params):
        Final_DY.paramOn(xframe, rt.RooFit.Layout(*params_box))
    xframe.getAttText().SetTextSize(0.025)
    xframe.getAttText().SetTextColor(linecolor)
    # xframe.getAttText().SetTextColor(color_line_corr)

    #     leg_text  = "#splitline{#mu = %.3f #pm %.3f}" % (mean.getVal(),  mean.getError())
    #     leg_text += "{#sigma = %.3f #pm %.3f}" % (sigma.getVal(),  sigma.getError())
        
    #     leg = rt.TLegend(0.03, 0.80, 0.20, 0.9)
    #     leg.AddEntry("gauss", leg_text, "lep")
    #     leg.Draw("same")

    fit_stats_dict = get_BWxCBplusExp_fit_stats(Mean, Sigma, CB_alpha, CB_exp, tau, fsig)
    return xframe, fit_stats_dict

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
        # Attributes below are useful when doing before/after pT corr.
        self.sigma_improv = None  # Fraction, not percent.
        self.sigma_improv_err = None  # Fraction, not percent
        self.mean_shift   = None  # MeV.
        self.mean_shift_err = None  # MeV

    def do_DSCB_fit(self, tree, canv, outfile_pdf,
        m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max, n_bins=100, method="AdHoc", year="2018",
        show_after_corr=False, make_new_page_after=True, mean_before_corr=0, sigma_before_corr=0, zoom=False):
        """Draw a DSCB fit to a canvas, store the info, and return the fit results.
        
        NOTE: Make sure you have an open TCanvas on which to draw.

        Parameters
        ----------
        tree : TTree
            Contains either set of branches: 
            - Set 1: mass4mu, mass4muErr, passedFullSelection, finalState
            - Set 2: m4mu, m4mu_corr
        canv : ROOT.TCanvas
            The canvas on which to draw the plots.
        outfile_pdf : str
            The full filepath to store the pdf.
        m4mu_min : float
            The starting point of your binning range.
        m4mu_max : float
            The ending point of your binning range.
        relm4muerr_min : float
            The min cut on mass4muErr.
        relm4muerr_max : float
            The max cut on mass4muErr.
        n_bins : int, default=100
            The number of bins in which to bin the mass4mu distribution.
            NOTE: It's actually an unbinned DSCB fit, but the binning
            is just for plotting purposes (can't see an unbinned fit!).
        method : str, default="AdHoc"
            Either "AdHoc" or "GeoFit".
        year : str, default="2018"
            The year of the sample.
        show_after_corr : bool, default=False
            If True, show a second DSCB fit on the same canvas.
            Controls positioning of labels in the C++ fit script.
        make_new_page_after : bool, default=True
            If True, then the canvas will save the page to the growing PDF
            and start a new page.
        mean_before_corr : float, default=0
            The mean of the previous DSCB fit which will be compared
            to a new fit.
            Useful when show_after_corr=True.
        sigma_before_corr : float, default=0
            The sigma of the previous DSCB fit which will be compared
            to a new fit.
            Useful when show_after_corr=True.
        zoom : bool, default=False
            If True, then x-axis range on plot will be same as fit range.
            Else, x-axis range will be [105, 145] GeV.

        Returns
        -------
        result : ROOT.RooFitResult
            A pointer which stores the results from a DSCB fit.
        """
        assert method in ["AdHoc", "GeoFit"]
        msg = "  * Performing DSCB fit: "
        if show_after_corr:
            msg = msg.replace(": ", "s - before and after pT corr:")
            color_line = rt.kRed
        else:
            color_line = rt.kBlue
        print(msg)
        self.info_printer(m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max, n_bins)
        # Perform selections and draw to canvas.
        integral = rt.vector('Double_t')()
        result = rt.fit_and_draw_DSCB(
                        tree, m4mu_min, m4mu_max, relm4muerr_min, relm4muerr_max,
                        integral, method, year, canv, outfile_pdf,
                        n_bins, color_line, show_after_corr, mean_before_corr, sigma_before_corr, zoom)
        if make_new_page_after:
            canv.Print(outfile_pdf)
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

    def parse_DSCBfit_result(self, res, name_mean="mean", name_sigma="sigma"):
        """
        Return a 4-tuple of the fit results from a DSCB fit.
        Store the fit result info.
        
        Parameters
        ----------
        res : RooFitResult* (pointer to a RooFitResult)
        name_mean : str
            The internal ROOT name of the RooRealVar mean.
        name_sigma : str
            The internal ROOT name of the RooRealVar sigma.
        """
        try:
            mean = res.floatParsFinal().find(name_mean).getVal()
            mean_err = res.floatParsFinal().find(name_mean).getError()
            sigma = res.floatParsFinal().find(name_sigma).getVal()
            sigma_err = res.floatParsFinal().find(name_sigma).getError()
        except AttributeError:
            mean = res.floatParsFinal().find("#mu").getVal()
            mean_err = res.floatParsFinal().find("#mu").getError()
            sigma = res.floatParsFinal().find("#sigma").getVal()
            sigma_err = res.floatParsFinal().find("#sigma").getError()
        self.mean = mean
        self.mean_err = mean_err
        self.sigma = sigma
        self.sigma_err = sigma_err
        assert all([x is not None for x in (mean, mean_err, sigma, sigma_err)])
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
        self.dscb_beforeaftercorr_ls = []

    def do_DSCB_fits_over_mass4mu_range(self, x_min, x_max, stepsize,
                                              max_lower_edge, min_upper_edge,
                                              tree, canv, outfile_pdf, n_bins,
                                              method, year,
                                              verbose=False,
                                              draw_beforeafter_corr=False,
                                              zoom=False):
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
        method : str
            Either "AdHoc" or "GeoFit".
        year : str
            The year of the sample.
        zoom : bool
            If True, then x-axis range on plot will be same as fit range.
            Else, x-axis range will be [105, 140] GeV.
        """
        # Make sure n_bins can be represented as an int.
        if not isinstance(n_bins, int):
            assert (n_bins).is_integer() 
        fit_range_edges = self.make_binedges_from_stepsize(x_min, x_max, stepsize)
        if verbose:
            print("This scanner made fit_range_edges:\n", fit_range_edges)
        for m4mu_min, m4mu_max in zip(fit_range_edges, fit_range_edges[::-1]):
            if (m4mu_min >= max_lower_edge) or (m4mu_max < min_upper_edge):
                break
            if draw_beforeafter_corr:
                # Drawing 2 DSCBs: one before pT corr and one after.
                dscb_before = DSCBFitter()
                dscb_after = DSCBFitter()
                # Plot m(4mu) before pT corr.
                res_before = dscb_before.do_DSCB_fit(tree, canv, outfile_pdf,
                          m4mu_min, m4mu_max, self.relm4muerr_min, self.relm4muerr_max, n_bins,
                          method=method, year=year,
                          show_after_corr=False, make_new_page_after=False, zoom=zoom)
                mean, mean_err, sigma, sigma_err = dscb_before.parse_DSCBfit_result(
                                                       res_before, name_mean="#mu", name_sigma="#sigma")
                # Plot m(4mu) after pT corr.
                # Also use the old sigma to calculate the improvement.
                # Also also, use the old mean to calculate the mass shift.
                res_after = dscb_after.do_DSCB_fit(tree, canv, outfile_pdf,
                          m4mu_min, m4mu_max, self.relm4muerr_min, self.relm4muerr_max, n_bins,
                          method=method, year=year,
                          show_after_corr=True, make_new_page_after=True, mean_before_corr=mean, sigma_before_corr=sigma, zoom=zoom)
                mean_corr, mean_corr_err, sigma_corr, sigma_corr_err = dscb_after.parse_DSCBfit_result(
                                                        res_after, name_mean="#mu^{corr.}", name_sigma="#sigma^{corr.}")
                dscb_after.sigma_improv = abs(sigma_corr - sigma) / float(sigma)
                dscb_after.sigma_improv_err = prop_err_on_dsigoversig(sigma, sigma_corr, sigma_err, sigma_corr_err)
                dscb_after.mean_shift = (mean_corr - mean) * 1000.0  # MeV.
                dscb_after.mean_shift_err = np.sqrt(mean_err**2 + mean_corr_err**2)
                self.dscb_beforeaftercorr_ls.append((dscb_before, dscb_after))
            else:
                # Draw one instance of a DSCB.
                dscb = DSCBFitter()
                res = dscb.do_DSCB_fit(tree, canv, outfile_pdf,
                          m4mu_min, m4mu_max, self.relm4muerr_min, self.relm4muerr_max, n_bins,
                          method=method, year=year,
                          show_after_corr=False, make_new_page_after=True, zoom=zoom)
                mean, mean_err, sigma, sigma_err = dscb.parse_DSCBfit_result(res)
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
        if not isinstance(bins, int):
            try:
                assert bins.is_integer()  # If float is close to int, then good.
                bins = int(bins)
            except AssertionError:
                msg = "  Cannot create int number of bins using:\n"
                msg += "  x_min={}, x_max={}, stepsize={}".format( x_min, x_max, stepsize)
                if allow_rebinning:
                    print("...Warning! You will not get the stepsize you requested")
                    print(msg)
                    bins = int(round(bins))
                else:
                    raise ValueError(msg)
        return bins
    
    def make_binedges_from_stepsize(self, x_min, x_max, stepsize):
        """Return an array of bin edges: [x_min, x_max, stepsize]."""
        bins = self.calc_num_bins(x_min, x_max, stepsize, allow_rebinning=True)
        # if not isinstance(bins, int):
        #     assert (bins).is_integer()  # If float is close to int, then good.
        return np.linspace(x_min, x_max, bins+1)

class DSCBFitPlotter:
    """Class to show how DSCB fit parameters change over different fit ranges."""
    
    def __init__(self):
        pass
        # self.relmass4mu_ls = relmass4mu_ls

    def plot_X_vs_GeVfitrange(self, var, scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=False):
        """Draw plots of var vs. different GeV fit ranges,
        where var is something like: mean(DSCB) and sigma(DSCB).

        Example:
            Fit 1 range: [105, 140] GeV
            Fit 2 range: [110, 135] GeV

        Parameters
        ----------
        var : str
            Attribute of DSCBFitter() object to be plotted as the y-coord.
            Supports: "mean", "sigma", "integral", "sigma_improv", "mean_shift"
        scanner : DSCBFitScanner obj
            A collection of DSCBFitter objects that fit over various ranges.
        corrected_dscb : bool
            If True, assumes you are using a DSCB with corrected pTs.

        Returns
        -------

        """
        # Get the data from the fits.
        dscb_tup_ls = scanner.dscb_beforeaftercorr_ls
        if len(dscb_tup_ls) == 0:
            # Not doing before and after pT corrections.
            assert len(np.shape(scanner.dscb_ls)) == 1
            x = [dscb.m4mu_min for dscb in scanner.dscb_ls]
            y_vals = [getattr(dscb, var) for dscb in scanner.dscb_ls]  # E.g., dscb.mean, dscb.sigma.
            y_vals_err = [getattr(dscb, var+"_err") for dscb in scanner.dscb_ls] if var not in ["integral"] else np.zeros_like(y_vals)
        else:
            # Before and after pT corrections! Two DSCBs per fit range.
            assert len(np.shape(dscb_tup_ls)) == 2
            x = [dscb_tup[0].m4mu_min for dscb_tup in dscb_tup_ls]
            which_dscb = 1 if corrected_dscb else 0
            # Set y vals.
            y_vals = [getattr(dscb_tup[which_dscb], var) for dscb_tup in dscb_tup_ls]  # E.g., mean, sigma, integral, sigma_improve, mean_shift.
            y_vals_err = [getattr(dscb_tup[which_dscb], var+"_err") for dscb_tup in dscb_tup_ls] if var not in ["integral"] else np.zeros_like(y_vals)
            if var in ["sigma_improv"]:
                # Convert to percent.
                y_vals = np.array(y_vals) * 100.0  # As %.
                y_vals_err = y_vals * 0.01
                # y_vals_err = np.array(y_vals_err) * 100.0  # As %.

                
                # y_vals = [dscb_tup[1].sigma_improv * 100.0 for dscb_tup in dscb_tup_ls]
                # y_vals = [dscb_tup[1].sigma_improv_err * 100.0 for dscb_tup in dscb_tup_ls]

                # for dscb_tup in dscb_tup_ls:
                #     dscb_bef, dscb_aft = dscb_tup[0], dscb_tup[1]
                #     y_err = prop_err_on_dsigoversig(dscb_bef.sigma, dscb_aft.sigma,
                #                                     dscb_bef.sigma_err, dscb_aft.sigma_err)
                #     y_vals_err.append(y_err * 100.0)  # As %.
                #     y_vals.append(dscb_aft.sigma_improv * 100.0)  # As %.
            # else:
                # Either mean, sigma, integral, or mean_shift.
                # if corrected_dscb:
                    # Corrected DSCB.
                #     y_vals = [getattr(dscb_tup[1], var) for dscb_tup in dscb_tup_ls]  # E.g., dscb.mean, dscb.sigma, dscb.integral, dscb.mean_shift.
                #     y_vals_err = [getattr(dscb_tup[1], var+"_err") for dscb_tup in dscb_tup_ls] if var not in ["integral"] else np.zeros_like(y_vals)
                # else:
                #     # Uncorrected DSCB.
                #     y_vals = [getattr(dscb_tup[0], var) for dscb_tup in dscb_tup_ls]  # E.g., dscb.mean, dscb.sigma, dscb.integral, dscb.mean_shift.
                #     y_vals_err = [getattr(dscb_tup[0], var+"_err") for dscb_tup in dscb_tup_ls] if var not in ["integral"] else np.zeros_like(y_vals)

        if var in ["mean"]:
            var_latex = r"#mu"
            y_min = 124.8 # 124.7
            y_max = 124.96
            y_units = "GeV"
        elif var in ["sigma"]:
            var_latex = r"#sigma"
            y_min = 0.85  # 0.85, min(y_vals) * 0.80  
            y_max = 1.25   # 1.25, max(y_vals) * 1.20    
            y_units = "GeV"
        elif var in ["integral"]:
            var_latex = r"integral"
            y_min = min(y_vals) * 0.90  # 37000  # 
            y_max = max(y_vals) * 1.10  # 50000  # 
            y_units = ""
        elif var in ["sigma_improv"]:
            var_latex = r"#frac{#Delta#sigma}{#sigma}"
            y_min = 0  # 3.5 # min(y_vals) * 0.90  # , 6.2 
            y_max = 10 # 7.5  # max(y_vals) * 1.10  # 7.6
            y_units = r"%"
        elif var in ["mean_shift"]:
            var_latex = r"#Deltam_{4#mu}"
            y_min = min(y_vals) * 0.9  # 3.5 # min(y_vals) * 0.90  # , 6.2 
            y_max = max(y_vals) * 1.1 # 7.5  # max(y_vals) * 1.10  # 7.6
            y_units = r"MeV"
        else:
            raise ValueError
        print(f"var ({var}) vals:\n{y_vals}")
        print(f"var ({var}) vals err:\n{y_vals_err}")

        relmin = scanner.relm4muerr_min
        relmax = scanner.relm4muerr_max
        # Prepare the plot decor.
        title_template = r"Variation of DSCB fit %s over different fit ranges" % var_latex
        # title_template = r"#splitline{Variation of DSCB fit  ? over different fit ranges}"
        # title_template += r"{%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%}" % (relmin, relmax)
        x_label = r"lower bound of fit range"
        y_label = r"%s(DSCB)" % var_latex
        x_min = min(x) - 1
        x_max = max(x) + 1  # = 102, 118
        x_units = "GeV"
        line_color = rt.kRed if corrected_dscb else rt.kBlack
        line_width = 1
        marker_color = rt.kRed if corrected_dscb else rt.kBlack
        marker_style = 20
        marker_size = 0.5
        # Make the plot and draw it to the canvas.
        gr = make_TGraphErrors(x, y_vals, x_err=None, y_err=y_vals_err,
                        use_binwidth_xerrs=False,
                        title=title_template,
                        x_label=x_label, x_min=x_min, x_max=x_max, x_units=x_units,
                        y_label=y_label, y_min=y_min, y_max=y_max, y_units=y_units,
                        line_color=line_color, line_width=line_width,
                        marker_color=marker_color, marker_style=marker_style, marker_size=marker_size)
        gr.x_vals = x
        gr.y_vals = y_vals
        gr.y_vals_err = y_vals_err
        return gr

        #--- Code below is for drawing the plot ---#
        # draw_options = "AP" if make_new_page_after else "AP same"
        # gr.Draw(draw_options)
        # if corrected_dscb:
        #     # Add corrected DSCB to uncorrected using a TMultiGraph.
        #     leg = make_TLegend(x_dim=(0,1), y_dim=(0,1),
        #                     screenshot_dim=(878,872), buffer_dim=(176,438,610,131))
        #     # leg = ROOT.TLegend(0.20, 0.70, 0.50, 0.85)

        #     # leg_text = r"%.2f < #deltam_{4#mu}/m_{4#mu} < %.2f%%" % (relmin, relmax)
        #     leg_text = r"After p_{T} corr."
        #     leg.SetTextSize(0.02)
        #     leg.AddEntry(gr, leg_text, "lpfe")
        #     # leg.SetLineWidth(3)
        #     leg.SetBorderSize(1)
        #     leg.Draw("same")

        # if make_new_page_after:
        #     canv.Print(outfile_pdf)







    # def make_plot_mean_vs_GeVfitrange(self, scanner, canv, outfile_pdf, make_new_page_after=True, corrected_dscb=False):
    #     """Return a TGraphErrors of the mean vs. changing fit range."""
    #     # Get the data from the fits.
    #     dscb_tup_ls = scanner.dscb_beforeaftercorr_ls
    #     y_vals = [getattr(dscb_tup[1], "mean") for dscb_tup in dscb_tup_ls]  # E.g., dscb.mean, dscb.sigma.
    #     y_vals_err = [getattr(dscb_tup[1], var+"_err") for dscb_tup in dscb_tup_ls] if var not in ["integral", "sigma_improv"] else np.zeros_like(y_vals)

    #     var_latex = r"#mu"
    #     x_min, x_max = 102, 118
    #     y_min = 124.8 # 124.7
    #     y_max = 124.96
    #     y_units = "GeV"

    #     relmin = scanner.relm4muerr_min
    #     relmax = scanner.relm4muerr_max
    #     # Prepare the plot decor.
    #     title_template = r"Variation of DSCB fit %s over different fit ranges" % var_latex
    #     # title_template = r"#splitline{Variation of DSCB fit  ? over different fit ranges}"
    #     # title_template += r"{%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%}" % (relmin, relmax)
    #     x_label = r"lower bound of fit range"
    #     y_label = r"%s(DSCB)" % var_latex
    #     x_units = "GeV"
    #     line_color=1
    #     line_width=2
    #     marker_color = rt.kRed if corrected_dscb else rt.kBlack
    #     marker_style=20
    #     marker_size=1
    #     # Make the plot and draw it to the canvas.
    #     gr = make_TGraphErrors(x, y_vals, x_err=None, y_err=y_vals_err,
    #                     use_binwidth_xerrs=False,
    #                     title=title_template,
    #                     x_label=x_label, x_min=x_min, x_max=x_max, x_units=x_units,
    #                     y_label=y_label, y_min=y_min, y_max=y_max, y_units=y_units,
    #                     line_color=line_color, line_width=line_width,
    #                     marker_color=marker_color, marker_style=marker_style, marker_size=marker_size,
    #                     )
    #     draw_options = "AP" if make_new_page_after else "AP same"
    #     gr.Draw(draw_options)

# def RooFit_hist_DSCB_fit(hist, canv, fit_range=None, count=1):
    
#     x_min = hist.GetBinLowEdge(0)
#     x_max = Root_Hist_GetLastBinRightEdge(hist)
#     x_avg = hist.GetMean()
#     x_std = hist.GetStdDev()
    
#     x = rt.RooRealVar("x","x", x_min, x_max)
#     mean = rt.RooRealVar("mean","Mean of Gaussian", x_avg, x_min, x_max)
#     sigma = rt.RooRealVar("sigma","Width of Gaussian", x_std, 0, x_max)
#     gauss = rt.RooGaussian("gauss","gauss(x,mean,sigma)", x, mean, sigma)

    
#     rt.RooRealVar("Sigma_DSCB", "Sigma_DSCB", 0.005, -0.03, 0.03)
#     rt.RooRealVar("AlphaL", "AlphaL", 1, 0, 5)
#     rt.RooRealVar("ExpL", "ExpL", 1, 0, 10)
#     rt.RooRealVar("AlphaR", "AlphaR", 1, 0, 5)
#     rt.RooRealVar("ExpR", "ExpR", 1, 0, 10)
    
#     ls = rt.RooArgList(x)
#     data = rt.RooDataHist("data", "data set with x", ls, hist)

#     xframe = x.frame(rt.RooFit.Title("Some title here?"))
#     data.plotOn(xframe, rt.RooLinkedList())
#     data.plotOn(xframe, rt.RooFit.MarkerColor(2))

#     hist.Draw("same")
#     # Modify fit range, if needed.
#     if fit_range is None:
#         fit_x_min = x_min
#         fit_x_max = x_max
#     else:
#         fit_x_min = fit_range[0]
#         fit_x_max = fit_range[1]
        
#     result = gauss.fitTo(data, 
#                          rt.RooFit.Range(fit_x_min, fit_x_max), 
#                          rt.RooFit.PrintLevel(-1)
#                         )
    
#     color = skip_black_yellow_fit_line_colors(count)
        
#     # Draw the fit. 
#     gauss.plotOn(xframe, rt.RooFit.LineColor(color))
#     # Put the fit params on the plot.
# #     text_y_min = 0.95 - 0.11*(count-1)  # In the top left corner of TCanvas("c","c",600,600).
# #     gauss.paramOn(xframe, rt.RooFit.Layout(0.16, 0.4, text_y_min))
#     text_y_min = 0.92 - 0.11*(count-1)
#     gauss.paramOn(xframe, rt.RooFit.Layout(0.19, 0.4, text_y_min))
#     xframe.getAttText().SetTextSize(0.020)
#     xframe.getAttText().SetTextColor(color)
    
#     xframe.Draw("same")
    
# #     pavelabel_x_start = ( float(x_max) + float(x_min) ) * 0.65
# #     title = rt.TPaveLabel( pavelabel_x_start, 300, x_max, 350, 'Come on dude!' )
# #     title.SetTextFont( 50 )
# #     title.Draw("same")

#     # Make a legend.
# #     leg_text  = "#splitline{#mu = %.5f #pm %.5f}" % (mean.getVal(),  mean.getError())
# #     leg_text += "{#sigma = %.5f #pm %.5f}" % (sigma.getVal(),  sigma.getError())
# #     leg = rt.TLegend(0.10, 0.75, 0.35, 0.9)
# #     leg = rt.TLegend()
# #     leg.AddEntry("hist", leg_text, "pel")
# #     leg.Draw("same")
    
#     fixOverlay()
#     canv.Update()
#     canv.Draw()
    
#     return mean, sigma
