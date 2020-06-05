import ROOT

from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import fixOverlay
from Utils_ROOT.ROOT_Plotting import Root_Hist_GetLastBinRightEdge
from Utils_ROOT.ROOT_fns import skip_black_yellow_fit_line_colors

from d0_Studies.d0_Utils.d0_dicts import color_dict_RooFit

def RooFit_gaus_fit_hist_binned(hist, x_roofit, xframe, fit_range=None, count=1, 
                                draw_stats=False, line_color=None, marker_color=None):
    """
    Return a list of the best-fit parameters obtained from fitting
    a histogram (hist) with a Gaussian function: [mean, mean_err, sigma, sigma_err].
    Also draws the histogram to the canvas.
    ###can make the canvas print the fitted histogram to a PDF.

    Parameters
    ----------
    hist : ROOT.TH1
        Histogram to be fit with a Gaussian fit. 
    canv : ROOT.TCanvas
        Canvas on which histogram will be painted.
    fit_range : 2-element list, optional
        [fit_x_min, fit_x_max]
        If None, then full x-axis range will be used. 
    count : int, optional
        Used to position the fit stats labels on the plot. 
        Also used to color the fit lines. 
        Example: count=1 is kBlack, count=2 is kRed

    Returns
    -------
    fit_stats_ls : 4-element list
        The best-fit parameters from the fit: [mean, mean_err, sigma, sigma_err] 
    """
    x_min = hist.GetBinLowEdge(0)
    x_max = Root_Hist_GetLastBinRightEdge(hist)
    x_avg = hist.GetMean()
    x_std = hist.GetStdDev()
    
    # x = ROOT.RooRealVar("x","x", x_min, x_max)
    mean = ROOT.RooRealVar("mean","Mean of Gaussian", x_avg, x_min, x_max)
    sigma = ROOT.RooRealVar("sigma","Width of Gaussian", x_std, 0, 999999)
    gauss = ROOT.RooGaussian("gauss","gauss(x,mean,sigma)", x_roofit, mean, sigma)

    ls = ROOT.RooArgList(x_roofit)
    data = ROOT.RooDataHist("data", "data set with x", ls, hist)

    # xframe = x.frame(ROOT.RooFit.Title("Some title here?"))

    if count == 1:
        # New histogram. Draw to a clean canvas. 
        hist.Draw()
        data.plotOn(xframe, ROOT.RooLinkedList())
        if marker_color is not None:
            data.plotOn(xframe, ROOT.RooFit.MarkerColor(marker_color))
        else:
            data.plotOn(xframe, ROOT.RooFit.MarkerColor(1))
    # else:
        # Don't need to keep drawing to canvas...
    #     # Probably using an iterative Gaus fit. Draw to previous canvas. 
    #     hist.Draw("same")

    # Modify fit range, if needed.
    if fit_range is None:
        fit_x_min = x_min
        fit_x_max = x_max
    else:
        fit_x_min = fit_range[0]
        fit_x_max = fit_range[1]
        
    # Find and fill the mean and sigma variables.
    result = gauss.fitTo(data, 
                         ROOT.RooFit.Range(fit_x_min, fit_x_max), 
                         ROOT.RooFit.PrintLevel(-1)
                        )
    
    color = color_dict_RooFit[count]
        
    # Draw the fit. 
    if line_color is not None:
        gauss.plotOn(xframe, ROOT.RooFit.LineColor(line_color))
    else: 
        gauss.plotOn(xframe, ROOT.RooFit.LineColor(color))

    # Put the fit params on the plot.
    if (draw_stats):
        text_y_min = 0.95 - 0.11*(count-1)  # In the top left corner of TCanvas("c","c",600,600).
        gauss.paramOn(xframe, ROOT.RooFit.Layout(0.16, 0.44, text_y_min))
#     text_y_min = 0.92 - 0.11*(count-1)
#     gauss.paramOn(xframe, ROOT.RooFit.Layout(0.19, 0.4, text_y_min))
        xframe.getAttText().SetTextSize(0.020)
        xframe.getAttText().SetTextColor(color)
    
    xframe.Draw("same")
    
#     pavelabel_x_start = ( float(x_max) + float(x_min) ) * 0.65
#     title = ROOT.TPaveLabel( pavelabel_x_start, 300, x_max, 350, 'Come on dude!' )
#     title.SetTextFont( 50 )
#     title.Draw("same")

    # Make a legend.
#     leg_text  = "#splitline{#mu = %.5f #pm %.5f}" % (mean.getVal(),  mean.getError())
#     leg_text += "{#sigma = %.5f #pm %.5f}" % (sigma.getVal(),  sigma.getError())
#     leg = ROOT.TLegend(0.10, 0.75, 0.35, 0.9)
#     leg = ROOT.TLegend()
#     leg.AddEntry("hist", leg_text, "pel")
#     leg.Draw("same")
    
#     leg_text  = "#splitline{#mu = %.3f #pm %.3f}" % (mean.getVal(),  mean.getError())
#     leg_text += "{#sigma = %.3f #pm %.3f}" % (sigma.getVal(),  sigma.getError())

    fixOverlay()

    fit_stats_ls = [mean.getVal(), mean.getError(), sigma.getVal(), sigma.getError()]

    return fit_stats_ls

def RooFit_iterative_gaus_fit(data, x_roofit, xframe, iters, fit_range=None, num_sigmas=2, 
                              binned_data=True, draw_stats=False, 
                              only_draw_last=False, line_color=None, marker_color=None):
    """
    Performs an iterative Gaussian fit on a histogram (hist), 
    changing the fit_range as described below. 
    Returns a dictionary of the fit statistics for each iterative fit.

    Performs either binned or unbinned fits. 

    Procedure:
        Fit 1 fit range = [data.mean - 2*data.std,   data.mean + 2*data.std]
        Fit 2 fit range = [Fit1.mu   - 2*Fit1.sigma, Fit1.mu   + 2*Fit1.sigma]

    Parameters
    ----------
    data : ROOT.TH1 or numpy.ndarray
        Histogram or array to be fit with a Gaussian curve. 
        Type will be verified when this function is called. 
    canv : ROOT.TCanvas
        Canvas on which histogram will be painted.
    iters : int
        The number of fit iterations to perform.
    fit_range : 2-element list, optional
        [fit_x_min, fit_x_max]
        If None, then full x-axis range will be used. 
    num_sigmas : int
        The number of sigmas to go away from the mean in each direction. 
        Used in determining the fit range. 
    binned_data : bool
        If True, must accept a TH1 parameter as `data` and will 
        perform a binned fit. 
    only_draw_last : bool
        Only draw very last fit on canvas. 

    Returns
    -------
    fit_stats_dict : dict
        A dictionary of the iterated fit statistics. 
        Key : str
        Val : list of each iteration. The last element is the last fit. 
            "mean_ls"     : [mean_fit1, mean_fit2, ...]
            "mean_err_ls" : [mean_err_fit1, mean_err_fit2, ...]
            "std_ls"      : [sigma_fit1, sigma_fit2, ...]
            "std_err_ls"  : [sigma_err_fit1, sigma_err_fit2, ...]
    """
    if (binned_data):
        assert isinstance(data, ROOT.TH1)
    else:
        # Doing unbinned fit, which requires an array 
        # or other methods which are not yet implemented...
        import numpy
        assert isinstance(data, numpy.ndarray)

    mean_ls = []
    mean_err_ls = []
    std_ls = []
    std_err_ls = []
    
    count = 0
    while count < iters:
        count += 1 
        
        # Determine fit range:
        if count == 1:
            if (binned_data):
                x_min_ = data.GetMean() - num_sigmas * data.GetRMS()
                x_max_ = data.GetMean() + num_sigmas * data.GetRMS()
            else:
                x_min_ = data.mean() - num_sigmas * data.std()
                x_max_ = data.mean() + num_sigmas * data.std()
        else:
            x_min_ = mean_ls[-1] - num_sigmas * std_ls[-1]
            x_max_ = mean_ls[-1] + num_sigmas * std_ls[-1]
            
        fit_range = [x_min_, x_max_]
        
        color = count
        if (only_draw_last):
            # Only show last iterated fit and stats.
            # canv.Clear()
            color = 1  # Choose 1 so that stats box is in top right. 

        if (binned_data):
            fit_stats_ls = RooFit_gaus_fit_hist_binned(data, x_roofit, xframe, fit_range=fit_range, 
                                                        count=color, draw_stats=draw_stats,
                                                        line_color=line_color, marker_color=marker_color)
        else:
            # Do unbinned fit. 
            fit_stats_ls = RooFit_gaus_fit_unbinned_from_array(data, x_roofit, xframe, fit_range=fit_range, 
                                                        count=color, draw_stats=draw_stats,
                                                        line_color=line_color, marker_color=marker_color)
        
        mean_ls.append(fit_stats_ls[0])
        mean_err_ls.append(fit_stats_ls[1])
        std_ls.append(fit_stats_ls[2])
        std_err_ls.append(fit_stats_ls[3])
        
    fit_stats_dict = {
        "mean_ls" : mean_ls,
        "mean_err_ls" : mean_err_ls,
        "std_ls" : std_ls,
        "std_err_ls" : std_err_ls,
    }
        
    return fit_stats_dict 

# def RooFit_hist_DSCB_fit(hist, canv, fit_range=None, count=1):
    
#     x_min = hist.GetBinLowEdge(0)
#     x_max = Root_Hist_GetLastBinRightEdge(hist)
#     x_avg = hist.GetMean()
#     x_std = hist.GetStdDev()
    
#     x = ROOT.RooRealVar("x","x", x_min, x_max)
#     mean = ROOT.RooRealVar("mean","Mean of Gaussian", x_avg, x_min, x_max)
#     sigma = ROOT.RooRealVar("sigma","Width of Gaussian", x_std, 0, x_max)
#     gauss = ROOT.RooGaussian("gauss","gauss(x,mean,sigma)", x, mean, sigma)

    
#     ROOT.RooRealVar("Sigma_DSCB", "Sigma_DSCB", 0.005, -0.03, 0.03)
#     ROOT.RooRealVar("AlphaL", "AlphaL", 1, 0, 5)
#     ROOT.RooRealVar("ExpL", "ExpL", 1, 0, 10)
#     ROOT.RooRealVar("AlphaR", "AlphaR", 1, 0, 5)
#     ROOT.RooRealVar("ExpR", "ExpR", 1, 0, 10)
    
#     ls = ROOT.RooArgList(x)
#     data = ROOT.RooDataHist("data", "data set with x", ls, hist)

#     xframe = x.frame(ROOT.RooFit.Title("Some title here?"))
#     data.plotOn(xframe, ROOT.RooLinkedList())
#     data.plotOn(xframe, ROOT.RooFit.MarkerColor(2))

#     hist.Draw("same")
#     # Modify fit range, if needed.
#     if fit_range is None:
#         fit_x_min = x_min
#         fit_x_max = x_max
#     else:
#         fit_x_min = fit_range[0]
#         fit_x_max = fit_range[1]
        
#     result = gauss.fitTo(data, 
#                          ROOT.RooFit.Range(fit_x_min, fit_x_max), 
#                          ROOT.RooFit.PrintLevel(-1)
#                         )
    
#     color = skip_black_yellow_fit_line_colors(count)
        
#     # Draw the fit. 
#     gauss.plotOn(xframe, ROOT.RooFit.LineColor(color))
#     # Put the fit params on the plot.
# #     text_y_min = 0.95 - 0.11*(count-1)  # In the top left corner of TCanvas("c","c",600,600).
# #     gauss.paramOn(xframe, ROOT.RooFit.Layout(0.16, 0.4, text_y_min))
#     text_y_min = 0.92 - 0.11*(count-1)
#     gauss.paramOn(xframe, ROOT.RooFit.Layout(0.19, 0.4, text_y_min))
#     xframe.getAttText().SetTextSize(0.020)
#     xframe.getAttText().SetTextColor(color)
    
#     xframe.Draw("same")
    
# #     pavelabel_x_start = ( float(x_max) + float(x_min) ) * 0.65
# #     title = ROOT.TPaveLabel( pavelabel_x_start, 300, x_max, 350, 'Come on dude!' )
# #     title.SetTextFont( 50 )
# #     title.Draw("same")

#     # Make a legend.
# #     leg_text  = "#splitline{#mu = %.5f #pm %.5f}" % (mean.getVal(),  mean.getError())
# #     leg_text += "{#sigma = %.5f #pm %.5f}" % (sigma.getVal(),  sigma.getError())
# #     leg = ROOT.TLegend(0.10, 0.75, 0.35, 0.9)
# #     leg = ROOT.TLegend()
# #     leg.AddEntry("hist", leg_text, "pel")
# #     leg.Draw("same")
    
#     fixOverlay()
#     canv.Update()
#     canv.Draw()
    
#     return mean, sigma

def RooFit_gaus_fit_unbinned_from_TTree(tree, canv, fit_range=None, count=1):
    """
    FIXME: Not yet implemented or tested.
    Do a Gaussian fit on unbinned data. 
    Returns the best-fit Gaus parameters [mu, mu_err, sigma, sigma_err].
    
    Notes:
      Since unbinned data don't have a "y-value" (i.e. a bin height),
      then there is no need for a scaling coeff factor in the Gaussian. 
    """
    ROOT.RooMsgService.instance().setStreamStatus(1,False)

    x_min = data.min()
    x_max = data.max()
    x_avg = data.mean()
    x_std = data.std()
    
#     c = ROOT.TCanvas()
#     c.Draw()
    
    x = ROOT.RooRealVar("x","x", x_min, x_max)
    mean = ROOT.RooRealVar("mean","Mean of Gaussian", x_avg, x_min, x_max)
    sigma = ROOT.RooRealVar("sigma","Width of Gaussian", x_std, 0, 999999)
    gauss = ROOT.RooGaussian("gauss","gauss(x,mean,sigma)", x, mean, sigma)

    dataset = ROOT.RooDataSet("dataset", "dataset", tree, ROOT.RooArgSet(x))

#     xframe = x.frame()
#     dataset.plotOn(xframe, ROOT.RooLinkedList())
#     gauss.plotOn(xframe)
#     xframe.Draw("same")
    
#     leg_text  = "#splitline{#mu = %.3f #pm %.3f}" % (mean.getVal(),  mean.getError())
#     leg_text += "{#sigma = %.3f #pm %.3f}" % (sigma.getVal(),  sigma.getError())
    
#     leg = ROOT.TLegend(0.03, 0.80, 0.20, 0.9)
#     leg.AddEntry("gauss", leg_text, "lep")
#     leg.Draw("same")

    # Find and fill the mean and sigma variables.
    result = gauss.fitTo(dataset, ROOT.RooFit.PrintLevel(-1))
    
    fit_stats_ls = [mean.getVal(), mean.getError(),
                    sigma.getVal(), sigma.getError()]
    
    return fit_stats_ls