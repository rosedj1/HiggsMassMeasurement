import ROOT
import ROOT as r
import numpy as np
from array import array

from Utils_Python.Plot_Styles_ROOT.tdrstyle_official import fixOverlay
from Utils_ROOT.ROOT_Plotting import Root_Hist_GetLastBinRightEdge
from Utils_ROOT.ROOT_fns import skip_black_yellow_fit_line_colors

from d0_Studies.d0_Utils.d0_dicts import color_dict_RooFit

def RooFit_gaus_fit_unbinned(data, verbose=False):
    """
    Perform an unbinned Gaussian fit to some data.
    Return a list of the best-fit parameters: [mean, mean_err, sigma, sigma_err].

    Parameters
    ----------
    data : list or array-like
        The data to be fit.
    
    Returns
    -------
    fit_stats_ls : 4-element list
        The best-fit parameters from the fit: [mean, mean_err, sigma, sigma_err] 
    """
    data = np.array(data)
    
    min_ = data.min()
    max_ = data.max()
    avg_ = data.mean()
    std_ = data.std()
    
    x = ROOT.RooRealVar("x","x", min_, max_)
    mean = ROOT.RooRealVar("mean","Mean of Gaussian", avg_, min_, max_)
    sigma = ROOT.RooRealVar("sigma","Width of Gaussian", std_, 0, 1000)
    gauss = ROOT.RooGaussian("gauss","gauss(x,mean,sigma)", x, mean, sigma)

    # Make RooDataSet.
    ptr = array('f', [0.])
    tree = ROOT.TTree("tree", "tree")
    tree.Branch("x", ptr, "x/F")
    for val in data:
        ptr[0] = val
        tree.Fill()
    rds = ROOT.RooDataSet("rds","dataset from tree", tree, ROOT.RooArgSet(x))
    
    # Modify fit range, if needed.
#     if fit_range is None:
#         fit_x_min = x_min
#         fit_x_max = x_max
#     else:
#         fit_x_min = fit_range[0]
#         fit_x_max = fit_range[1]
        
    # Find and fill the mean and sigma variables.
    if (verbose):
        result = gauss.fitTo(rds)
    else:
        result = gauss.fitTo(rds, ROOT.RooFit.PrintLevel(-1))
#                          ROOT.RooFit.Range(fit_x_min, fit_x_max), 
        
    fit_stats_ls = [mean.getVal(), mean.getError(), sigma.getVal(), sigma.getError()]

    return fit_stats_ls

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
                                   linecolor=r.kBlue, markercolor=r.kBlue):
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
#     ROOT.RooMsgService.instance().setStreamStatus(1,False)

    BW_MEAN_PDG = 91.19
    BW_SIGMA_PDG = 2.44
    
    x_min = x_lim[0]
    x_max = x_lim[1]
    
    x = r.RooRealVar("x", "Mass (GeV/c^{2})", x_min, x_max)
    h_Zboson = r.RooDataHist("h_Zboson", "h_Zboson", r.RooArgList(x), hist)
    
    # Prepare the fit model.
    ## BW
    BW_mean_DY = r.RooRealVar("BW_mean_DY", "BW_mean_DY", BW_MEAN_PDG)
    BW_sigma_DY = r.RooRealVar("BW_sigma_DY", "BW_sigma_DY", BW_SIGMA_PDG)
    BW_DY = r.RooBreitWigner("BW_DY", "BW_DY", x, BW_mean_DY, BW_sigma_DY)
    
    ## CB
    Mean = r.RooRealVar("Mean", "Mean", 0.0, -5.0, 5.0)
    Sigma = r.RooRealVar("Sigma", "Sigma", 1, 0.0, 10.0)
    CB_alpha = r.RooRealVar("CB_alpha", "CB_alpha", 1, 0., 10)
    CB_exp = r.RooRealVar("CB_exp", "CB_exp", 5, 0., 30)
    CB = r.RooCBShape("CB", "CB", x, Mean, Sigma, CB_alpha, CB_exp)
    
    ## expo
    tau = r.RooRealVar("tau", "tau", 0, -1, 1)
    bkg = r.RooExponential("bkg","bkg", x, tau)
    fsig = r.RooRealVar("fsig","signal fraction", 0.7, 0.5, 1.2)
    
    # Combine models.
    BWxCB = r.RooFFTConvPdf("BWxCB","BWxCB", x, BW_DY, CB)
    Final_DY = r.RooAddPdf("Final_DY","Final_DY", r.RooArgList(BWxCB, bkg), r.RooArgList(fsig))

    # Plot the data and the fit.
    xframe = x.frame(r.RooFit.Title("BW x CB + exp"))
    h_Zboson.plotOn(xframe, r.RooFit.MarkerColor(markercolor))
    # h_Zboson_corr.plotOn(xframe, r.RooFit.LineColor(r.kGreen+2))
    
    # Draw the BWxCB.
    if fit_range is None:
        Final_DY.fitTo(h_Zboson)
    else:
        Final_DY.fitTo(h_Zboson, r.RooFit.Range(fit_range[0], fit_range[1]))
    Final_DY.plotOn(xframe, r.RooFit.LineColor(linecolor),
                            r.RooFit.MarkerColor(markercolor),
                            r.RooFit.MarkerSize(0.05))
    # Draw the bkg.
    Final_DY.plotOn(xframe, r.RooFit.Components("bkg"), 
                            r.RooFit.LineColor(r.kBlue), 
                            r.RooFit.LineStyle(r.kDashed))
    if (show_params):
        Final_DY.paramOn(xframe, r.RooFit.Layout(*params_box))
    xframe.getAttText().SetTextSize(0.025)
    xframe.getAttText().SetTextColor(linecolor)
    # xframe.getAttText().SetTextColor(color_line_corr)

#     leg_text  = "#splitline{#mu = %.3f #pm %.3f}" % (mean.getVal(),  mean.getError())
#     leg_text += "{#sigma = %.3f #pm %.3f}" % (sigma.getVal(),  sigma.getError())
    
#     leg = r.TLegend(0.03, 0.80, 0.20, 0.9)
#     leg.AddEntry("gauss", leg_text, "lep")
#     leg.Draw("same")

    fit_stats_dict = get_BWxCBplusExp_fit_stats(Mean, Sigma, CB_alpha, CB_exp, tau, fsig)
    return xframe, fit_stats_dict

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
    
    x = ROOT.RooRealVar("x","The Independent Variable", x_min, x_max)
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