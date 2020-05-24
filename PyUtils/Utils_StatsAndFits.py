import warnings
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning
from d0_Utils.d0_fns import centers_of_binning_array, get_subset_mask, print_header_message
from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict

from PyUtils.Utils_Plotting import make_stats_legend_for_gaus_fit
#--- Fitting Functions ---#
def linear_func(x, b, m):
    """
    Calculate the y-value for a given x-value along a straight line, 
    with b = y-intercept, and m = slope.     
    """
    return b + float(m) * x

def gaussian_func(x, coeff, mu, sigma):#, normalize=False):
    """
    Calculate the y-value for a given x-value along a Gaussian curve, 
    with mean = mu, and stdev = sigma. 
    
    For reference, here is ROOT's docstring for the Gaus function:
    
    Double_t TMath::Gaus(Double_t x, Double_t mean, Double_t sigma, Bool_t norm)
{
   // Calculate a gaussian function with mean and sigma.
   // If norm=kTRUE (default is kFALSE) the result is divided
   // by sqrt(2*Pi)*sigma.

   if (sigma == 0) return 1.e30;
   Double_t arg = (x-mean)/sigma;
   Double_t res = TMath::Exp(-0.5*arg*arg);
   if (!norm) return res;
   return res/(2.50662827463100024*sigma); //sqrt(2*Pi)=2.50662827463100024
}
    """
    rel_deviation = (x - mu) / float(sigma)
    expon = -1 / 2. * np.float_power(rel_deviation, 2)
#    if (normalize):
#        coeff = 1 / np.sqrt( 2 * np.pi) / sigma  
    
    return coeff * np.exp(expon)

def fit_with_line(x_vals, y_vals, 
                  guess_params=[0,1], 
                  param_bounds=([-np.inf,-np.inf], [np.inf,np.inf]),
                  absolute_sigma=True):
    """
    Fit a line to a set of data. Returns the best (intercept, slope) which fit the data.
    
    FIXME: NOT UP-TO-DATE
    
    Parameters
    ----------
    x_vals : array-like
        The x values of the data. 
    y_vals : array-like
        The y values of the data.
    guess_params : array-like
        Initial guess for the parameters of the fitted Gaussian.: [coeff, mean, sigma].
        Putting guess values can speed up fit time and can make fits converge where they would otherwise fail.
        
    Returns
    -------
    popt : 3-element array
        A 3-element array of the the optimized parameters (coefficient, mean, sigma) of the fitted Gaussian: [coeff_best, mu_best, sigma_best]
        From the scipy.optimize.curve_fit docstring:
            "Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized."
    pcov : 2d array
        The covariance of popt.
        From the scipy.optimize.curve_fit docstring:
            "To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov))"
    """
    popt, pcov = curve_fit(linear_func, x_vals, y_vals, p0=guess_params, bounds=param_bounds)
    
    # FIXME: For some strange reason, sigma can turn out to be negative...
    # Turns out SciPy already has an optional parameter to fix this!
    # ...still getting negative sigma values... so take abs().
#     popt[2] = np.abs(popt[2])
    
    popt_err = np.sqrt(np.diag(pcov))

    return popt, popt_err, pcov


def fit_with_gaussian(x_vals, y_vals, 
                      guess_params=[1,0,1], 
                      param_bounds=([0,-np.inf,-np.inf], [np.inf,np.inf,np.inf]),
                      absolute_sigma=True):
    """
    Fit a Gaussian curve to a set of data. Returns the best (mu, sigma) which fit the data.
    
    Parameters
    ----------
    x_vals : array-like
        The x values of the data. 
    y_vals : array-like
        The y values of the data.
    guess_params : array-like
        Initial guess for the parameters of the fitted Gaussian.: [coeff, mean, sigma].
        Putting guess values can speed up fit time and can make fits converge where they would otherwise fail.
        
    Returns
    -------
    popt : 3-element array
        A 3-element array of the the optimized parameters (coefficient, mean, sigma) of the fitted Gaussian: [coeff_best, mu_best, sigma_best]
        From the scipy.optimize.curve_fit docstring:
            "Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized."
    pcov : 2d array
        The covariance of popt.
        From the scipy.optimize.curve_fit docstring:
            "To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov))"
    """
    popt, pcov = curve_fit(gaussian_func, x_vals, y_vals, p0=guess_params, bounds=param_bounds)
    
    # FIXME: For some strange reason, sigma can turn out to be negative...
    # Turns out SciPy already has an optional parameter to fix this!
    # ...still getting negative sigma values... so take abs().
    popt[2] = np.abs(popt[2])
    
    popt_err = np.sqrt(np.diag(pcov))

    return popt, popt_err, pcov

def iterative_fit_gaus(iterations, bin_edges, bin_vals, 
                       param_guess=[1,0,1], 
                       param_bounds=([-np.inf,-np.inf,-np.inf], [np.inf,np.inf,np.inf]),
                       num_sigmas=2,
                       ax=None, draw_on_axes=True, verbose=True, skip_bad_fit=False):
    """
    Fit a Gaussian function to the core of a distribution. 
    Do this iteratively to improve the fit of the core. 
    Can draw the fit onto the plot and print debug info. 
    
    Parameters
    ----------
    iterations : int
        Number of fits to perform.
    bin_edges : array
        An array of the edges of the bins of the histogram being fit. 
    bin_vals : array 
        The y-values of each bin, where the first element contains the y-value of first bin, etc.
        Equivalently, the number of entries in each bin.
        Should have length = len(bin_edges) - 1.
    param_guess : list or array-like, optional
        A 3-element list of your best guess at the optimum Gaussian fit parameters.
        It helps the fit converge faster (or sometimes - at all!).
    param_bounds : 2-tuple of lists
        ([param1_min, param2_min, ...], [param1_max, param2_max, ...])
    num_sigmas : int
        The number of sigmas to look in both directions away from the mean.
        Helps determine the range of the next fit. 
    ax : axes obj
        The axes obj on which to draw the plot.
    draw_on_axes : bool 
        If True, draw a plot.
    verbose : bool
        If True, print debug info.
    skip_bad_fit : bool
        If True, then don't stop analysis when a fit error is thrown. 
        
    Returns 
    -------
    stats_dict : dict
        Dictionary of the fit statistics, including errors.
        Each value is a list of the iterative fit stats, in order of the fit performed
        (i.e. first element in list is first fit, last element is last fit).
            Note: len(list) == iterations
         Keys (str) : Val
            'coeff_ls' : best-fit coefficient of Gaussian fit
            'mean_ls' : best-fit mean (mu)
            'stdev_ls' : best-fit standard deviation (sigma)
            'coeff_err_ls' : error on best-fit coefficient
            'mean_err_ls' : error on best-fit mean
            'stdev_err_ls' : error on best-fit standard deviation
            'fit_converged_ls' : bool list, True means fit did not throw an error
            'error_ls' : type of error returned; default is None
    ax : axes object
        The original axes object, but now may have Gaussian fits drawn to it.
    """
    if (verbose):
        msg = "Performing {} iterative Gaussian fits".format(iterations)
        if (iterations == 1):
            msg = msg.replace("fits", "fit") 
        print(msg)

    bin_centers = centers_of_binning_array(bin_edges)

    stats_dict = {
        'coeff_ls' : [],
        'coeff_err_ls' : [],
        'mean_ls' : [],
        'mean_err_ls' : [],
        'stdev_ls' : [],
        'stdev_err_ls' : [],
        'fit_converged_ls' : [],
        'error_ls' : [],
    }
    
    count = 0
    fit_converged = True
    error = None
    while count < iterations:
        count += 1
        
        if (count == 1) or not (fit_converged):
            # First fit: use original histogram's mean and stdev to choose a fit range.
            this_coeff = param_guess[0]
            this_mean  = param_guess[1]
            this_stdev = param_guess[2]
        else:
            # Otherwise use the last fit's optimized parameters.
            this_coeff = popt[0]
            this_mean  = popt[1]
            this_stdev = popt[2]
        
        best_params_so_far = np.array([this_coeff, this_mean, this_stdev])
        bounds_min = param_bounds[0]
        bounds_max = param_bounds[1]

        # Note: Jake, if you try to implement this function below,
        # You may have to do something like np.abs(param) first. 
        # Just think about it!
#        def calc_min_max_param_bounds(param, perc_shift=200):
#            param_min = param * (1 - perc_shift / 100.)
#            param_max = param * (1 + perc_shift / 100.)
#            return param_min, param_max

        def check_within_bounds(arr, bound_min, bound_max):
            if any(arr <= bound_min):
                print("[WARNING] Optimized parameters reached minimum bound.")
            if any(arr >= bound_max):
                print("[WARNING] Optimized parameters reached maximum bound.")

#        bounds_min, bounds_max = calc_min_max_param_bounds(best_params_so_far)
        check_within_bounds(best_params_so_far, bounds_min, bounds_max)
        
        this_x_min = this_mean - float(num_sigmas) * this_stdev
        this_x_max = this_mean + float(num_sigmas) * this_stdev
        
        # Make new fit range.
        mask = get_subset_mask(bin_centers, x_min=this_x_min, x_max=this_x_max)
        new_bin_centers = bin_centers[mask]
        new_bin_vals = bin_vals[mask]
        
        if (verbose):
            print("Fit {}, scanning +-{} sigmas from start point: {:.4E}".format(count, num_sigmas, this_mean))
            print("  this_x_min: {}".format(this_x_min))
            print("  this_x_max: {}".format(this_x_max))
            print("  bounds_min: {}".format(bounds_min))
            print("  bounds_max: {}".format(bounds_max))
            print("  Parameters guessed: {}".format(best_params_so_far))

        # Try to fit, but may encounter errors.
        with warnings.catch_warnings():
            warnings.simplefilter("error", OptimizeWarning)
            try:
                popt, popt_err, pcov = fit_with_gaussian(new_bin_centers, new_bin_vals, 
                                                        guess_params=best_params_so_far,  # Use previous params as a guess.
                                                        param_bounds=(bounds_min, bounds_max),
                                                        absolute_sigma=True)
            except ValueError:
                # SciPy docs: 
                # if either ydata or xdata contain NaNs, or if incompatible options are used.
                error = "ValueError"
                msg = "[WARNING] {} encountered during fit. Continuing on...".format(error)
                print_header_message(msg, pad_char="!", n_center_pad_chars=10)
                if (skip_bad_fit):
                    fit_converged = False
                else:
                    raise ValueError
            except RuntimeError:
                # SciPy docs: 
                # if the least-squares minimization fails.
                error = "RuntimeError"
                msg = "[WARNING] {} encountered during fit. Continuing on...".format(error)
                print_header_message(msg, pad_char="!", n_center_pad_chars=10)
                if (skip_bad_fit):
                    fit_converged = False
                else:
                    raise RuntimeError
            except OptimizeWarning:
                # SciPy docs: 
                # if the least-squares minimization fails.
                error = "OptimizeWarning"
                msg = "[WARNING] {} encountered during fit. Continuing on...".format(error)
                print_header_message(msg, pad_char="!", n_center_pad_chars=10)
                if (skip_bad_fit):
                    fit_converged = False
                else:
                    raise OptimizeWarning

        if (verbose):
            print("Best fit Gaussian parameters:")
            print("    coeff = {:.3E}".format(popt[0]) )
            print("    mean  = {:.3E}".format(popt[1]) )
            print("    sigma = {:.3E}".format(popt[2]) + "\n")
        
        # Update stats dict.
        try:
            stats_dict['coeff_ls'].append(popt[0])
            stats_dict['mean_ls'].append(popt[1])
            stats_dict['stdev_ls'].append(popt[2])
            stats_dict['coeff_err_ls'].append(popt_err[0])
            stats_dict['mean_err_ls'].append(popt_err[1])
            stats_dict['stdev_err_ls'].append(popt_err[2])
            stats_dict['fit_converged_ls'].append(fit_converged)
            stats_dict['error_ls'].append(error)
        except UnboundLocalError:
            # Then probably "local variable 'popt' referenced before assignment"
            stats_dict['coeff_ls'].append(None)
            stats_dict['mean_ls'].append(None)
            stats_dict['stdev_ls'].append(None)
            stats_dict['coeff_err_ls'].append(None)
            stats_dict['mean_err_ls'].append(None)
            stats_dict['stdev_err_ls'].append(None)
            stats_dict['fit_converged_ls'].append(fit_converged)
            stats_dict['error_ls'].append(error)
        
        if (draw_on_axes):
            if (fit_converged):
                gaus_vals_x = np.linspace(new_bin_centers[0], new_bin_centers[-1], 500)
                gaus_vals_y = gaussian_func(gaus_vals_x, *popt)
                leg_label_fit = make_stats_legend_for_gaus_fit(popt, popt_err)
                leg_label_fit = leg_label_fit.replace("Fit", "Gaus Fit {}:".format(count))
            else:
                gaus_vals_x = 0
                gaus_vals_y = 0
                leg_label_fit = "Fit did not converge due to {}".format(error)
    
            if ax is None:
                f, ax = plt.subplots(figsize=(10,8))
            # Plot the Gaussian curve across its optimal x-range.
            ax.plot(gaus_vals_x, gaus_vals_y, color=color_dict[count], label=leg_label_fit, linestyle='-', marker="")
            ax.legend(loc="upper right")
        
    return stats_dict, ax


def iterative_fit_gaus_unbinned(num_iters, data,
                                bin_edges=[], 
                                bin_vals=[],
                                num_sigmas=2,
                                ax=None, draw_on_axes=True, verbose=True):
                                                                        
    """
    Fit a Gaussian pdf to the core of an unbinned distribution of data. 
    Do this iteratively to improve the fit of the core, 
    using a fit range = [prev_fit_mu - num_sigmas * prev_fit_sigma,
                         prev_fit_mu + num_sigmas * prev_fit_sigma]
    The starting fit range (Fit 1) is the entire data range. 
    
    Parameters
    ----------
    num_iters : int
        Number of iterative fits to perform.
    data : array
        The original data to be fit.
        During the iterations, the fit range over the data will change 
        depending on the values from the previous fit.
    bin_edges : array
        An array of the edges of the bins of the histogram being fit. 
        Only used for display purposes, since an unbinned fit cannot be shown 2-D. 
    bin_vals : array
        The y-values of each bin, where the first element contains the y-value of first bin, etc.
        Equivalently, the number of entries in each bin.
        Should have length = len(bin_edges) - 1.    
    num_sigmas : int
        The number of sigmas to look in both directions away from the mean.
        Helps determine the range of the next fit. 
    ax : axes obj
        The axes obj on which to draw the plot.
    draw_on_axes : bool 
        If True, draw a plot.
    verbose : bool
        If True, print debug info.
        
    Returns 
    -------
    stats_dict : dict
        Dictionary of the fit statistics, including errors.
        Each value is a list of the iterative fit stats, in order of the fit performed
        (i.e. first element in list is first fit, last element is last fit).
            Note: len(list) == iterations performed
         Keys (str) : Val
            'mean_ls' : best-fit mean (mu)
            'stdev_ls' : best-fit standard deviation (sigma)
            'mean_err_ls' : error on best-fit mean
            'stdev_err_ls' : error on best-fit standard deviation
    ax : axes object
        The original axes object, but now may have Gaussian fits drawn to it.
    """
    if (draw_on_axes):
        assert len(bin_edges) >= 2
        assert len(bin_vals) >= 1
        
    if (verbose):
        msg = "Performing {} iterative UNBINNED Gaussian fits".format(num_iters)
        if (num_iters == 1):
            msg = msg.replace("fits", "fit") 
        print(msg)
        print("original size of data:",len(data))

    stats_dict = {
        'mean_ls' : [],
        'mean_err_ls' : [],
        'stdev_ls' : [],
        'stdev_err_ls' : [],
    }
    
    count = 0
    while count < num_iters:
        count += 1
        
        if count == 1:
            x_min = data.min()
            x_max = data.max()
        elif count > 1:
            # Use the num_sigmas, previous fit sigma, and mean, to determine next fit range. 
            x_min = bestfit_mean - num_sigmas * bestfit_stdev
            x_max = bestfit_mean + num_sigmas * bestfit_stdev
        else: 
            raise ValueError("[ERROR] Invalid value for `count`.")
            
        # Do unbinned fit.
        trimmed_data = data[(x_min <= data) & (data <= x_max)]
        stats_ls = unbinned_gaus_fit_RooFit(trimmed_data)
        
        bestfit_mean = stats_ls[0]
        bestfit_mean_err = stats_ls[1]
        bestfit_stdev = stats_ls[2]
        bestfit_stdev_err = stats_ls[3]
        
        if (verbose):
            print("Fit {} has len(trimmed_data) = {}".format(count, len(trimmed_data)))
            if count == 1:
                print("Fit {}, scanning over whole data range (x_min={:.2E}, x_max={:.2E})".format(count, x_min, x_max))
            else:
                print("Fit {}, scanning +-{} sigmas from start point: {:.4E}".format(count, num_sigmas, bestfit_mean))
            print("  this_x_min: {}".format(x_min))
            print("  this_x_max: {}".format(x_max))
            print("Best fit Gaussian parameters:")
            print("   mean  = {:.5E}".format(bestfit_mean) )
            print("   sigma = {:.5E}".format(bestfit_stdev) + "\n")
        
        # Update stats dict.
        stats_dict['mean_ls'].append(bestfit_mean)
        stats_dict['stdev_ls'].append(bestfit_stdev)
        stats_dict['mean_err_ls'].append(bestfit_mean_err)
        stats_dict['stdev_err_ls'].append(bestfit_stdev_err)
        
        if (draw_on_axes):
            gaus_vals_x = np.linspace(x_min, x_max, 500)

            # Make new fit range.
            bin_centers = centers_of_binning_array(bin_edges)
            mask = get_subset_mask(bin_centers, x_min=x_min, x_max=x_max)
            new_bin_centers = bin_centers[mask]
            new_bin_vals = bin_vals[mask]

            def gaussian_func_fixed_mu_sigma(x, C):
                """
                This is a Gaussian fit with mean and sigma FIXED.
                The best-fit mean and sigma have already been obtained from an unbinned fit.
                Now we just need to find the best-fit coeff for plotting purposes. 

                Hooray it works!!!
                """
                return C * np.exp( -1/2. * ((x - bestfit_mean)/bestfit_stdev)**2 )
            
            # Coeff. of Gaussian usually has value close to y-max. 
            param_guess=[new_bin_vals.max()]
            
            bestfit_coeff, pcov = curve_fit(gaussian_func_fixed_mu_sigma, 
                                   new_bin_centers, new_bin_vals, 
                                   p0=param_guess, bounds=([0], [np.inf]))
            
            bestfit_coeff_err = np.sqrt(np.diag(pcov))
            
            # bestfit_coeff and bestfit_coeff_err are 1 element arrays.
            popt = [bestfit_coeff[0], bestfit_mean, bestfit_stdev]
            popt_err = [bestfit_coeff_err[0], bestfit_mean_err, bestfit_stdev_err]

            gaus_vals_y = gaussian_func(gaus_vals_x, *popt)
            
            leg_label_fit = make_stats_legend_for_gaus_fit(popt, popt_err)
            leg_label_fit = leg_label_fit.replace("Fit", "Unbinned Gaus Fit {}:".format(count))
            
            if ax is None:
                f, ax = plt.subplots(figsize=(10,8))
            # Plot the Gaussian curve across its optimal x-range.
            ax.plot(gaus_vals_x, gaus_vals_y, 
                    color=color_dict[count], label=leg_label_fit, linestyle='-', marker="")
            ax.legend(loc="upper right")
        
    return stats_dict, ax

def unbinned_gaus_fit_RooFit(data):
    import ROOT
    ROOT.RooMsgService.instance().setStreamStatus(1,False)
    """
    Do a Gaussian fit on unbinned data. 
    Returns the best-fit Gaus parameters [mu, mu_err, sigma, sigma_err].
    
    Notes:
      Since unbinned data don't have a "y-value" (i.e. a bin height),
      then there is no need for a scaling coeff factor in the Gaussian. 
    """
    x_min = data.min()
    x_max = data.max()
    x_avg = data.mean()
    x_std = data.std()
    
    # Make and fill a TTree with unbinned data.
    tree = ROOT.TTree("tree","tree")
    x = np.zeros(1,dtype=float)
    tree.Branch("x",x,'x/D')
    for i in range(len(data)):
        x[0] = data[i]
        tree.Fill()

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

#--------------------------------------#
#----- Error propagation formulae -----#
#--------------------------------------#
def prop_err_x_div_y(x, y, dx, dy):
    """
    Return the ratio of two numbers (r = x/y) and the 
    corresponding uncertainty (dr), depending on (x, y, dx, dy).
    by  from x (numerator) and y (denominator).

    The error propagation formula is:
        (dr)^2 = (dr/dx)^2 * (dx)^2 + (dr/dy)^2 * (dy)^2 + 2 * dr/dx * dr/dy * dx*dy
        but we will ignore the final cross-term (correlation factor).
            Newton says: 
            dr/dx = 1/y
            dr/dy = -x/(y^2)
        So:
            dr = sqrt( (dx/y)^2 + (-x/(y)^2 * dy)^2 )

    *** This function been verified by an online calculator.
    """
    r = x / float(y)
    dr = np.sqrt((dx / float(y))**2 + (x / float(y)**2 * dy)**2)
    return r, dr