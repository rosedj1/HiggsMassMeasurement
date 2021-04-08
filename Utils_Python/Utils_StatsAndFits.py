import warnings
import numpy as np
from scipy.optimize import curve_fit, OptimizeWarning

# Local imports. 
from d0_Utils.d0_fns import centers_of_binning_array, get_subset_mask
from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict
# from Utils_ROOT.ROOT_StatsAndFits import RooFit_gaus_fit  # FIXME: Giving circular import error.
from Utils_Python.Utils_Plotting import make_stats_legend_for_gaus_fit
from Utils_Python.Utils_Physics import perc_diff
from Utils_Python.printing import print_header_message

#-----------------------------#
#----- Fitting Functions -----#
#-----------------------------#
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

def crystal_ball_func(x, coeff, alpha, n, mu, sigma, normalize=False):
    """
    FIXME: 
      [ ] DEPRECATED - NEEDS TO BE SYNCHED UP WITH crystal_ball_doublesided_func
      [ ] Doesn't work with arrays yet.

    Calculate the y-value for a given x-value along a Crystal Ball function. 
    A Crystal Ball function is a Gaussian Core combined with a left power-law tail.
    
    Parameters
    ----------
    x : float
        Independent variable.
    alpha : float
        Describes where the Gaussian-to-power-law switch takes place. 
    n : float
        The power of the power-law function.
        The greater n is, the more of a tail there will be.
    mu : float
        The mean of the Gaussian core.
    sigma : float
        The stdev of the Gaussian core. 
    
    For reference, here is ROOT's docstring for the Crystal Ball function:
    
    Double_t RooCBShape::evaluate() const {
       Double_t t = (m-m0)/sigma;
       if (alpha < 0) t = -t;

       Double_t absAlpha = fabs((Double_t)alpha);

       if (t >= -absAlpha) {
         return exp(-0.5*t*t);
       }
       else {
         Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
         Double_t b= n/absAlpha - absAlpha;

         return a/TMath::Power(b - t, n);
       }
     }
     https://root.cern.ch/doc/master/RooCBShape_8cxx_source.html
    """
    # Follow Wikipedia: https://en.wikipedia.org/wiki/Crystal_Ball_function
#     C = 
#     if (normalize):
#         coeff = 1 / np.sqrt( 2 * np.pi) / sigma  
#         return
    if (alpha < 0):
        raise ValueError("alpha should not be negative")
        
    dev = (x - mu) / sigma
    
    if dev > -1*alpha:
        return coeff * np.exp(-0.5 * dev * dev )
    
    else:
        absalpha = abs(alpha)
    
        A = (n / absalpha)**n * np.exp(-0.5 * absalpha * absalpha)
        B = (n / absalpha) - absalpha
        
        return coeff * A * (B - dev)**(-n)

def crystal_ball_doublesided_func(x_arr, coeff, alphaL, nL, alphaR, nR, mu, sigma):
    """
    Return an array of double-sided Crystal Ball (DSCB) y-vals for each x-val in an array.
    A DSCB is a Gaussian Core "stitched together" with a left power-law tail and right power-law tail.
    
    Parameters
    ----------
    x_arr : list or array-like
        Independent variable.
    coeff : float
        Scaling coefficient for the entire function. 
        Usually has a value around the peak of the DSCB.
    alphaL : float
        Describes where the Gaussian-to-LEFT-tail-power-law switch takes place. 
        Gives the number of standard deviations when the switch happens.
    nL : float
        The power of the power-law function.
        The greater n is, the more of a tail there will be.
    alphaR : float
        Describes where the Gaussian-to-RIGHT-tail-power-law switch takes place. 
        Gives the number of standard deviations when the switch happens.
    nR : float
        The power of the power-law function.
        The greater n is, the more of a tail there will be.
    mu : float
        The mean of the Gaussian core.
    sigma : float
        The stdev of the Gaussian core. 
    
    For reference, here is ROOT's docstring for the Crystal Ball function:
    
    Double_t RooCBShape::evaluate() const {
       Double_t t = (m-m0)/sigma;
       if (alpha < 0) t = -t;

       Double_t absAlpha = fabs((Double_t)alpha);

       if (t >= -absAlpha) {
         return exp(-0.5*t*t);
       }
       else {
         Double_t a =  TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
         Double_t b= n/absAlpha - absAlpha;

         return a/TMath::Power(b - t, n);
       }
     }
     https://root.cern.ch/doc/master/RooCBShape_8cxx_source.html
    """
    x_arr = np.array(x_arr)

    # Follow Wikipedia: https://en.wikipedia.org/wiki/Crystal_Ball_function
    if (alphaL < 0) or (alphaR < 0):
        raise ValueError("alphaL or alphaR should not be negative")
        
    # Make array of relative deviations.
    dev = (x_arr - mu) / sigma

    # Use a Gaussian, if: -alphaL <= dev <= alphaR.
    gaus_reg = (dev >= -1*alphaL) & (dev <= alphaR)  # Makes a bool array.
    gaus_val = coeff * np.exp(-0.5 * dev * dev )
    y_vals = np.where(gaus_reg, gaus_val, x_arr)
    
    # Use a left tail power-law, if: dev < -alphaL.
    left_reg = (dev < -1*alphaL)
    absalphaL = abs(alphaL)
    AL = (nL / absalphaL)**nL * np.exp(-0.5 * absalphaL * absalphaL)
    BL = (nL / absalphaL) - absalphaL
    left_val = coeff * AL * (BL - dev)**(-nL)
    y_vals = np.where(left_reg, left_val, y_vals)

    # Use a right tail power-law, if: dev > alphaR.
    right_reg = (dev > alphaR)
    absalphaR = abs(alphaR)
    AR = (nR / absalphaR)**nR * np.exp(-0.5 * absalphaR * absalphaR)
    BR = (nR / absalphaR) - absalphaR
    right_val = coeff * AR * (BR + dev)**(-nR)
    y_vals = np.where(right_reg, right_val, y_vals)
        
    # Make sure all of the x_arr values got transformed.
    # Possible bug if x_arr == 0 == y?
    assert all([a != b for a,b in zip(x_arr,y_vals)])
    return y_vals
        
def exp_gaus_exp_func(x_arr, coeff, kL, kR, mu, sigma):
    """
    Calculate the y-value for a given x-value along exp-Gaus-exp function. 
    This function is a Gaussian Core stitched together with a left exponential tail
    and a right exponential tail.
    
    A UF postdoc, Souvdik Das, along with Jaco Konigsberg invented this function.
    
    Parameters
    ----------
    x_arr : list or array-like, float
        Array of x-vals for which this function will be evaluated.
        The independent variable.
    kL : float
        Describes the exponential decay of the left tail.
    kR : float
        Describes the exponential decay of the right tail.
    mu : float
        The mean of the Gaussian core.
    sigma : float
        The stdev of the Gaussian core. 
        
    https://arxiv.org/pdf/1603.08591.pdf
    """
    x_arr = np.array(x_arr)

    if (kL < 0) or (kR < 0):
        raise ValueError("kL or kR should not be negative")
        
    # Make array of relative deviations.
    dev = (x_arr - mu) / sigma

    # Use a Gaussian, if: -kL < dev <= kR.
    gaus_reg = (-kL < dev) & (dev <= kR)  # Makes a bool array.
    gaus_val = coeff * np.exp(-0.5 * dev * dev )
    # gaus_val = np.exp(-0.5 * dev * dev)
    y_vals = np.where(gaus_reg, gaus_val, x_arr)
    
    # Use a left tail power-law, if: dev <= -kL.
    left_reg = (dev <= -kL)
    left_exp = 0.5 * kL**2 + kL * dev
    left_val = coeff * np.exp(left_exp)
    # left_val = np.exp(left_exp)
    y_vals = np.where(left_reg, left_val, y_vals)

    # Use a right tail power-law, if: dev < -alphaL.
    right_reg = (kR < dev)
    right_exp = 0.5 * kR**2 - kR * dev
    right_val = coeff * np.exp(right_exp)
    # right_val = np.exp(right_exp)
    y_vals = np.where(right_reg, right_val, y_vals)

    # Make sure all of the x_arr values got transformed.
    # Possible bug if x_arr == 0 == y?
    assert all([a != b for a,b in zip(x_arr,y_vals)])

    return y_vals

#----------------#
#----- FITS -----#
#----------------#
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

# tmp imports
from scipy.optimize import curve_fit

def fit_with_crystal_ball_doublesided(x_vals, y_vals, 
                      guess_params=[1,1,1,1,1,0,1], 
                      param_bounds=([0,0,0,0,0,0,0], 
#                           [-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf,-np.inf],
                                    [np.inf,np.inf,np.inf,np.inf,np.inf,np.inf,np.inf]),
                      absolute_sigma=True):
    """
    x, coeff, alphaL, nL, alphaR, nR, mu, sigma
    Fit a double-sided crystal ball curve to a set of data. 
    Returns the best (coeff, alphaL, nL, alphaR, nR, mu, sigma) which fit the data.
    
    Parameters
    ----------
    x_vals : list or array-like
        The x values of the data. 
    y_vals : list or array-like
        The y values of the data.
    guess_params : list or array-like
        Initial guess for the parameters of the fitted DSCB.: [coeff, alphaL, nL, alphaR, nR, mu, sigma].
        Putting guess values can speed up fit time and can make fits converge where they would otherwise fail.
    param_bounds : 2-tuple of lists
        First element of tuple is a list of the minimum-allowed values for each parameter.
        Second element of tuple is a list of the maximum-allowed values.
        Each list should have length == number of parameters
    absolute_sigma : bool
        I think this ensures that sigma > 0, but this has not been verified.
        Check the scipy docs instead.
        
    Returns
    -------
    popt : 7-element array
        The optimized parameters of the fitted DSCB.
        From the scipy.optimize.curve_fit docstring:
            "Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized."
    popt_err : 7-element array 
        The uncertainties on the optimized parameters.
    pcov : 2d array
        The covariance of popt.
        From the scipy.optimize.curve_fit docstring:
            "To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov))"
    """
    popt, pcov = curve_fit(crystal_ball_doublesided_func, x_vals, y_vals, p0=guess_params, bounds=param_bounds)
    
    # FIXME: For some strange reason, sigma can turn out to be negative...
    # Turns out SciPy already has an optional parameter to fix this!
    # ...still getting negative sigma values... so take abs().
    popt[6] = np.abs(popt[6])
    
    popt_err = np.sqrt(np.diag(pcov))
    
    return popt, popt_err, pcov

def fit_with_exp_gaus_exp(x_vals, y_vals, 
                      guess_params=[1,1,1,0,1], 
                      param_bounds=([0,0,0,0,0], 
                                    [np.inf,np.inf,np.inf,np.inf,np.inf]),
                      absolute_sigma=True):
    """
    Fit an exp-gaus-exp curve to a set of data. 
    Returns the best (kL, kR, mu, sigma) which fit the data.
    
    Parameters
    ----------
    x_vals : list or array-like
        The x values of the data. 
    y_vals : list or array-like
        The y values of the data.
    guess_params : list or array-like
        Initial guess for the parameters of the fitted ExpGausExp.: [kL, kR, mu, sigma].
        Putting guess values can speed up fit time and can make fits converge where they would otherwise fail.
    param_bounds : 2-tuple of lists
        First element of tuple is a list of the minimum-allowed values for each parameter.
        Second element of tuple is a list of the maximum-allowed values.
        Each list should have length == number of parameters
    absolute_sigma : bool
        I think this ensures that sigma > 0, but this has not been verified.
        Check the scipy docs instead.
        
    Returns
    -------
    popt : 4-element array
        The optimized parameters of the fitted ExpGausExp.
        From the scipy.optimize.curve_fit docstring:
            "Optimal values for the parameters so that the sum of the squared residuals of f(xdata, *popt) - ydata is minimized."
    popt_err : 4-element array 
        The uncertainties on the optimized parameters.
    pcov : 2d array
        The covariance of popt.
        From the scipy.optimize.curve_fit docstring:
            "To compute one standard deviation errors on the parameters use perr = np.sqrt(np.diag(pcov))"
    """
    popt, pcov = curve_fit(exp_gaus_exp_func, x_vals, y_vals, p0=guess_params, bounds=param_bounds)
    
    # FIXME: For some strange reason, sigma can turn out to be negative...
    # Turns out SciPy already has an optional parameter to fix this!
    # ...still getting negative sigma values... so take abs().
    popt[4] = np.abs(popt[4])
    
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
                                ax=None, draw_on_axes=True, framealpha=1.0, verbose=True):
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
            # x_min = data.min()
            # x_max = data.max()
            x_min = data.mean() - num_sigmas * abs(data.std())
            x_max = data.mean() + num_sigmas * abs(data.std())
        elif count > 1:
            # Use the num_sigmas, previous fit sigma, and mean, to determine next fit range. 
            x_min = bestfit_mean - num_sigmas * abs(bestfit_stdev)
            x_max = bestfit_mean + num_sigmas * abs(bestfit_stdev)
        else: 
            raise ValueError("[ERROR] Invalid value for `count`.")
            
        # Do unbinned fit.
        trimmed_data = data[(x_min <= data) & (data <= x_max)]
        stats_ls, xframe = RooFit_gaus_fit(trimmed_data, 
                                           binned_fit=True, fit_range=None, xframe=None, 
                                           count=1, 
                                           x_label="Independent Var", 
                                           draw=draw_on_axes, verbose=verbose,
                                           n_bins=100,
                                           line_color=4, marker_color=1)
        
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
            ax.legend(loc="upper right", framealpha=framealpha)
        
    return stats_dict, ax

#--------------------------------------#
#----- Error Propagation Formulae -----#
#--------------------------------------#
def prop_err_x_plus_y(x, y, dx, dy):
    """
    Return the sum of two numbers (z = x + y) and the 
    corresponding uncertainty (dz), depending on (x, y, dx, dy).

    The error propagation formula is:
        (dz)^2 = (dz/dx)^2 * (dx)^2 + (dz/dy)^2 * (dy)^2 + 2 * dz/dx * dz/dy * dx*dy
        but we will ignore the final cross-term (correlation factor).
            Newton says: 
            dz/dx = 1
            dz/dy = 1
        So:
            dz = sqrt( (dx)^2 + (dy)^2 )
    """
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    dx = np.array(dx, dtype=float)
    dy = np.array(dy, dtype=float)
    
    z = x + y
    dz = np.sqrt(dx**2 + dy**2)
    return z, dz

def prop_err_x_times_y(x, y, dx, dy):
    """
    Return a 2-tuple:
        tup[0]: product of two numbers (r = x * y)
        tup[1]: corresponding uncertainty (dr), depending on (x, y, dx, dy).

    The error propagation formula is:
        (dr)^2 = (dr/dx)^2 * (dx)^2 + (dr/dy)^2 * (dy)^2 + 2*dr/dx*dr/dy * dx*dy
        but we will ignore the final cross-term (correlation factor).
            Newton says: 
            dr/dx = y
            dr/dy = x
        So:
            dr = sqrt((y * dx)^2  + (x * dy)^2)

    *** This function been verified by an online calculator.
    """
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    dx = np.array(dx, dtype=float)
    dy = np.array(dy, dtype=float)
    
    r = x * y
    dr = np.sqrt((y * dx)**2  + (x * dy)**2)
    return r, dr

def prop_err_x_div_y(x, y, dx, dy):
    """
    Return a 2-tuple:
        tup[0]: ratio of two numbers (r = x/y)
        tup[1]: corresponding uncertainty (dr), depending on (x, y, dx, dy).

    The error propagation formula is:
        (dr)^2 = (dr/dx)^2 * (dx)^2 + (dr/dy)^2 * (dy)^2 + 2*dr/dx*dr/dy * dx*dy
        but we will ignore the final cross-term (correlation factor).
            Newton says: 
            dr/dx = 1/y
            dr/dy = -x/(y^2)
        So:
            dr = sqrt( (dx/y)^2 + (-x/(y)^2 * dy)^2 )

    *** This function been verified by an online calculator.
    """
    x = np.array([x])#, dtype=float)
    y = np.array([y])#, dtype=float)
    dx = np.array([x])#, dtype=float)
    dy = np.array([y])#, dtype=float)
    assert all(arr.size > 0 for arr in [x, y, dx, dy])
    
    undef = np.ones_like(x) * np.inf
    # Where denom is not 0, do the division. Elsewhere put inf.
    r = np.true_divide(x, y, out=undef, where=y!=0)
    # r = x / y
    dr = np.sqrt((dx / y)**2 + (x / y**2 * dy)**2)
    return r, dr

def prop_err_on_dsigoversig(sig1, sig2, sig_err1, sig_err2):
    """
    Returns the error on (sig2 - sig1) / sig1 (as a fraction, not percentage).

    If we let: 
        r = (n - b) / b = n/b - 1
        n +- dn
        b +- db
    Then the final propagation formula is:
        dr = n/b * sqrt[ (dn/n)^2 + (db/b)^2 ]

    Derivation:
        The error propagation formula is:
        (dr)^2 = (dr/dn)^2 * (dn)^2 + (dr/db)^2 * (db)^2 + 2*dr/dn*dr/db * dn*db
            Newton says: 
            dr/dn = 1/b
            dr/db = -n/(b^2)
        Plugging these in and ignoring the cross-term (correlation factor):
        (dr)^2 = (1/b)^2 * (dn)^2 + [-n/(b^2)]^2 * (db)^2
               = (1/b)^2 * (dn)^2 + 1/b^2*(n/b)^2 * (db)^2
               = (1/b)^2 * [(dn)^2*(n/n)^2 + (n/b)^2*(db)^2]
               = (n/b)^2 * [(dn/n)^2 + (db/b)^2]

    NOTE: This formula should be able to be derived from prop_err_x_div_y()
          above, but I haven't figured it out yet.

    Parameters
    ----------
    sig1 : float or array-like
        The initial value, used as a reference for sig2.
    sig2 : float or array-like
        The final value, which is compared to sig1.
    sig_err1 : float or array-like
        The error on sig1.
    sig_err2 : float or array-like
        The error on sig2.
    """
    sig1 = np.array([sig1], dtype=float) if not isinstance(sig1, np.ndarray) else np.array(sig1, dtype=float)
    sig2 = np.array([sig2], dtype=float) if not isinstance(sig2, np.ndarray) else np.array(sig2, dtype=float)
    sig_err1 = np.array([sig_err1], dtype=float) if not isinstance(sig_err1, np.ndarray) else np.array(sig_err1, dtype=float)
    sig_err2 = np.array([sig_err2], dtype=float) if not isinstance(sig_err2, np.ndarray) else np.array(sig_err2, dtype=float)
    
    relsig1 = sig_err1 / sig1
    relsig2 = sig_err2 / sig2
    return (sig2 / sig1) * np.sqrt(relsig1**2 + relsig2**2)

#-------------------------#
#----- Stats Helpers -----#
#-------------------------#
def get_standarderrorofmean(arr):
    """Return the standard error of the mean."""
    if not isinstance(arr, np.ndarray):
        arr = np.array(arr)
    return np.std(arr) / np.sqrt(len(arr))
    
def check_fit_convergence(stat_ls, max_perc_diff=5, compare_to_last=3, alarm_level="warning"):
    """Raise a warning or error if a list of subsequent fit values did not converge.
    
    NOTE: 
        The idea here is to see if the final converged fit value differs 
        too much from the average of the last few fit values.
        This is not the most sophisticated way to check, probably assuming some 
        kind of Gaussian distribution among the mean values and comparing the final
        fit value to this mean would be better.
        Another idea could be to see how far the final fit val differs from the 
        sigma of such a Gauss dist.

    Parameters
    ----------
    stat_ls : list
        An ORDERED list of fit values (e.g., iterative Gaussian fit means).
    max_perc_diff : float, optional
        The maximum allowed percent difference (%) between: 
            - the average of the last `compare_to_last` entries of `stat_ls`, and
            - the final fit value (last element in stat_ls)
    compare_to_last : int, optional
        Take the mean of this many elements from the back of `stat_ls`
        and compare that mean to the very last element to check for convergence.
    alarm_level : str, optional
        What level of complaining to print out, in case convergence seems to fail.
        Can choose: ["warning", "error"].
        Choosing `"error"` will raise a ValueError.
    """
    # Sanity checks.
    assert compare_to_last >= 1
    while len(stat_ls) < compare_to_last:
        compare_to_last -= 1
    alarm = alarm_level.upper()
    assert alarm in ["WARNING", "ERROR"]
    # Compare to the elements specified.
    arr = np.array(stat_ls)
    last_few_elem = arr[-compare_to_last:]
    final_converge_val = arr[-1]
    avg = np.mean(last_few_elem)
    pdiff = perc_diff(final_converge_val, avg)  # Returns a percentage, not a fraction.
    if pdiff > max_perc_diff:
        print(
            f"[{alarm}] A percent difference of {pdiff:.3f}% between the final convergence value ({final_converge_val})\n"
            f"  and the mean ({avg:.3f}) of the last {compare_to_last} fit values was found.\n"
            f"  However, a maximum percent difference of {max_perc_diff:.3f}% was specified.\n"
            f"  Subsequent fit values list:\n"
            f"    {stat_ls}\n"
            )
        if alarm in "ERROR":
            raise ValueError

def get_status_redchi2_fit(red_chi2, min_rc2, max_rc2):
    """Return True if: min_rc2 < reduced chi^2 value < max_rc2."""
    if (red_chi2 < min_rc2):
        print(f"Reduced chi^2 ({red_chi2}) is less than minimum allowed: {min_rc2}")
        return False
    elif (red_chi2 > max_rc2):
        print(f"Reduced chi^2 ({red_chi2}) is greater than minimum allowed: {max_rc2}")
        return False
    else:
        return True

def get_bestfit_vals_from_statsdict(d, check_convergence=False):
    """Return the best-fit mean, mean_err, std, std_err from a dict."""
    if check_convergence:
        raise RuntimeError("Jake, update this section!")
        # check_fit_convergence(FIXME: put stuff here)
    bf_mean     = d["mean_ls"][-1]
    bf_mean_err = d["mean_err_ls"][-1]
    bf_std      = d["std_ls"][-1]
    bf_std_err  = d["std_err_ls"][-1]
    return (bf_mean, bf_mean_err, bf_std, bf_std_err)