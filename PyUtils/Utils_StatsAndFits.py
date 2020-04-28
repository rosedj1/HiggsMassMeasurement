import numpy as np
from scipy.optimize import curve_fit

from d0_Utils.d0_fns import shift_binning_array, get_subset_mask
from d0_Utils.d0_dicts import color_dict, label_LaTeX_dict

from PyUtils.Utils_Plotting import make_stats_legend_for_gaus_fit

#--- Fitting Functions ---#
def gaussian(x, coeff, mu, sigma):#, normalize=False):
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
    rel_deviation = (x - mu) / sigma
    expon = -1 / 2 * np.float_power(rel_deviation, 2)
#    if (normalize):
#        coeff = 1 / np.sqrt( 2 * np.pi) / sigma  
    
    return coeff * np.exp(expon)

def fit_with_gaussian(x_vals, y_vals, guess_params=[1,0,1]):
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
    popt, pcov = curve_fit(gaussian, x_vals, y_vals, p0=guess_params)
    
    # FIXME: For some strange reason, sigma can turn out to be negative...
    popt[2] = np.abs(popt[2])
    
    popt_err = np.sqrt(np.diag(pcov))

    return popt, popt_err, pcov

def iterative_fit_gaus(iterations, bin_edges, bin_vals, first_mean=0, first_stdev=1, ax=None, draw_on_axes=True, verbose=True):
    """
    Fit a Gaussian function to the core of a distribution. 
    Do this iteratively to improve the fit of the core. 
    Draws the Gaus
    
    Parameters
    ----------
    iterations : int
        Number of fits to perform.
    bin_edges : array
        An array of the edges of the bins of the histogram being fit. 
    bin_vals : array 
        The y-values of each bin, where the first element contains the y-value of first bin, etc.
        Should have length = len(bin_edges) - 1.
    first_mean : float
        A guess at the mean of the first Gaussian fit.
    first_stdev : float
        A guess at the standard deviation of the first Gaussian fit.
        
    Returns 
    -------
    ax : axes object
        The original axes object, but now may have Gaussian fits drawn to it.
    """
    if (verbose):
        msg = "Performing {} iterative Gaussian fits".format(iterations)
        if (iterations == 1):
            msg = msg.replace("fits", "fit") 
        print(msg)

    bin_centers = shift_binning_array(bin_edges)

    count = 0
    popt = np.zeros(3)
    stats_dict = {
        'coeff_ls' : [],
        'coeff_err_ls' : [],
        'mean_ls' : [],
        'mean_err_ls' : [],
        'stdev_ls' : [],
        'stdev_err_ls' : [],
    }
    
    while count < iterations:
        count += 1
        
        if count == 1:
            # First fit: use original histogram's mean and stdev to choose a fit range.
            this_mean  = first_mean
            this_stdev = first_stdev
        else:
            # Otherwise use the last fit's optimized parameters.
            this_mean  = popt[1]
            this_stdev = popt[2]
        this_x_min = this_mean - 2*this_stdev
        this_x_max = this_mean + 2*this_stdev

        # Create 
        mask = get_subset_mask(bin_centers, x_min=this_x_min, x_max=this_x_max)

        new_bin_centers = bin_centers[mask]
        new_bin_vals = bin_vals[mask]

        popt, popt_err, pcov = fit_with_gaussian(new_bin_centers, new_bin_vals, guess_params=[1,this_mean,this_stdev])
        
        # Update stats dict.
        stats_dict['coeff_ls'].append(popt[0])
        stats_dict['mean_ls'].append(popt[1])
        stats_dict['stdev_ls'].append(popt[2])
        stats_dict['coeff_err_ls'].append(popt_err[0])
        stats_dict['mean_err_ls'].append(popt_err[1])
        stats_dict['stdev_err_ls'].append(popt_err[2])
        
        if (draw_on_axes):
            # Get the y-vals of the Gaussian fit for plotting
            gaus_y_vals = gaussian(new_bin_centers, *popt)

            leg_label_fit = make_stats_legend_for_gaus_fit(popt, popt_err)
            leg_label_fit = leg_label_fit.replace("Fit", "Gaus Fit {}:".format(count))
    
            if ax is None:
                f, ax = plt.subplots(figsize=(10,8))
            ax.plot(new_bin_centers, gaus_y_vals, color=color_dict[count+1], label=leg_label_fit, linestyle='-', marker="")
            ax.legend()
        
    return stats_dict, ax