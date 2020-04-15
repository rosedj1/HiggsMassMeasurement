import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
import sys
sys.path.append('/Users/Jake/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/')
sys.path.append('/Users/Jake/HiggsMassMeasurement/d0_Studies/')
from PyUtils.Utils_Files import makeDirs, make_str_title_friendly

def change_cmap_bkg_to_white(colormap, n=256):
    """
    Make any matplotlib colormap have a white background (i.e. any bins with 0 entries are colored white).
    This makes it easy to distinguish between bins with 0 entries and bins with only a few entries.
    
    Parameters
    ----------
    colormap : Colormap or str, optional
        A `.colors.Colormap` instance.  If not set, use rc settings.
    n : int
        The number of divisions of the color bar. 
        Use a higher number for a smoother transition of colors.
    """
    colors = cm.get_cmap(colormap, n)
    newcolors = colors(np.linspace(0, 1, n))
    white = np.array([1, 1, 1, 1])    # White background (Red, Green, Blue, Alpha).
    newcolors[0, :] = white    # Only change bins with 0 entries.
    newcmp = ListedColormap(newcolors)
    
    return newcmp

def save_plots_to_outpath(save_plot=False, outpath="", file_name="DEFAULT_NAME", save_as_png=False, verbose=False):
    if (save_plot):
        file_name = make_str_title_friendly(file_name)
        
        fullpath = os.path.join(outpath, file_name)
        
        makeDirs(outpath)
        
        plt.savefig(fullpath + '.pdf')
        if (verbose): print(f"Figure saved at: {fullpath + '.pdf'}")
        if (save_as_png):
            plt.savefig(fullpath + '.png')
            if (verbose): print(f"Figure saved at: {fullpath + '.png'}")
