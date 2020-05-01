import numpy as np
import pandas as pd

#--------------------------------------------#
#--- Essential Particle Physics functions ---#
#--------------------------------------------#
def theta2pseudorap(theta):
    """Returns the pseudorapidity (eta), given the polar angle (theta)."""
    return -np.log( np.tan(theta/2) )

def pseudorap2theta(eta):
    """Returns polar angle (theta) between (-pi/2, pi/2) ."""
    return 2 * np.arctan( np.exp(-eta) )

def calc_dphi(phi2, phi1):
    """
    Calculate: phi2 - phi1 in the interval: [-pi, pi).
    Measure from lep1 to lep2. 
    Clockwise difference from 1 to 2 gives dphi < 0.
    Anticlockwise difference from 1 to 2 gives dphi > 0.
    
    Code converted from ROOT: https://root.cern.ch/doc/master/TVector2_8cxx_source.html#l00101
    
    FIXME: Implement arrays, in addition to Series and floats.
    
    Parameters
    ----------
    phi2 : float (or array-like)
        The azimuthal angle of, say, lepton 2. 
    phi1 : float (or array-like)
        The azimuthal angle of, say, lepton 1. 
        
    Returns
    -------
    dphi : float (or array-like)
        The difference in azimuthal angles in the range specified above.
    """
    # Using a Series. 
    if (isinstance(phi1, pd.core.series.Series) or 
        isinstance(phi2, pd.core.series.Series) 
       ):
        dphi = phi2 - phi1
        dphi = dphi.mask(dphi >= np.pi, dphi - 2*np.pi)
        dphi = dphi.mask(dphi < -np.pi, dphi + 2*np.pi)
        return dphi
    
    # Using floats.
    else:
        dphi = phi2 - phi1
        while (dphi >= np.pi): dphi -= 2*np.pi  
        while (dphi < np.pi): dphi += 2*np.pi
        return dphi

def calc_dR(deta, dphi):
    return np.sqrt( deta**2 + dphi**2 )
