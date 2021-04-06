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
    
    # Using floats or vaex expressions.
    else:
        dphi = phi2 - phi1
        while (dphi >= np.pi): dphi -= 2*np.pi  
        while (dphi < -np.pi): dphi += 2*np.pi
        return dphi

def calc_dR(deta, dphi):
    return np.sqrt( deta**2 + dphi**2 )

def perc_diff(num, ref):
    """
    Return the signed percent difference between num and ref:

    Returns: (num - ref) / ref * 100.0
    
    Parameters
    ----------
    num : float or array
        New number to which you want to compare to the reference number (ref).
        If num < ref, then will return a negative percent difference.
    ref : float or array
        Reference number. Goes in the denominator.
    """
    return (num - ref) / ref * 100.0

def calc_Hmass(mu_ls):
    """
    Return the invariant mass of the muons in mu_ls (len must be 2 or 4).

    NOTE: 
    - This method relies on calls like: mu.Pt(), mu.Eta(), etc. 
      So you must change the muon kinematics by doing: mu.SetPt(val), e.g.
    - The particles must have type ROOT.Math.LorentzVector.
    """
    len_ = len(mu_ls)
    if len_ == 2:
        H = mu_ls[0] + mu_ls[1]
    elif len_ == 4:
        H = mu_ls[0] + mu_ls[1] + mu_ls[2] + mu_ls[3]
    else:
        msg = f"mu_ls must have length 2 or 4, but has length {len_}."
        return ValueError(msg)
    return H.M()