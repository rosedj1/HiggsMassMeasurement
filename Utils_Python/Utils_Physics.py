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
    Return the signed percent difference between two numbers, relative to one of them.
    
    Parameters
    ----------
    num : float or array
        New number to which you want to compare to the reference number (ref).
        If num < ref, then will return a negative percent difference.
    ref : float or array
        Reference number. Goes in the denominator.
    """
    return (num - ref) / ref * 100.

# FIXME: Below is not yet tested.
class Particle:
    def __init__(self):
        pass
    
class Muon(Particle):
    def __init__(self, ID):
        self.ID = ID
        self.charge = self.charge_from_ID(ID)
    
    def SetPtEtaPhiMass(self, pT, eta, phi, mass):
        self.pT = pT
        self.eta = eta
        self.phi = phi
        self.mass = mass
        
    def SetGENPtEtaPhiMass(self, gpT, geta, gphi, gmass):
        self.gen_pT = gpT
        self.gen_eta = geta
        self.gen_phi = gphi
        self.gen_mass = gmass
                           
    def charge_from_ID(self, ID):
        if ID == 13:
            return -1 
        elif ID == -13: 
            return +1 
        else:
            raise TypeError("Based on the ID ({}), this is not a muon!".format(ID))