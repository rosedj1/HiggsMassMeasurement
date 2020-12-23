"""Ad Hoc pT Correction Analyzer

FIXME: This doc string is old!

Purpose:
    Select "good" H->2mu events.
    Put each muon into User-specified (eta, pT, q*d0) bins,
    Make dpT/pT and q*d0 distributions in each bin.
    Perform iterated Gaussian fits of dpT/pT distributions to describe sigma of core. 
    all kinematic plots into a single PDF.
Syntax:
    python script.py > output.txt  (recommended)
    python script.py 
Notes:
    This code runs on a H->2mu sample.
    Make sure to check all the parameters in "User Parameters".
    Should be used with Python 3.X.

    This code makes the following distributions, before/after pT corr:
        (1) muon pT dist.
        (2) m4mu dist. (DSCB fit?)
        (3) eta dist., for each eta bin
        (4) muon qd0 dist, per (eta, pT) bin
    The following plots are produced:
Author:  Jake Rosenzweig
Created: 2020-08-24
Updated: 2020-12-17
"""
import pickle
import os
import ROOT as r
import numpy as np
# Local imports.
from Utils_Python.Utils_Files import check_overwrite, makeDirs
from Utils_ROOT.Printer import CanvasPrinter
from ParticleCollections import MyMuonCollection
from Particles import MyMuon
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum)
from d0_Studies.d0_Utils.d0_cls import KinBin2D
from d0_Studies.d0_Utils.d0_fns import find_bin_edges_of_value
#----- User Parameters -----#
# Input.
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2018.root"
inpkl_pT_corr_factor_dict = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Pickles/CorrFactorPkls/MC2018_Xunwu_pT_corrfactors_fromhismacro_5pT1000.pkl"
# Output.
outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018ggH_m4muvals_fullstats_usingXunwuHmumucorrfactorsfrommacro_withoutFSR.root"

pT_ls = [5.0, 1000.0]
overwrite = False
verbose = True

max_n_evts = -1
print_out_every = 50000

#----- Main -----#
# def main():
# Prep your area.
with open(inpkl_pT_corr_factor_dict, "rb") as p:
    pT_corr_factor_dict = pickle.load(p)

check_overwrite(outpath_rootfile, overwrite)
# Begin analysis.
muon_collection = MyMuonCollection()
muon_collection.extract_muons_from_H4mu_file(infile_path, n_evts=max_n_evts, print_out_every=print_out_every,
                                    eta_min=0.0, eta_max=2.4,
                                    pT_min=pT_ls[0], pT_max=pT_ls[1],
                                    d0_max=1,
                                    do_mu_pT_corr=True, 
                                    pT_corr_factor_dict=pT_corr_factor_dict,
                                    use_GeoFit_algo=True,
                                    verbose=True)

# Make root file with muon info.
muon_collection.write_m4muinfo_to_rootfile(outpath_rootfile, overwrite=overwrite) #write_geofit_vals=True, )