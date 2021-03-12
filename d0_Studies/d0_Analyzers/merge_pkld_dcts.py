"""Pickled dictionary merger

This code finds '.pkl' files with similar names,
opens up the pickled dictionary contained within,
and merges all these dictionaries into a single one
saving the final 'combined' pickled dict.

NOTE: infile_name_template should have globbing characters
(e.g. '*') to grab all matching files.

Author: Jake Rosenzweig
Created: <2021-02-05
Updated: 2021-03-10
"""
import pickle
import os
from glob import glob
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite

#-----------------------#
#--- User Parameters ---#
#-----------------------#
infile_name_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/MC2016_DY_dpTOverpT_1E6events/MC2016_DY_dpTOverpT_1E6events*sigs.pkl"
outdir = None    # If None, then put new dict in same place as src dicts.
namebase = "MC2016_DY_dpTOverpT_1E6events_12qd0reg_1000minmu_3itergausfit_2p5sigs.pkl"  # If None, then name of combined pickle will be: 'combined_dcts.pkl'
overwrite = 0
#------------------------#
#--- Script Functions ---#
#------------------------#
def make_full_filepath(f, outname=None, outdir=None):
    """Return the absolute file path of where combined pickle will be placed."""
    name = "combined_dcts.pkl" if outname is None else outname
    name = f"{name}.pkl" if not name.endswith(".pkl") else name
    outdir = os.path.dirname(f) if outdir is None else outdir
    abs_file_path = os.path.join(outdir, name)
    return abs_file_path

if __name__ == "__main__":
    fullpath = make_full_filepath(infile_name_template, outname=namebase, outdir=outdir)
    check_overwrite(fullpath, overwrite)
    print(f"[INFO] Globbing files at:\n  {infile_name_template}")
    pkl_file_ls = glob(infile_name_template)
    # Make combined dict of KinBin2Ds.
    comb = {}
    for f in pkl_file_ls:
        dct = open_pkl(f)
        comb.update(dct)
    save_to_pkl(comb, fullpath)