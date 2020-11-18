"""Pickled dictionary merger

This code finds '.pkl' files with similar names,
opens up the pickled dictionary contained within,
and merges all these dictionaries into a single one
saving the final 'combined' pickled dict.

NOTE: infile_name_template should have the globbing characters
(e.g. '*') to grab all matching files.
"""
import pickle
import os
from glob import glob

#-----------------------#
#--- User Parameters ---#
#-----------------------#
infile_name_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/2017/ggH_m4mu/BestSoFar/MC2017_ggH_dpTOverpT_fullstats_*_5itergausfit_2p5sigs.pkl"
outdir = None    # If None, then put new dict in same place as src dicts.
namebase = "combined_muon_dct"  # If None, then name of combined pickle will be: 'combined_dcts.pkl'

#------------------------#
#--- Script Functions ---#
#------------------------#
def parse_paths(f, outname=None, outdir=None):
    """Return the absolute file path of where combined pickle will be placed."""
    name = "combined_dcts.pkl" if outname is None else outname
    outdir = os.path.split(f)[0] if outdir is None else outdir
    if not name.endswith(".pkl"):
        name = f"{name}.pkl"
    abs_file_path = os.path.join(outdir, name)
    return abs_file_path

def get_dct_from_pkl(f):
    """Return a single dictionary that has been pickled in file `f`."""
    with open(f, "rb") as inpkl:
        dct = pickle.load(inpkl)
    return dct

def save_pickle(obj, fullpath):
    """Save obj as a pickle at fullpath/name."""
    with open(fullpath, "wb") as pkl:
        print(f"[INFO] Writing pickle to:\n  {fullpath}")
        pickle.dump(obj, pkl, protocol=2)

# def sort_dct():
#     """
if __name__ == "__main__":
    fullpath = parse_paths(infile_name_template, outname=namebase, outdir=outdir)
    print(f"[INFO] Globbing files at:\n  {infile_name_template}")
    pkl_file_ls = glob(infile_name_template)
    # Make combined dict of KinBin2Ds.
    comb = {}
    for f in pkl_file_ls:
        dct = get_dct_from_pkl(f)
        comb.update(dct)
    # sorted_comb = sort_dct(comb)
    save_pickle(comb, fullpath)