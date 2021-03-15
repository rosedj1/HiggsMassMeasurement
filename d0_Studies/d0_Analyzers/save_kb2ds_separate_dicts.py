"""KB2D Splitter

This script opens up a MyMuonCollection and saves each KB2D
as its own standalone dict.
Useful when a MyMuonCollection is too large to hold in RAM.

Syntax: python this_script.py

Author: Jake Rosenzweig
Created: 2021-03-14 # Happy pi day!
Updated: 

NOTE:
- kb2d.KinBin3D_dict  # Some dicts are empty!

TODO:
- This script should be implemented as a method in MyMuonCollection.
"""
import os
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, make_str_title_friendly, make_dirs

infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_fullstats_muoncoll_withkb3dbins.pkl"
subdir_name = "kb2d_individ_dicts"

file_suffix = ""
overwrite = 0

def make_outpath(infile_path, subdir_name, file_suffix, key, overwrite):
    """Return the absolute path of the outfile and make subdir if needed."""
    outdir_path, infile_name = os.path.split(infile_path)
    new_outdir_path = os.path.join(outdir_path, subdir_name)
    make_dirs(new_outdir_path, verbose=True)
    infile_name = infile_name.rstrip('.pkl')
    new_name = f"{infile_name}_{file_suffix}_{make_str_title_friendly(key)}.pkl"
    outpath = os.path.join(new_outdir_path, new_name)
    check_overwrite(outpath, overwrite=overwrite)
    return outpath

def main():
    muon_coll = open_pkl(infile_path)
    for ct, kb2d in enumerate(muon_coll.KinBin2D_dict.values()):
        # if ct == 1:
        #     break
        key = kb2d.get_bin_key(title_friendly=True)
        outpath = make_outpath(infile_path, subdir_name, file_suffix, key, overwrite)
        save_to_pkl(kb2d, outpath)

if __name__ == "__main__":
    main()