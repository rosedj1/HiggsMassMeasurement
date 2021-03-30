"""KB2D Splitter

This script opens up a MyMuonCollection and saves each KB2D
as its own standalone dict.
Useful when a MyMuonCollection is too large to hold in RAM.

NOTE:
- If MyMuonCollection is stored inside a large dict,
then start a dev session on HPG first:
    srun --partition=bigmem --mem=128gb --ntasks=1 --cpus-per-task=8 --time=08:00:00 --pty bash -i
- kb2d.KinBin3D_dict  # Some dicts may be empty!

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-14 # Happy pi day!
Updated: 2021-03-29
"""
from Utils_Python.Utils_Files import open_pkl, make_dirs

year = "2016"
infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p0_d0_0p01.pkl"
# infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2018/DY/MC2018DY_fullstats_muoncoll_withkb3dbins.pkl"
outdir      = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dictsnofitinfo"
file_prefix = "MC2016DY_fullstats_muoncoll_withkb3dbins_nogenmatching_0p0_d0_0p01"
overwrite = 0
verbose = 1

if __name__ == "__main__":
    assert all(year in path for path in (infile_path, outdir, file_prefix))
    make_dirs(outdir)
    muon_coll = open_pkl(infile_path)
    muon_coll.save_KB2Ds_separate_dcts(outdir=outdir, file_prefix=file_prefix, overwrite=overwrite, verbose=verbose)