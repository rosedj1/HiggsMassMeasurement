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
Updated: 2021-04-09
"""
import argparse
import os
from Utils_Python.Utils_Files import open_pkl, make_dirs

#--- User Parameters ---#
# Will be superceded by passed-in arguments.
infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/pickles/MC2016DY2mu_skim_fullstats_verify.pkl"
year = "2016"
overwrite = 0
verbose = 1

parser = argparse.ArgumentParser(description='Separate KB2Ds.')
parser.add_argument('-v', '--verbose', dest='verbose', help="Print debug info.", action='store_true', default=verbose)
parser.add_argument('-o', '--overwrite', dest='overwrite', help="Overwrite existing files.", action='store_true', default=overwrite)
parser.add_argument('--year', dest='year', type=str, help='Year of sample.', default=year)
parser.add_argument('--infile', dest='infile', type=str, help='Path to input file.', default=infile)
args = parser.parse_args()

def main(args):
    infile = args.infile
    outdir = os.path.join(os.path.dirname(infile), "kb2d_dicts")
    assert all(year in path for path in (infile, outdir))
    make_dirs(outdir)
    print(f"...Opening pkl:\n  {infile}")
    muon_coll = open_pkl(infile)
    print(f"...Saving KinBin2Ds to individual pkls...")
    muon_coll.save_KB2Ds_separate_dcts(outdir=outdir, file_prefix="", overwrite=args.overwrite, verbose=args.verbose)

# infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p01_d0_1000p0_withGeoFitcorr.pkl"
# infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2018/DY/MC2018DY_fullstats_muoncoll_withkb3dbins.pkl"
# outdir      = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p01_d0_1000p0"
# file_prefix = "" #"MC2016DY_fullstats_muoncoll_withkb3dbins_nogenmatching_0p0_d0_0p01"

if __name__ == "__main__":
    main(args)