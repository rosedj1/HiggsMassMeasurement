"""Extract Muons from Sample

Skim a sample that was produced from the HZZ4LAnalyzer.
Store the muon info as a pickled MyMuonCollection (.pkl).
Muons are stored as MyMuons.

Author:  Jake Rosenzweig
Created: 2021-03-11
Updated: 
"""
import os
# Local imports.
from Utils_Python.Utils_Files import check_overwrite, make_dirs, save_to_pkl
from ParticleCollections import MyMuonCollection
#----- User Parameters -----#
# Input.
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/DY_JPsi_Upsilon/DY_2017.root"
prod_mode = "DY2mu"
# Output.
outdir_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2017/DY/test/"
outfile_prefix = "MC2017DY_skim_fullstats"

overwrite = 1
verbose = 1
max_n_evts = -1
print_out_every = 1E6
#----- Functions -----#
def prep_area(outfile_prefix, outdir_pkl, overwrite=False):
    """Return the file paths for the produced files."""
    new_file_name = f"{outfile_prefix}"
    outpath_pkl = os.path.join(outdir_pkl, f"{new_file_name}.pkl")
    # See if these files already exist and overwrite if it's OK to do so.
    make_dirs(outdir_pkl, verbose)
    for f in [outpath_pkl]:
        check_overwrite(f, overwrite)
    return (outpath_pkl)

def main():
    # Prep your area.
    outpath_pkl = prep_area(outfile_prefix, outdir_pkl, overwrite=overwrite)
    # Begin analysis.
    mu_coll = MyMuonCollection(prod_mode=prod_mode)
    mu_coll.extract_muons(infile_path, prod_mode=prod_mode, n_evts=max_n_evts,
                                  print_out_every=print_out_every, eta_min=0.0, eta_max=2.4,
                                  pT_min=5, pT_max=200,
                                  d0_max=1, dR_max=0.002,
                                  do_mu_pT_corr=False,
                                  force_zero_intercept=False,
                                  pT_corr_factor_dict=None,
                                  use_GeoFit_algo=False,
                                  verbose=verbose)
    save_to_pkl(mu_coll, outpath_pkl)

if __name__ == "__main__":
    main()