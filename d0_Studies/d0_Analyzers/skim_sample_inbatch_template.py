"""Extract Muons from Sample

Skim a sample that was produced from the HZZ4LAnalyzer.
Store the muon info as a pickled MyMuonCollection (.pkl).
Muons are stored as MyMuons.

Author:  Jake Rosenzweig
Created: 2021-03-11
Updated: 2021-03-18
"""
import os
from ROOT import TFile
# Local imports.
from Utils_Python.Utils_Files import check_overwrite, make_dirs, save_to_pkl
from ParticleCollections import MyMuonCollection
#----- User Parameters -----#
# Input.
inpath_file = "INFILE_PATH"
prod_mode = "PROD_MODE"
# Output.
outdir_pkl = "OUTPKL_DIR"
filename_full = "FILENAME_FULL"

overwrite = OVERWRITE
verbose = VERBOSE
max_n_evts = None      # max_n_evts == -1 takes precedence.
n_evt_beg = N_EVT_BEG  # Then specifying a range takes precedence.
n_evt_end = N_EVT_END  
print_out_every = 1E6
#----- Functions -----#
def prep_area(filename_full, outdir_pkl, overwrite=False):
    """Return the file paths for the produced files."""
    new_file_name = f"{filename_full}"
    outpath_pkl = os.path.join(outdir_pkl, f"{new_file_name}.pkl")
    # See if these files already exist and overwrite if it's OK to do so.
    make_dirs(outdir_pkl, verbose)
    for f in [outpath_pkl]:
        check_overwrite(f, overwrite)
    return (outpath_pkl)

def main():
    f = TFile(inpath_file)
    t = f.Get("Ana/passedEvents")
    total_evts = t.GetEntries()
    # Prep your area.
    outpath_pkl = prep_area(filename_full, outdir_pkl, overwrite=overwrite)
    # Begin analysis.
    mu_coll = MyMuonCollection(prod_mode=prod_mode)

    mu_coll.extract_muons(inpath_file, prod_mode=prod_mode, n_evts=max_n_evts,
                                  n_evt_beg=n_evt_beg, n_evt_end=n_evt_end,
                                  print_out_every=print_out_every, eta_min=0.0, eta_max=2.4,
                                  pT_min=5, pT_max=1000,
                                  d0_max=0.01, dR_max=0.002,
                                  do_mu_pT_corr=False,
                                  force_zero_intercept=False,
                                  pT_corr_factor_dict=None,
                                  use_GeoFit_algo=False,
                                  verbose=verbose)
    save_to_pkl(mu_coll, outpath_pkl)

if __name__ == "__main__":
    main()