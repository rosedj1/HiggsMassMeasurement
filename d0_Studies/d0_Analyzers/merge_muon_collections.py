"""MyMuonCollection Merger and Muon Sorter

This code globs all MyMuonCollections stored in pickles,
merges their muon_ls (MyMuons) and m4mu_ls together into
a large MyMuonCollection (mu_coll).

Then mu_coll makes hists of all inclusive muon kinematics.
NOTE: This script doesn't draw plots!

Finally, mu_coll sorts all MyMuons into their respective KB2Ds
based on the provided eta and pT bin edges.
MyMuonCollection.muon_ls gets overwritten,
but muons are stored in each KB2D.

KB3Ds do not get made yet.

Syntax: python <this.py>
Author: Jake Rosenzweig
Created: 2021-03-18
Updated: 2021-03-29
"""
import os
from glob import glob
from pprint import pprint
from ParticleCollections import MyMuonCollection
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite
from d0_Studies.kinematic_bins import (equal_entry_bin_edges_eta_mod1_wholenum,
                                       bin_edges_pT_sevenfifths_to1000GeV_wholenum,
                                       eta_bins_geofitvsVXBS, pT_bins_geofitvsVXBS)

year = "2016"
prod_mode = "DY2mu"
eta_ls = eta_bins_geofitvsVXBS #equal_entry_bin_edges_eta_mod1_wholenum
pT_ls = pT_bins_geofitvsVXBS #bin_edges_pT_sevenfifths_to1000GeV_wholenum

d0_lim = [0.0, 0.01]  # Only keep muons in this d0 bin.
# d0_lim = [0.01, 1000.0]  # Only keep muons in this d0 bin.

n_bins_dpTOverpT = 100
x_lim_dpTOverpT = [-0.4, 0.4]
n_bins_qd0 = 100
x_lim_qd0 = [-0.01, 0.01]
verbose = 1
overwrite = 0

indir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching"
# indir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2018/DY/MC2018DY_skim_fullstats_inbatch"
glob_template = "MC2016DY_skim_fullstats_nogenmatching_skim_sample_inbatch_template_copy_*.pkl"
outfile = "MC2016DY_skim_fullstats_nogenmatching_0p0_d0_0p01"
# outfile = "MC2016DY_skim_fullstats_nogenmatching_0p01_d0_1000p0"

assert all(year in f for f in (indir, outfile, glob_template))
if ".pkl" not in outfile:
    outfile += ".pkl"
outpkl_path = os.path.join(indir, outfile)
check_overwrite(outpkl_path, overwrite=overwrite)
# Glob files.
pkl_file_ls = glob(os.path.join(indir, glob_template))
assert len(pkl_file_ls) > 0, "No files were globbed."
print(f"Globbed {len(pkl_file_ls)} files:")
pprint(pkl_file_ls)
# Empty muon collection to which all others will add.
mu_coll = MyMuonCollection(prod_mode)
for ct, p in enumerate(pkl_file_ls, 1):
    print(f"...Opening file {ct}:\n{  p}")
    mu_coll_part = open_pkl(p)
    for ct, mu in enumerate(mu_coll_part.muon_ls):
        # Only look for muons in |d0| bin.
        if not (d0_lim[0] < abs(mu.d0) and abs(mu.d0) < d0_lim[1]):
            continue
        mu_coll.muon_ls.extend([mu])
        # Only grab the nth invariant mass of the event.
        # n = number of muons in final state.
        if (ct % 2) == 0:
            assert "DY" in prod_mode
            mu_coll.m4mu_ls.extend([mu.inv_m_event])
    # mu_coll.m4mu_ls.extend(mu_coll_part.m4mu_ls)
    # mu_coll.muon_ls.extend(mu_coll_part.muon_ls)
    del mu_coll_part

# Can also make inclusive kinematic plots:
mu_coll.make_inclusive_kinematic_plots()
# Sort muons into specified (eta, pT) bins.
# This overwrites MyMuonCollection's muon_ls.
mu_coll.sort_muons(eta_ls, pT_ls, pT_corr_factor_dict=None,
                   n_bins_dpTOverpT=n_bins_dpTOverpT, x_lim_dpTOverpT=x_lim_dpTOverpT,
                   n_bins_qd0=n_bins_qd0, x_lim_qd0=x_lim_qd0, verbose=verbose)
save_to_pkl(mu_coll, outpkl_path=os.path.join(indir, outfile), overwrite=overwrite)
