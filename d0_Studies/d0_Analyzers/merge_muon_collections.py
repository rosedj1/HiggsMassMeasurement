import os
from glob import glob
from pprint import pprint
from ParticleCollections import MyMuonCollection
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum, bin_edges_pT_sevenfifths_to1000GeV_wholenum

eta_ls = equal_entry_bin_edges_eta_mod1_wholenum
pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum

n_bins_dpTOverpT = 100
x_lim_dpTOverpT = [-0.4, 0.4]
n_bins_qd0 = 100
x_lim_qd0 = [-0.01, 0.01]
verbose = True
overwrite = 0

indir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2018/DY/MC2018DY_skim_fullstats_inbatch"
outfile = "MC2018DY_skim_fullstats"
glob_template = "MC2018DY_skim_fullstats_skim_sample_inbatch_template_copy*.pkl"

if ".pkl" not in outfile:
    outfile += ".pkl"
outpkl_path = os.path.join(indir, outfile)
check_overwrite(outpkl_path, overwrite=overwrite)
# Glob files.
pkl_file_ls = glob(os.path.join(indir, glob_template))
print(f"Globbed {len(pkl_file_ls)} files:")
pprint(pkl_file_ls)
# Empty muon collection to which all others will add.
mu_coll = MyMuonCollection("DY2mu")
for ct, p in enumerate(pkl_file_ls, 1):
    print(f"...Opening file {ct}:\n{  p}")
    mu_coll_part = open_pkl(p)
    mu_coll.m4mu_ls.extend(mu_coll_part.m4mu_ls)
    mu_coll.muon_ls.extend(mu_coll_part.muon_ls)
    del mu_coll_part

# Can also make inclusive kinematic plots:
mu_coll.make_inclusive_kinematic_plots()
# Sort muons into specified (eta, pT, qd0) bins.
mu_coll.sort_muons(eta_ls, pT_ls, pT_corr_factor_dict=None,
                   n_bins_dpTOverpT=n_bins_dpTOverpT, x_lim_dpTOverpT=x_lim_dpTOverpT,
                   n_bins_qd0=n_bins_qd0, x_lim_qd0=x_lim_qd0, verbose=verbose)
save_to_pkl(mu_coll, outpkl_path=os.path.join(indir, outfile))
