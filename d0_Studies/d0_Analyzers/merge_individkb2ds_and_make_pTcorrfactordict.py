"""Submit Iterated Gaussian Fit to SLURM

Purpose: This script opens a pickled KinBin2D with MyMuons
sorted into all the KB3Ds.
Iterated Gaussian fits (IGFs) are performed on the dpT/pT dist of each KB3D.
The fit results are saved to the KB3D.
The processed KB2D is then pickled into a new file.

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-14  # Happy pi day!
Updated: 2021-03-15
"""
import os
from pprint import pprint
from glob import glob
from Utils_ROOT.Printer import CanvasPrinter
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, make_dirs
from d0_Studies.d0_Utils.d0_fns import print_header_message
from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
from ParticleCollections import MyMuonCollection
from d0_Studies.kinematic_bins import equal_entry_bin_edges_eta_mod1_wholenum #, bin_edges_pT_sevenfifths_to1000GeV_wholenum

overwrite = 0

eta_range = equal_entry_bin_edges_eta_mod1_wholenum
# Grab all pTs.
inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/MC2016DY_individKB2D_withitergaussfitsonKB3Ds/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_*eta*_*pT*.pkl"

# eta_name = make_name_from_ls(eta_range, val="eta")
filename_base = "MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds"
# full_filename_base = f"{filename_base}_{eta_name}"
full_filename_base = f"{filename_base}"
outpkl_path_mucoll     = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/{full_filename_base}.pkl"
outpkl_path_pTcorrfact = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/{full_filename_base}_pTcorrfactors.pkl"
outpdf_path            = f"/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/final_muon_coll/{full_filename_base}.pdf"

def prep_area(outpkl_path, overwrite=False):
    """Make dirs, check overwrite."""
    outdir = os.path.dirname(outpkl_path)
    make_dirs(outdir)
    check_overwrite(outpkl_path, overwrite)

# def main():
printer = CanvasPrinter(show_plots=0)
printer.make_plots_pretty()
for f in (outpkl_path_mucoll, outpkl_path_pTcorrfact, outpdf_path):
    prep_area(f, overwrite=overwrite)
mu_coll = MyMuonCollection("DY2mu")
# inpkl_path_ls = glob(inpkl_path_template.replace("ETAPART", eta_name))
inpkl_path_ls = glob(inpkl_path_template)
# print("Found files:")
# pprint(inpkl_path_ls)
for pklfile in inpkl_path_ls:
    kb2d = open_pkl(pklfile)
    print(f"Opened pkl:\n{  pklfile}")
    mu_coll.KinBin2D_dict[kb2d.get_bin_key()] = kb2d
# mu_coll.make_KinBin2D_graphs()  # Should already be made.
print("Filled MyMuonCollection with KB2Ds.")
mu_coll.make_pT_corr_dict()  # Needs each kb2d in self.KinBin2D_dict to have interc_and_err and slope_and_err
pprint(f"Making dicts:\n  {[outpkl_path_pTcorrfact, outpkl_path_mucoll]}")
save_to_pkl(mu_coll.pT_corr_factor_dict, outpkl_path_pTcorrfact)
save_to_pkl(mu_coll, outpkl_path_mucoll)

    # mu_coll.kinbin3d_iterfitplot_ls.extend([kb3d.frame_dpTOverpT for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
    # mu_coll.get_all_plots()
    # mu_coll.kinbin3d_hist_ls.extend([kb3d.h_qd0 for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
# if __name__ == "__main__":
#     main()