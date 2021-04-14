"""KB2D Merger

Finds pickled KB2Ds and merges them together into a MyMuonCollection.

NOTE:
- User has the option to save the pT correction factor dict
as a standalone pickled dict.
    asdf

Syntax: python this_script.py
Author: Jake Rosenzweig
Created: 2021-03-14  # Happy pi day!
Updated: 2021-04-12
"""
import os
from pprint import pprint
from glob import glob
from Utils_ROOT.Printer import CanvasPrinter
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, check_overwrite, make_dirs
from d0_Studies.d0_Utils.d0_fns import get_pT_part, parse_etapT_key
from Utils_Python.printing import print_header_message
from d0_Studies.d0_Analyzers.slurm_inbatch_derive_pTcorrfactors import make_name_from_ls
from ParticleCollections import MyMuonCollection
from natsort import natsorted

# Grab all pTs.
year = "2016"
prod_mode = "DY2mu"
inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/ApplyCorr/MC2016DY2mu/individKB2Ds_beforeaftercorr/pickles/*_withkb2dfits.pkl"
# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc_redo/pickles/*.pkl"
# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc/pickles/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc_*.pkl"
# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p01_d0_1000p0/unbinnedfit_widerwindow_fitwholerangefirstiter*.pkl"
# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dicts_beforeafterGeoFitcorr__0p0_d0_0p01/unbinnedfit_*.pkl"
# inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_0p0_d0_0p01/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_0p0_d0_0p01_*eta*_*pT*.pkl"
outpkl_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/"
overwrite = 0
# Make a separate pkl with the pT correction factor dict.
make_pTcorrdict = 0

outfilename_base = "muoncoll_withKB2Dfits_beforeaftercorr"
# outfilename_base = "muoncoll_itergaussfitsonKB3Ds_redo"
# outfilename_base = "muoncoll_itergaussfitsonKB2Ds_0p01_d0_1000p0_unbinned_widerwindow_fitwholerangefirstiter"
# outfilename_base = "MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds_0p0_d0_0p01"

def prep_area(d, overwrite=False):
    """Make directory `d` and check overwrite."""
    outdir = os.path.dirname(d)
    make_dirs(outdir)
    check_overwrite(d, overwrite)

def get_eta_min(key):
    """Return the leading eta value in the str `key` as a float.

    Example:
    key = "beforeaftercorr/pickles/0p6eta0p8_38p0pT50p0_withkb2dfits.pkl"
    returns: `0.6`
    """
    base = os.path.basename(key)
    nameparts_tup = base.split('_')
    for part in nameparts_tup:
        if "eta" in part:
            eta_part = part
            break
    eta_min = eta_part.split("eta")[0]
    return float(eta_min.replace("p", "."))

def main():
    outpkl_path_mucoll      = os.path.join(outpkl_dir, f"{outfilename_base}.pkl")
    outpkl_path_pTcorrfact  = os.path.join(outpkl_dir, f"{outfilename_base}_pTcorrfactors.pkl")
    outpdf_path             = os.path.join(outpkl_dir, f"{outfilename_base}.pdf")
    printer = CanvasPrinter(show_plots=0)
    printer.make_plots_pretty()
    for f in (outpkl_path_mucoll, outpkl_path_pTcorrfact, outpdf_path):
        prep_area(f, overwrite=overwrite)
    mu_coll = MyMuonCollection(prod_mode)
    # inpkl_path_ls = glob(inpkl_path_template.replace("ETAPART", eta_name))
    inpkl_path_ls = glob(inpkl_path_template)
    n_files = len(inpkl_path_ls)
    assert n_files > 0, '[ERROR] No files to glob!'
    # path_ls_etasorted = sorted(inpkl_path_ls)
    inpkl_path_ls_etapTsorted = natsorted(inpkl_path_ls)
    # May not be perfectly sorted still:
    # 1p5 comes before 1p25. Manually sort by eta vals.
    inpkl_path_ls_etapTsorted = sorted(inpkl_path_ls_etapTsorted, key=get_eta_min)
    print(f"Found {n_files} files:")
    pprint(inpkl_path_ls_etapTsorted)
    for pklfile in inpkl_path_ls_etapTsorted:
        kb2d = open_pkl(pklfile)
        print(f"Opened pkl:\n{  pklfile}")
        mu_coll.KinBin2D_dict[kb2d.get_bin_key()] = kb2d
    # mu_coll.make_KinBin2D_graphs()  # Should already be made.
    print("Filled MyMuonCollection with KB2Ds.")
    print("Making dicts:")
    if make_pTcorrdict:
        mu_coll.make_pT_corr_dict()  # Needs each kb2d in self.KinBin2D_dict to have interc_and_err and slope_and_err
        print("...Saving pT corr factor dict at:")
        save_to_pkl(mu_coll.pT_corr_factor_dict, outpkl_path_pTcorrfact, overwrite=overwrite)
    print("...Saving pickled MyMuonCollection at:")
    save_to_pkl(mu_coll, outpkl_path_mucoll, overwrite=overwrite)

if __name__ == "__main__":
    main()

    # mu_coll.kinbin3d_iterfitplot_ls.extend([kb3d.frame_dpTOverpT for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])
    # mu_coll.get_all_plots()
    # mu_coll.kinbin3d_hist_ls.extend([kb3d.h_qd0 for kb2d in self.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()])