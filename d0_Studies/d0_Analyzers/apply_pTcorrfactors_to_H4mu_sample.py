"""Ad Hoc pT Correction Analyzer

Purpose:
    Select "good" H->4mu events.
    Apply AdHoc/GeoFit pT correction procedure, per muon.
    Store the uncorr and corr muon pT in a '.root' file. 
Syntax:
    python script.py 
Notes:
    This code runs on a H->4mu sample.
    Should be used with Python 3.X.
Author:  Jake Rosenzweig
Created: 2020-08-24
Updated: 2021-02-17
"""
from ParticleCollections import MyMuonCollection
from Utils_Python.Utils_Files import check_overwrite
import pickle

overwrite = 0
verbose = 1
max_n_evts = -1
print_out_every = 100000

use_GeoFit_algo = True  # Will use Hmumu group's GeoFit factors and method.
force_zero_intercept = False

infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2018.root"
inpkl_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/GeoFit/GeoFitcorrfact_derivedfromMC2018_3etabins_0p0eta2p4.pkl"
# inpkl_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/AdHoc/AdHocpTcorrfact_derivedfromMC2018DYJpsi_13etabins12pTbins_0p0eta2p4_5p0pT1000p0.pkl"
# inpkl_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/MC2018_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0_derivedfromDYJpsi_final.pkl"
outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromGeoFitfactors_fullstats.root"
# outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromAdHocfactors_fullstats.root"
# outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromGeoFitfactors_fullstats_noFSR_zerointerc.root"

if __name__ == "__main__":
    if use_GeoFit_algo:
        assert "GeoFit" in inpkl_path
        assert "AdHoc" not in outpath_rootfile
    else:
        assert "AdHoc" in inpkl_path
        assert "GeoFit" not in outpath_rootfile
    with open(inpkl_path, "rb") as p:
        pT_corr_factor_dict = pickle.load(p)
    check_overwrite(outpath_rootfile, overwrite)
    mucoll = MyMuonCollection(["ggH"])
    mucoll.extract_muons_from_H4mu_file(infile_path, n_evts=max_n_evts, print_out_every=print_out_every,
                                        eta_min=0.0, eta_max=2.4,
                                        pT_min=5, pT_max=1000,
                                        d0_max=1,
                                        do_mu_pT_corr=True,
                                        force_zero_intercept=force_zero_intercept,
                                        pT_corr_factor_dict=pT_corr_factor_dict,
                                        use_GeoFit_algo=use_GeoFit_algo,
                                        verbose=verbose)
    mucoll.write_m4muinfo_to_rootfile(outpath_rootfile, overwrite=overwrite)  #write_geofit_vals=True, )
    print("Done.")