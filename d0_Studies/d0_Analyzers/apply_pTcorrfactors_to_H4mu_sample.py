"""AdHoc pT Correction Applier

Purpose:
    Select "good" H->4mu events.
    Apply AdHoc/GeoFit pT correction procedure, per muon.
    Store the uncorr and corr muon pT in a '.root' file. 
Syntax:
    python script.py 
Notes:
    This code runs on a H->4mu sample.
    Should be used with Python 3.X.
    `inpkl_path_pTcorrfactors` should be a pickle of pT correction factors with the
    following structure:
    {
        'dpTOverpT_vs_qd0' : {
            '0.0eta0.2_10.0pT14.0' : {
                'intercept'     : -0.0004421981787913626,
                'intercept_err' : 4.7910538089331156e-05,
                'slope'         : 1.1863059538824643,
                'slope_err'     : 0.02705983833855872
            },
            '0.0eta0.2_100.0pT150.0' : {
                'intercept'     : ...,
            },
            ...,
        },
        'dpTOverpTscaled_vs_qd0' : {
            '0.0eta0.2_10.0pT14.0' : {
                ...
            },
            ...,
        }
        'dpTOverpTtimesavgOf1divpT_vs_qd0' : {
            '0.0eta0.2_10.0pT14.0' : {
                ...
            },
            ...,
        }
        'dpTOverpTtimesmuOf1divpT_vs_qd0' : {
            '0.0eta0.2_10.0pT14.0' : {
                ...
            },
            ...,
        }
    }

Author:  Jake Rosenzweig
Created: 2020-08-24
Updated: 2021-03-31
"""
from ParticleCollections import MyMuonCollection
from Utils_Python.Utils_Files import check_overwrite
import pickle

year = "2016"
overwrite = 0
verbose = 1
max_n_evts = -1
print_out_every = 100000

use_GeoFit_algo = True  # Will use Hmumu group's GeoFit factors and method.
force_zero_intercept = True

eta_lim = [0.0, 2.4]
pT_lim = [5, 200]
d0_lim = [0, 9999]
inv_m_lim = [0, 200]
# inv_m_lim = [105.0, 145.0]
dR_max = None #0.002

infile_path = "/cmsuf/data/store/user/t2/users/ferrico/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_2016.root"

# inpkl_path_pTcorrfactors = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/GeoFit/GeoFitTcorrfact_derivedfromMC2016_3etabins_0p0eta2p4_newestformat.pkl"
inpkl_path_pTcorrfactors = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/GeoFit/GeoFitTcorrfact_derivedfromMC2016WRONGactually2018_3etabins_0p0eta2p4.pkl"
# inpkl_path_pTcorrfactors = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY/pickles/final_muon_coll/MC2018DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds_new_pTcorrfactors.pkl"
# inpkl_path_pTcorrfactors = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/AdHoc/AdHocpTcorrfact_derivedfromMC2016DY_12etabins12pTbins_0p0eta2p4_5p0pT1000p0_newformat.pkl"
# inpkl_path_pTcorrfactors = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/MC2016_DY_dpTOverpT_1E6events/MC2016_DY_dpTOverpT_1E6events.pkl"
# inpkl_path_pTcorrfactors = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/GeoFit/GeoFitcorrfact_derivedfromMC2018_3etabins_0p0eta2p4.pkl"
# inpkl_path_pTcorrfactors = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/AdHoc/AdHocpTcorrfact_derivedfromMC2018DYJpsi_13etabins12pTbins_0p0eta2p4_5p0pT1000p0.pkl"
# inpkl_path_pTcorrfactors = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/CorrFactors/MC2018_d0_pT_corrfactors_0p0eta2p4_5p0pT1000p0_derivedfromDYJpsi_final.pkl"

# outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2016_m4mu_m4mucorrfromGeoFitfactors_fullstats_nodRcut.root"
outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2016_m4mu_m4mucorrfromWRONG2018GeoFitfactors_fullstats_zerointerc_chasebadval03.root"
# outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromGeoFitfactors_fullstats.root"
# outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromAdHocfactors_fullstats.root"
# outpath_rootfile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromGeoFitfactors_fullstats_noFSR_zerointerc.root"

if __name__ == "__main__":
    assert all(year in f for f in (infile_path, inpkl_path_pTcorrfactors, outpath_rootfile))
    if use_GeoFit_algo:
        assert "GeoFit" in inpkl_path_pTcorrfactors
        assert "AdHoc" not in outpath_rootfile
    else:
        # assert "AdHoc" in inpkl_path_pTcorrfactors
        assert "GeoFit" not in outpath_rootfile
    with open(inpkl_path_pTcorrfactors, "rb") as p:
        pT_corr_factor_dict = pickle.load(p)
    check_overwrite(outpath_rootfile, overwrite)
    mucoll = MyMuonCollection(["ggH"])
    mucoll.extract_muons(infile_path, prod_mode="H4mu",
                         n_evts=max_n_evts, n_evt_beg=None, n_evt_end=None,
                         print_out_every=print_out_every,
                         inv_m_lim=inv_m_lim,
                         eta_lim=eta_lim,
                         pT_lim=pT_lim,
                         d0_lim=d0_lim,
                         dR_max=dR_max,
                         do_mu_pT_corr=True,
                         force_zero_intercept=force_zero_intercept,
                         pT_corr_factor_dict=pT_corr_factor_dict,
                         correction_type=None,
                         use_GeoFit_algo=use_GeoFit_algo,
                         verbose=verbose)
    mucoll.write_m4muinfo_to_rootfile(outpath_rootfile, overwrite=overwrite)  #write_geofit_vals=True, )
    print("Done.")