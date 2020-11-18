"""Muon pT Corrector

NOTE: This script has been superceded by
MuonCollection.extract_muons_from_Higgs_file().
Per muon pT correction can be performed inside that function.

Use a dictionary of pT correction factors
to correct the muon pT for a dictionary of MyMuon objects.

The updated dictionary is pickled in a new '.pkl'.
"""
import pickle
# from ParticleCollections import MyMuonCollection

inpkl_mu = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/Testing/combined_muon_dct.pkl"
inpkl_pT_corr_factor_dict = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/BestSoFar/combined_pT_corr_factor_dct.pkl"
verbose = False

eta_binedge_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
standard_pT_binedge_ls  = [5.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 1000.0]
truncated_pT_binedge_ls = [5.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 1000.0]

def open_pkls(inpkl_mu, inpkl_pT_corr_factor_dict):
    """Open and return dict of muon info and dict of pT corr factors."""
    with open(inpkl_mu, "rb") as p:
        mu_dct = pickle.load(p)
    with open(inpkl_pT_corr_factor_dict, "rb") as p:
        pT_corr_factor_dct = pickle.load(p)
    return (mu_dct, pT_corr_factor_dct)

def save_pkl(outpkl_mu, obj):
    """Save a dict of muon info with corr pT info into a new pickle."""
    outfile = outpkl_mu.replace(".pkl", "_withcorrections.pkl")
    with open(outfile, "wb") as p:
        print(f"Saving new pickle at:\n  {outfile}")
        pickle.dump(obj, p, protocol=2)

if __name__ == "__main__":
    mu_dct, pT_corr_factor_dct = open_pkls(inpkl_mu, inpkl_pT_corr_factor_dict)
    for count_x,kb2d in enumerate(mu_dct.values()):
        if verbose:
            print(f"Working on KinBin2D #{count_x}, eta_range={kb2d.eta_range}, pT_range={kb2d.pT_range}:")
        for count_y,kb3d in enumerate(kb2d.KinBin3D_dict.values()):
            print(f"  Working on KinBin3D #{count_x, count_y}, eta_min={kb3d.eta_min}:")
            # |eta| > 2 has fewer stats and so required different pT bins.
            if abs(kb3d.eta_min) >= 2.0:
                pT_binedge_ls = truncated_pT_binedge_ls
            else:
                pT_binedge_ls = standard_pT_binedge_ls
            # Correct muon pTs.
            corr_mu_ls = kb3d.correct_muon_ls_pT(pT_corr_factor_dct,
                                                 eta_binedge_ls, pT_binedge_ls, verbose=verbose)
            kb3d.muon_ls = corr_mu_ls
    save_pkl(inpkl_mu, mu_dct)