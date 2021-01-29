"""Ad hoc pT correction factor printer

"""
class PtCorrectionOrg:
    """Organizer for ad hoc pT correction factors."""

    def __init__(self, pT_corr_dct, sort_keys=True):
        """Initializes a pT correction organizer based on the dictionary values passed in.
        
        If sort_keys is True, then sorted_key_ls contains the keys
        sorted by eta and then pT. This does not affect the pT_corr_dct at all.

        Example:
          '1.2eta1.5_10.0pT20.0'
          '1.2eta1.5_5.0pT10.0' 
          '0.2eta0.4_5.0pT10.0' 
        becomes:
          '0.2eta0.4_5.0pT10.0' 
          '1.2eta1.5_5.0pT10.0' 
          '1.2eta1.5_10.0pT20.0'
        """
        self.pT_corr_dct = pT_corr_dct
        self.sorted_key_ls = self.sort_keys(pT_corr_dct) if sort_keys else list(pT_corr_dct.keys())
        self.sort_keys = sort_keys
        
    def sort_keys(self, d):
        """Return a list of keys sorted by eta and pT, numerically.
        
        Example key: '1.2eta1.5_10.0pT20.0'
        """
        srt_pT = lambda x : float(x.split("_")[1].split("pT")[0])
        srt_eta = lambda x : float(x.split("_")[0].split("eta")[0])
        # After some testing, it seems best to work *backwards*.
        # I.e. first sort pTs, then sort etas.
        sorted_by_pT = sorted(d, key=srt_pT)
        sorted_by_eta_and_pT = sorted(sorted_by_pT, key=srt_eta)
        return sorted_by_eta_and_pT

    # def parse_eta_pT_ranges(self):
    #     """Return a list of 2-tuples of all eta bins based on pT_corr_dct keys."""
    #     eta_bin_edges = []
    #     pT_bin_edges = []
    #     for key in self.pT_corr_dct:
    #         eta_part, pT_part = key.split("_")
    #         eta_min, eta_max = eta_part.split("eta")
            
        
    #     return tuple(set([(eta.)]))

    # def get_eta_ls(self):
    #     """Return a list of all eta bin values."""


import pickle
from pprint import pprint

infile_pkl = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Pickles/RerunEndcap_InclusiveEta_5pT1000w9bins_5itergausfit_2p5sigs/combined_dcts.pkl"

with open(infile_pkl,"rb") as inpkl:
    dct = pickle.load(inpkl)

ptcorr_org = PtCorrectionOrg(dct)
print("Sorted Dict:")
for key in ptcorr_org.sorted_key_ls:
    print(f"  {key} : {ptcorr_org.pT_corr_dct[key]}")

# pprint([for kb2d in pT_corr_dct.values() for kb3d in kb2d.KinBin3D_dict.values()])
#     : 
#         print(f"  eta_range={kb3d.eta_range}") 