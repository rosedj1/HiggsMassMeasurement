import os
import pickle
from pprint import pprint
from glob import glob
from Utils_Python.Utils_Files import open_pkl, save_to_pkl, save_to_json, check_overwrite
from natsort import natsorted

overwrite = 0
inpkl_ls = glob("/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/AdHoc/AdHocpTcorrfact_derivedfromMC*.pkl")
pprint(inpkl_ls)

new_dct_format = {
            "dpTOverpT_vs_qd0"                 : {},
            "dpTOverpTscaled_vs_qd0"           : {},
            "dpTOverpTtimesavgOf1divpT_vs_qd0" : {},
            "dpTOverpTtimesmuOf1divpT_vs_qd0"  : {},
}

for inpkl in inpkl_ls:
    new_outfile_path = inpkl.replace('.pkl', '_newerformat.pkl')
    check_overwrite(new_outfile_path, overwrite)
    with open(inpkl, "rb") as p:
        dct = pickle.load(p)
    print(f"...Opened {inpkl}.")
    print("Original dct:")
    pprint(dct)
    print("Sorted dct:")
    srt_key_ls = natsorted(dct)
    
    pprint(srt_dct)
    if "dpTOverpT_vs_qd0" not in srt_dct.keys():
        print(f"dpTOverpT_vs_qd0 not found in keys. Updating dict.")
        updated_dct = {
            "dpTOverpT_vs_qd0"                 : srt_dct,
            "dpTOverpTscaled_vs_qd0"           : {},
            "dpTOverpTtimesavgOf1divpT_vs_qd0" : {},
            "dpTOverpTtimesmuOf1divpT_vs_qd0"  : {},
        }
        save_to_pkl(updated_dct, new_outfile_path, overwrite=overwrite)
        save_to_json(updated_dct, new_outfile_path.replace(".pkl", ".json"), overwrite=overwrite)
        break