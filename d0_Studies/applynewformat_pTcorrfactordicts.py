import pickle
from pprint import pprint
from glob import glob

inpkl_ls = glob("/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/AdHoc/*.pkl")
pprint(inpkl_ls)

for inpkl in inpkl_ls:
    with open(inpkl, "rb") as p:
        dct = pickle.load(p)
    print(f"{inpkl}:")
    print(list(dct.keys())[0])
    print(list(dct.values())[0])
#     if "dpTOverpT_vs_qd0" not in dct.keys():
#         old_dct_copy = dct.copy()
#         new_dct = { 'dpTOverpT_vs_qd0': {'0.0eta0.2_100.0pT150.0': {'intercept': 2.022759454757394e-06,
#    'intercept_err': 7.701238592076959e-07,
#    'slope': 0.03106886975917445,
#    'slope_err': 0.00044907629479497816},
#   '2.3eta2.4_27.0pT38.0': {'intercept': -0.00012372292642212856,
#    'intercept_err': 2.8639250676959453e-06,
#    'slope': 0.23176579952234125,
#    'slope_err': 0.0010737226458023294}
#         }
#         # Add the new style dict
#         new_filename = f"{inpkl.rstrip('.pkl')}_newerformat.pkl"

