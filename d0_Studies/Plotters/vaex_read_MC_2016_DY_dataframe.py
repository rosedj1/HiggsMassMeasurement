import vaex  # The holy grail.
import numpy as np

infile_path_MC_2016_hdf_trimmed = "/Users/Jake/Desktop/MC_2016_trimmed.hdf5"
vdf = vaex.open(infile_path_MC_2016_hdf_trimmed)

# Make a concatenated VDF with muon1 and muon2 treated independently.
# Get names of columns for each muon.
col_ls = vdf.column_names
col_mu1_ls = [x for x in col_ls if x[-1] == "1"]
col_mu2_ls = [x for x in col_ls if x[-1] == "2"]
col_common_ls = [x for x in col_ls if x[-1] not in ["1","2"]]

msg = "[WARNING] There is a problem with number of columns or column naming!"
assert len(col_mu1_ls) + len(col_mu2_ls) + len(col_common_ls) == len(col_ls), msg

# Drop all the columns that end with "1" or "2".
vdf_mu1 = vdf.drop(col_mu1_ls + col_mu2_ls)
vdf_mu2 = vdf.drop(col_mu1_ls + col_mu2_ls)

#----- NOTE! -----#
# It is important to not add any new columns 
# to the VDF by this point. 

# Gives vdf_mu1 all the muon1 values and similarly for muon2.
for col_mu1, col_mu2 in zip(col_mu1_ls, col_mu2_ls):
    new_col_name = col_mu1[:-1]
    vdf_mu1[new_col_name] = getattr(vdf, col_mu1).values
    vdf_mu2[new_col_name] = getattr(vdf, col_mu2).values

# Put muon1 and muon2 info together into single vaex DataFrame.
vdf_concat = vaex.concat([vdf_mu1, vdf_mu2])
vdf_concat.drop("index", inplace=True)

#------------------------------------------#
#--- Save and manipulate new variables. ---#
#------------------------------------------#
# GEN info. 
eta_gen_ser = vdf_concat['genLep_eta']  # Good!
phi_gen_ser = vdf_concat['genLep_phi']  # Good!
pT_gen_ser  = vdf_concat['genLep_pt']  # Good!

# RECO info.
eta_rec_ser = vdf_concat['eta']  # Good!
phi_rec_ser = vdf_concat['phi']  # Good!
pT_rec_ser  = vdf_concat['pT']  # Good!

# Store other variables.
vdf_concat['delta_eta'] = deta_ser = eta_rec_ser - eta_gen_ser  # Good!
# Remember that delta_phi requires special treatment:
# -pi < delta_phi < pi
# vdf_concat['delta_phi'] = dphi_ser = vdf_concat.apply(calc_dphi, (phi_rec_ser, phi_gen_ser))  # Good!
vdf_concat['delta_phi'] = dphi_ser = phi_rec_ser - phi_gen_ser  # Good!
# vdf_concat['delta_R'] = dR_ser = vdf_concat.apply(calc_dR, (deta_ser, dphi_ser))
vdf_concat['delta_R'] = dR_ser = np.sqrt(deta_ser**2 + dphi_ser**2)
vdf_concat['delta_pT'] = dpT = pT_rec_ser - pT_gen_ser  # Good!
vdf_concat['delta_pToverGenpT'] = dpTratioGen = dpT / pT_gen_ser  # Good!
vdf_concat['delta_pToverRecpT'] = dpTratioRec = dpT / pT_rec_ser  # Good!
vdf_concat['delta_pToverGenpTsqred'] = dpTratioGen / pT_gen_ser  # Good! 
vdf_concat['delta_pToverRecpTsqred'] = dpTratioRec / pT_rec_ser  # Good!  
vdf_concat['qd0BS'] = vdf_concat['d0BS'] * vdf_concat['Id'] / -13.  # Good!
vdf_concat['qd0PV'] = vdf_concat['d0PV'] * vdf_concat['Id'] / -13.  # Good!

print("[INFO] Successfully retrieved vdf_concat!")
print("[INFO]   ...Analyzing MC 2016 DY...")