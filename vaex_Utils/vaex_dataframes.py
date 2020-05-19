"""
A file to store vaex functions and filepaths.
"""
import vaex  # The holy grail.
import numpy as np

vdf_MC_2016_DY = vaex.open("/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/NTuples/MC/2016/MC_2016_DY.hdf5")
vdf_MC_2017_DY = vaex.open("/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/NTuples/MC/2017/MC_2017_DY.hdf5")
vdf_MC_2017_Jpsi = vaex.open("/Users/Jake/Desktop/Research/Higgs_Mass_Measurement/NTuples/MC/2017/MC_2017_Jpsi.hdf5")

def df_drop_cols(df, col_keep_ls, inplace=True):
    import pandas
    """
    Drop all columns in a vaex or pandas DataFrame (df) 
    that are NOT specified in the given list (col_keep_ls).
    
    Parameters
    ----------
    df : vaex.hdf5.dataset.Hdf5MemoryMapped or pandas.core.frame.DataFrame
        DataFrame whose columns will be dropped.
    col_keep_ls : list of str
        The names of the columns that will NOT be dropped. 
    inplace : bool
        If True, then permanently drop the columns from df. 
    """
    if isinstance(df, vaex.hdf5.dataset.Hdf5MemoryMapped):
        all_col_names = set(df.column_names)
    elif isinstance(df, pandas.core.frame.DataFrame):
        all_col_names = set(df.columns)
        
    col_keep_set = set(col_keep_ls)
    col_drop_set = all_col_names - col_keep_set
    
    for col in col_drop_set:
        if isinstance(df, vaex.hdf5.dataset.Hdf5MemoryMapped):
            df.drop(col, inplace=inplace)
        elif isinstance(df, pandas.core.frame.DataFrame):
            df.drop(col, axis=1, inplace=inplace)

def prepare_vaex_df(vdf):    
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

    # Give each row a number. 
    vdf_concat["index"] = np.arange(vdf_concat.count())
    
    print("[INFO] Successfully retrieved vdf_concat!")
    print("[INFO] vdf_concat has the following columns:\n{}\n".format(vdf_concat.column_names))
    
    return vdf_concat

def vaex_apply_masks(vdf, eta_minmax, pT_minmax, qd0_minmax, massZ_minmax, dR_max):
    print("Applying the following cuts to the VDF:")
    cut_str  = '{} < vdf["eta"].abs()) & (vdf["eta"].abs() < {})\n'.format(eta_minmax[0], eta_minmax[1])
    cut_str += '{} < vdf["pT"]) & (vdf["pT"] < {})\n'.format(pT_minmax[0], pT_minmax[1])
    cut_str += '{} < vdf["qd0BS"]) & (vdf["qd0BS"] < {})\n'.format(qd0_minmax[0], qd0_minmax[1])
    cut_str += '{} < vdf["massZ"]) & (vdf["massZ"] < {})\n'.format(massZ_minmax[0], massZ_minmax[1])
    cut_str += 'vdf["delta_R"] < {})'.format(dR_max)
    print(cut_str)
    
    mask_eta = (eta_minmax[0] < vdf["eta"].abs()) & (vdf["eta"].abs() < eta_minmax[1])
    mask_pT = (pT_minmax[0] < vdf["pT"]) & (vdf["pT"] < pT_minmax[1])
    mask_qd0 = (qd0_minmax[0] < vdf["qd0BS"]) & (vdf["qd0BS"] < qd0_minmax[1])
    mask_massZ = (massZ_minmax[0] < vdf["massZ"]) & (vdf["massZ"] < massZ_minmax[1])
    mask_dR = (vdf["delta_R"] < dR_max)

    all_masks = mask_eta & mask_pT & mask_qd0 & mask_massZ & mask_dR
    
    return all_masks

def convert_npz_to_hdf5(inpath_npz, outpath_hdf5, col_keep_ls):
    """
    FIXME: 
      [ ] Test this function first!
      [ ] Get this to run on a remote machine. 
    
    Convert a .npz file to a vaex-friendly .hdf5 file.
    This is useful for vaex to instantly read .hdf5 files into 
    a vaex DataFrame. 
    
    NOTES: 
      - This function goes through a couple DataFrame conversions and 
        because of this consumes a lot of memory. I haven't found a
        way around it yet. 
      - For a 6 GB .npz file, this function should take about 2 hrs to finish. 
      
    Parameters
    ----------
    inpath_npz : .npz file
        Absolute file path to .npz file. 
    outpath_hdf5 : .hdf5 file
        Absolute file path where you want the newly created .hdf5 to be stored.
    col_keep_ls : list of str
        The names of the columns that will be stored in a DataFrame. 
    """
    df = pd.DataFrame(np.load(inpath_npz)['arr_0'])
    print("[INFO] Loaded pandas DataFrame from npz.")
    df_drop_cols(df, col_keep_ls)
    print("[INFO] Dropped pandas DataFrame columns not specified in col_keep_ls.")
    vdf = vaex.from_pandas(df)
    print("[INFO] Converted pandas DF to vaex DF.")
    del df
    print("[INFO] Deleted pandas DF from memory.")
    vdf.export_hdf5(outpath_hdf5)
    print("[INFO] Created hdf5 file:\n  {}".format(outpath_hdf5))
    del vdf
    print("[INFO] Deleted vaex DF from memory.")
    print("* * * Conversion complete. * * *")