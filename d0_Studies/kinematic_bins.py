# Equal-entry bins: 
# Determined by counting total number of entries (N) in dataset,
# ordering values of (eta, pT or q*d0) from least to greatest, 
# and determining bin edge where N/K is the same for all K regions.
# NOTE:
#--- Anything labelled with "_mod" is not guaranteed to 
#--- have equal entries across all bins.
equal_entry_bin_edges_eta = [0.0, 0.2, 0.395, 0.6, 0.805, 1.02, 1.25, 1.49, 1.75, 2.05, 2.4]
equal_entry_bin_edges_eta_mod1 = [0.0, 0.2, 0.395, 0.6, 0.805, 1.02, 1.25, 1.49, 1.75, 2.0, 2.1, 2.2, 2.3, 2.4]
equal_entry_bin_edges_eta_mod1_wholenum = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
equal_entry_bin_edges_eta_mod2 = [0.0, 0.2, 0.395, 0.6, 0.805, 1.02, 1.25, 1.49, 1.75, 2.0, 2.2, 2.4]
bin_edges_eta_outer_strip_layer = [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]

equal_entry_bin_edges_pT  = [5.0, 23.01, 29.06, 33.45, 37.05, 40.18, 43.04, 45.99, 50.46, 60.77, 200.0]
equal_entry_bin_edges_pT_mod1 = [5.0, 23.01, 29.06, 33.45, 37.05, 40.18, 43.04, 45.99, 50.46, 60.77, 200.0]
equal_entry_bin_edges_pT_mod2 = [5.0, 14, 17, 20, 25, 35, 40, 45, 50, 55, 60, 85, 200.0]
equal_entry_bin_edges_pT_sevenfifths = [5.0, 7.0, 9.8, 13.720000000000002, 19.208000000000006, 26.89120000000001, 37.647680000000015, 52.706752000000016, 73.78945280000002, 103.30523392000002]
equal_entry_bin_edges_pT_sevenfifths_mod = [5.0, 7.0, 9.8, 13.720000000000002, 19.208000000000006, 26.89120000000001, 37.647680000000015, 52.706752000000016, 73.78945280000002, 103.30523392000002]
equal_entry_bin_edges_pT_sevenfifths_to1000GeV = [5.0, 7.0, 9.8, 13.720000000000002, 19.208000000000006, 26.89120000000001, 37.647680000000015, 52.706752000000016, 73.78945280000002, 103.30523392000002, 144.627, 202.478, 1000]
bin_edges_pT_sevenfifths_to1000GeV_wholenum = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000]

equal_entry_bin_edges_qd0 = [-0.20014, -0.00257, -0.00166, -0.00103, -0.0005, -1e-05, 0.00048, 0.001, 0.00163, 0.00254, 0.19913]





#--- Stuff below was never fully developed. ---#
#   n_evts = 10000
#   
#   dR_cut=dR_cut = 0.02
#   
#   eta_bin_BARREL  = [0.0, 0.3]
#   eta_bin_OVERLAP = [0.8, 1.1]
#   eta_bin_ENDCAP  = [2.1, 2.4]
#   
#   pT_bin_LOW1 = [5, 7]
#   pT_bin_LOW2 = [7, 10]
#   pT_bin_LOW3 = [10, 15]
#   pT_bin_LOW4 = [15, 20]
#   pT_bin_LOW5 = [20, 25]
#   pT_bin_LOW6 = [25, 30]
#   
#   pT_bin_MED1 = [30, 35]
#   pT_bin_MED2 = [35, 40]
#   pT_bin_MED3 = [40, 45]
#   pT_bin_MED4 = [45, 50]
#   pT_bin_MED5 = [50, 55]
#   pT_bin_MED6 = [55, 60]
#   
#   pT_bin_HIGH1 = [60, 65]
#   pT_bin_HIGH2 = [65, 70]
#   pT_bin_HIGH3 = [70, 75]
#   pT_bin_HIGH4 = [75, 80]
#   pT_bin_HIGH5 = [80, 90]
#   pT_bin_HIGH6 = [90, 100]
#   
#   #--------------#
#   #--- Barrel ---#
#   #--------------#
#   # Low pT.
#   kbin_eta_BARREL_pT_bin_LOW1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_LOW1, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_LOW2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_LOW2, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_LOW3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_LOW3, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_LOW4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_LOW4, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_LOW5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_LOW5, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_LOW6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_LOW6, dR_cut=dR_cut)
#   
#   # Med pT.
#   kbin_eta_BARREL_pT_bin_MED1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_MED1, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_MED2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_MED2, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_MED3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_MED3, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_MED4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_MED4, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_MED5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_MED5, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_MED6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_MED6, dR_cut=dR_cut)
#   
#   # High pT.
#   kbin_eta_BARREL_pT_bin_HIGH1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_HIGH1, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_HIGH2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_HIGH2, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_HIGH3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_HIGH3, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_HIGH4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_HIGH4, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_HIGH5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_HIGH5, dR_cut=dR_cut)
#   kbin_eta_BARREL_pT_bin_HIGH6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_BARREL, pT_cut_ls=pT_bin_HIGH6, dR_cut=dR_cut)
#   
#   
#   #---------------#
#   #--- Overlap ---#
#   #---------------#
#   # Low pT.
#   kbin_eta_OVERLAP_pT_bin_LOW1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_LOW1, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_LOW2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_LOW2, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_LOW3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_LOW3, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_LOW4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_LOW4, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_LOW5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_LOW5, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_LOW6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_LOW6, dR_cut=dR_cut)
#   
#   # Med pT.
#   kbin_eta_OVERLAP_pT_bin_MED1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_MED1, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_MED2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_MED2, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_MED3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_MED3, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_MED4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_MED4, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_MED5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_MED5, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_MED6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_MED6, dR_cut=dR_cut)
#   
#   # High pT.
#   kbin_eta_OVERLAP_pT_bin_HIGH1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_HIGH1, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_HIGH2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_HIGH2, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_HIGH3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_HIGH3, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_HIGH4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_HIGH4, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_HIGH5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_HIGH5, dR_cut=dR_cut)
#   kbin_eta_OVERLAP_pT_bin_HIGH6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_OVERLAP, pT_cut_ls=pT_bin_HIGH6, dR_cut=dR_cut)
#   
#   #---------------#
#   #--- Endcap ---#
#   #---------------#
#   # Low pT.
#   kbin_eta_ENDCAP_pT_bin_LOW1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_LOW1, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_LOW2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_LOW2, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_LOW3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_LOW3, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_LOW4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_LOW4, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_LOW5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_LOW5, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_LOW6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_LOW6, dR_cut=dR_cut)
#   
#   # Med pT.
#   kbin_eta_ENDCAP_pT_bin_MED1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_MED1, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_MED2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_MED2, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_MED3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_MED3, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_MED4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_MED4, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_MED5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_MED5, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_MED6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_MED6, dR_cut=dR_cut)
#   
#   # High pT.
#   kbin_eta_ENDCAP_pT_bin_HIGH1 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_HIGH1, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_HIGH2 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_HIGH2, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_HIGH3 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_HIGH3, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_HIGH4 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_HIGH4, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_HIGH5 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_HIGH5, dR_cut=dR_cut)
#   kbin_eta_ENDCAP_pT_bin_HIGH6 = KinemBinnedEtaPt(df_MC_2016, n_evts=n_evts, eta_cut_ls=eta_bin_ENDCAP, pT_cut_ls=pT_bin_HIGH6, dR_cut=dR_cut)
#   
#   # LOW pT.
#   kbin_ls_eta_INCLUSIVE_pT_LOW1 = [
#       [kbin_eta_BARREL_pT_bin_LOW1, kbin_eta_OVERLAP_pT_bin_LOW1, kbin_eta_ENDCAP_pT_bin_LOW1],
#       [kbin_eta_BARREL_pT_bin_LOW2, kbin_eta_OVERLAP_pT_bin_LOW2, kbin_eta_ENDCAP_pT_bin_LOW2],
#       [kbin_eta_BARREL_pT_bin_LOW3, kbin_eta_OVERLAP_pT_bin_LOW3, kbin_eta_ENDCAP_pT_bin_LOW3],
#   ]
#   
#   kbin_ls_eta_INCLUSIVE_pT_LOW2 = [
#       [kbin_eta_BARREL_pT_bin_LOW4, kbin_eta_OVERLAP_pT_bin_LOW4, kbin_eta_ENDCAP_pT_bin_LOW4],
#       [kbin_eta_BARREL_pT_bin_LOW5, kbin_eta_OVERLAP_pT_bin_LOW5, kbin_eta_ENDCAP_pT_bin_LOW5],
#       [kbin_eta_BARREL_pT_bin_LOW6, kbin_eta_OVERLAP_pT_bin_LOW6, kbin_eta_ENDCAP_pT_bin_LOW6],
#   ]
#   
#   # MED pT.
#   kbin_ls_eta_INCLUSIVE_pT_MED1 = [
#       [kbin_eta_BARREL_pT_bin_MED1, kbin_eta_OVERLAP_pT_bin_MED1, kbin_eta_ENDCAP_pT_bin_MED1],
#       [kbin_eta_BARREL_pT_bin_MED2, kbin_eta_OVERLAP_pT_bin_MED2, kbin_eta_ENDCAP_pT_bin_MED2],
#       [kbin_eta_BARREL_pT_bin_MED3, kbin_eta_OVERLAP_pT_bin_MED3, kbin_eta_ENDCAP_pT_bin_MED3],
#   ]
#   kbin_ls_eta_INCLUSIVE_pT_MED2 = [
#       [kbin_eta_BARREL_pT_bin_MED4, kbin_eta_OVERLAP_pT_bin_MED4, kbin_eta_ENDCAP_pT_bin_MED4],
#       [kbin_eta_BARREL_pT_bin_MED5, kbin_eta_OVERLAP_pT_bin_MED5, kbin_eta_ENDCAP_pT_bin_MED5],
#       [kbin_eta_BARREL_pT_bin_MED6, kbin_eta_OVERLAP_pT_bin_MED6, kbin_eta_ENDCAP_pT_bin_MED6],
#   ]
#   
#   # HIGH pT.
#   kbin_ls_eta_INCLUSIVE_pT_HIGH1 = [
#       [kbin_eta_BARREL_pT_bin_HIGH1, kbin_eta_OVERLAP_pT_bin_HIGH1, kbin_eta_ENDCAP_pT_bin_HIGH1],
#       [kbin_eta_BARREL_pT_bin_HIGH2, kbin_eta_OVERLAP_pT_bin_HIGH2, kbin_eta_ENDCAP_pT_bin_HIGH2],
#       [kbin_eta_BARREL_pT_bin_HIGH3, kbin_eta_OVERLAP_pT_bin_HIGH3, kbin_eta_ENDCAP_pT_bin_HIGH3],
#   ]
#   kbin_ls_eta_INCLUSIVE_pT_HIGH2 = [
#       [kbin_eta_BARREL_pT_bin_HIGH4, kbin_eta_OVERLAP_pT_bin_HIGH4, kbin_eta_ENDCAP_pT_bin_HIGH4],
#       [kbin_eta_BARREL_pT_bin_HIGH5, kbin_eta_OVERLAP_pT_bin_HIGH5, kbin_eta_ENDCAP_pT_bin_HIGH5],
#       [kbin_eta_BARREL_pT_bin_HIGH6, kbin_eta_OVERLAP_pT_bin_HIGH6, kbin_eta_ENDCAP_pT_bin_HIGH6],
#   ]