# Equal-entry bins: 
# Determined by counting total number of entries (N) in dataset,
# ordering values of (eta, pT or q*d0) from least to greatest, 
# and determining bin edge where N/K is the same for all K regions.
# NOTE:
#--- Anything labelled with "_mod" is not guaranteed to 
#--- have equal entries across all bins.

# eta bins
equal_entry_bin_edges_eta = [0.0, 0.2, 0.395, 0.6, 0.805, 1.02, 1.25, 1.49, 1.75, 2.05, 2.4]
equal_entry_bin_edges_eta_mod1 = [0.0, 0.2, 0.395, 0.6, 0.805, 1.02, 1.25, 1.49, 1.75, 2.0, 2.1, 2.2, 2.3, 2.4]
equal_entry_bin_edges_eta_mod1_wholenum = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
equal_entry_bin_edges_eta_mod2 = [0.0, 0.2, 0.395, 0.6, 0.805, 1.02, 1.25, 1.49, 1.75, 2.0, 2.2, 2.4]
bin_edges_eta_outer_strip_layer = [1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9]

# pT bins
equal_entry_bin_edges_pT  = [5.0, 23.01, 29.06, 33.45, 37.05, 40.18, 43.04, 45.99, 50.46, 60.77, 200.0]
equal_entry_bin_edges_pT_mod1 = [5.0, 23.01, 29.06, 33.45, 37.05, 40.18, 43.04, 45.99, 50.46, 60.77, 200.0]
equal_entry_bin_edges_pT_mod2 = [5.0, 14, 17, 20, 25, 35, 40, 45, 50, 55, 60, 85, 200.0]
equal_entry_bin_edges_pT_sevenfifths = [5.0, 7.0, 9.8, 13.720000000000002, 19.208000000000006, 26.89120000000001, 37.647680000000015, 52.706752000000016, 73.78945280000002, 103.30523392000002]
equal_entry_bin_edges_pT_sevenfifths_mod = [5.0, 7.0, 9.8, 13.720000000000002, 19.208000000000006, 26.89120000000001, 37.647680000000015, 52.706752000000016, 73.78945280000002, 103.30523392000002]
equal_entry_bin_edges_pT_sevenfifths_to1000GeV = [5.0, 7.0, 9.8, 13.720000000000002, 19.208000000000006, 26.89120000000001, 37.647680000000015, 52.706752000000016, 73.78945280000002, 103.30523392000002, 144.627, 202.478, 1000]
bin_edges_pT_sevenfifths_to1000GeV_wholenum = [5.0, 7.0, 10.0, 14.0, 20.0, 27.0, 38.0, 50.0, 75.0, 100.0, 150.0, 200.0, 1000.0]

# q*d0 bins
equal_entry_bin_edges_qd0 = [-0.20014, -0.00257, -0.00166, -0.00103, -0.0005, -1e-05, 0.00048, 0.001, 0.00163, 0.00254, 0.19913]
binedges_qd0_equalentry_smallest_qd0RMS = [-0.01034, -0.00153, -0.00098, -0.00059, -0.00027, 3e-05, 0.00032, 0.00064, 0.00103, 0.00158, 0.0102]
binedges_qd0_equalentry_smallest_qd0RMS_clipped = [-0.005, -0.00153, -0.00098, -0.00059, -0.00027, 3e-05, 0.00032, 0.00064, 0.00103, 0.00158, 0.005]
binedges_qd0_tracker_res = [-0.005, -0.004, -0.003, -0.002, -0.001, 0, 0.001, 0.002, 0.003, 0.004, 0.005]
binedges_qd0_tracker_res_mod = [-0.012, -0.004, -0.003, -0.002, -0.001, 0, 0.001, 0.002, 0.003, 0.004, 0.012]