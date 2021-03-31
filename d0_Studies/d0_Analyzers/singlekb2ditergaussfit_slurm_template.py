from Utils_Python.Utils_Files import open_pkl, save_to_pkl

overwrite = OVERWRITE
inpkl_path_kb2d = "INPKL_PATH"
outpkl_path_kb2d = "OUTPKL_PATH"
delete_kb2d_muon_ls = DELETE_KB2D_MUON_LS
switch_to_binned_fit = SWITCH2BINNED

kb2d = open_pkl(inpkl_path_kb2d)
# Without pT corr.
kb2d.do_itergausfit(bins_dpTOverpT=1000, bins_qd0=100, 
                    x_lim_dpTOverpT=[-0.5, 0.5], x_lim_qd0=[-0.01, 0.01], 
                    fit_whole_range_first_iter=True, 
                    iters=5, num_sigmas=2.5, marker_color=None, line_color=None, 
                    switch_to_binned_fit=switch_to_binned_fit, verbose=True, alarm_level="warning", 
                    use_mu_pT_corr=False, only_draw_last=False) 
# WITH pT corr.
kb2d.do_itergausfit(bins_dpTOverpT=1000, bins_qd0=100, 
                    x_lim_dpTOverpT=[-0.5, 0.5], x_lim_qd0=[-0.01, 0.01], 
                    fit_whole_range_first_iter=True, 
                    iters=5, num_sigmas=2.5, marker_color=None, line_color=None, 
                    switch_to_binned_fit=switch_to_binned_fit, verbose=True, alarm_level="warning", 
                    use_mu_pT_corr=True, only_draw_last=False)
kb2d.overwrite_muon_info(delete_all=True)
save_to_pkl(kb2d, outpkl_path_kb2d, overwrite=overwrite)