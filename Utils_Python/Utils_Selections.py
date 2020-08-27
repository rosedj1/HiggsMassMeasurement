class Selector():
    """
    A class to organize the selections for different samples. 
    """
    def __init__(self, sample_name, verbose=False):
        """
        sample_name : str
            "Jpsi", "DY", "Upsilon", "Higgs"
        verbose : bool, optional
            Show debug info.
        """
#         pass_overall=None
#         pass_overall : bool  
#             If False, failed at least 1 selection. 
#             If True, has passed all selections.
#             If None, has not had any selections applied.
#         self.pass_overall = pass_overall
            
        self.sample_name = sample_name
        self.verbose = verbose
        
#         self.passed_dR_cut = False
#         self.passed_m2l_cut = False
#         self.passed_pT_cut = False
#         self.passed_eta_cut = False
        
        self.cut_dict = {
            "Jpsi"    : {"m2l_cut" : [2.9, 3.3], "dR_cut" : 0.005},  # Filippo: 0.005
            "DY"      : {"m2l_cut" : [60, 120],  "dR_cut" : 0.002},  # Filippo: 0.002
            "Upsilon" : {"m2l_cut" : [8.5, 11],  "dR_cut" : 0.005},  # Filippo: 0.005
            }
        
#     def check_passed_all_selections(self):
#         self.passed_all = all( (self.passed_dR_cut, 
#                                 self.passed_m2l_cut,
#                                 self.passed_pT_cut,
#                                 self.passed_eta_cut)
#                              )
    
    def pass_selection_val(self, val, val_min=None, val_max=None):
        """
        Return True if val is within val_min, val_max.
        Otherwise, returns False.
        
        Useful for applying pT cuts but is used in other methods. 
        """
        if (val_min is not None) and (val_max is not None):
            return True if (val > val_min) and (val < val_max) else False
        elif (val_min is None) and (val_max is not None):
            # Check against val_max only. 
            return True if (val < val_max) else False
        elif (val_min is not None) and (val_max is None):
            # Check against val_min only. 
            return True if (val > val_min) else False
        else:
            msg = "[WARNING] You performed a cut, but didn't specify any bounds."
            raise RuntimeWarning(msg)
        
    def pass_selection_abseta(self, eta, eta_min=None, eta_max=None):
        """
        Return True if abs(eta) is within eta_min, eta_max bounds.
        Otherwise, returns False.
        """
        return self.pass_selection_val(abs(eta), eta_min, eta_max)

    def pass_selection_m2l(self, m2l):
        """Return True if m2l passes selections, based on sample_name."""
        m2l_cut_ls = self.cut_dict[self.sample_name]["m2l_cut"]
        return self.pass_selection_val(val=m2l, val_min=m2l_cut_ls[0], val_max=m2l_cut_ls[1])
    
    def pass_selection_dR(self, dR):
        """Return True if dR passes selections, based on sample_name."""
        dR_max = self.cut_dict[self.sample_name]["dR_cut"]
        return self.pass_selection_val(val=dR, val_max=dR_max)