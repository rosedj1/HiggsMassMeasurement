"""FIXME: Incomplete. Still deciding if this file is necessary.

Functions to create and manipulate ROOT.TTree objects in Python."""

def add_branch_pTcorr_fromadhocmethod(tree, pT_corr_factor_dict, outpath_root):
    """
    Add a branch to a clone of a TTree.
    The branch is filled with pT corrected values, using ad hoc correction
    factors from pT_corr_factor_dict.
    The cloned tree is written to a new '.root' file (outpath_root).
    """
    outfile_root = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/"
    newfile = ROOT.TFile(outfile_path, "recreate")
    treeclone = tree.CloneTree(0)
    ptr = array('f', [0.])
    treeclone.Branch("pTcorr_adhoc", ptr, "m4mu/F")


def write_m4muinfo_to_rootfile(self, outpath_rootfile, overwrite=False):#, write_geofit_vals=False, ):
    """Write m4mu and m4mu_corr info to root file for DSCB fits, etc."""
    # New file, TTree, and TH1Fs.
    check_overwrite(outpath_rootfile, overwrite)
    print(f"  Writing m4mu and m4mu_corr vals to root file:")
    print(f"  {outpath_rootfile}")
    outf = r.TFile(outpath_rootfile, "recreate")
    newtree = r.TTree("tree", "tree_m4mu_vals")

    # Set pointers for values to be stored in root file.
    ptr_m4mu = array('f', [0.])
    ptr_m4mu_corr = array('f', [0.])
    newtree.Branch("m4mu", ptr_m4mu, "m4mu/F")
    newtree.Branch("m4mu_corr", ptr_m4mu_corr, "m4mu_corr/F")

    # if write_geofit_vals:
    #     ptr_m4mu_geofit_corr = array('f', [0.])
    #     newtree.Branch("m4mu_geofit_corr", ptr_m4mu_geofit_corr, "m4mu_geofit_corr/F")
    
    # Loop over stored values inside MyMuonCollection.
    print(f"...Filling m4mu and m4mu_corr values...")
    start = time.perf_counter()
    assert len(self.m4mu_ls) == len(self.m4mu_corr_ls) != 0
    # if write_geofit_vals:
    #     assert len(self.m4mu_ls) == len(self.m4mu_corr_geofit_ls)

    # if write_geofit_vals:
    #     for ct,(m4mu,m4mu_corr,m4mu_corr_geofit) in enumerate(zip(self.m4mu_ls, self.m4mu_corr_ls, self.m4mu_corr_geofit_ls)):
    #         m4mu_diff = m4mu_corr - m4mu
    #         m4mu_diff_geofit = m4mu_corr_geofit - m4mu
    #         rel_diff = m4mu_diff / m4mu
    #         rel_diff_geofit = m4mu_diff_geofit / m4mu

    #         if (abs(rel_diff) > 0.05) or (abs(rel_diff_geofit) > 0.05):
    #             print(f"event {ct}:")
    #             print(f"  m4mu={m4mu}, m4mu_corr={m4mu_corr}, rel_diff={rel_diff}, rel_diff_geofit={rel_diff_geofit}")

    #         ptr_m4mu[0] = m4mu
    #         ptr_m4mu_corr[0] = m4mu_corr
    #         ptr_m4mu_geofit_corr[0] = m4mu_corr_geofit
    #         newtree.Fill()
    # else:
    for ct,(m4mu,m4mu_corr) in enumerate(zip(self.m4mu_ls, self.m4mu_corr_ls)):
        m4mu_diff = m4mu_corr - m4mu
        rel_diff = m4mu_diff / m4mu

        if abs(rel_diff) > 0.05:
            print(f"event {ct}:")
            print(f"  m4mu={m4mu}, m4mu_corr={m4mu_corr}, rel_diff={rel_diff}")

        ptr_m4mu[0] = m4mu
        ptr_m4mu_corr[0] = m4mu_corr
        newtree.Fill()
    print(f"  Number of good events found: {len(self.m4mu_ls)}")
    end = time.perf_counter()
    print(f"  (Took {end - start:.4f} seconds.)")

    # Save tree and close file.
    newtree.Write()
    outf.Close()