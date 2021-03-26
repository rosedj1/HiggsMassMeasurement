#!/bin/bash

year="2018"

base_dir="/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies"  # No '/'.
inpkl_pTcorrfactors="/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/pTCorrFactors/AdHoc/AdHocpTcorrfact_derivedfromMC2018DYJpsi_13etabins12pTbins_0p0eta2p4_5p0pT1000p0.pkl"
outpath_root_beforeaftercorr="/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromAdHocfactors_fullstats_zerointerc_new100muonsperregion.root"
outpath_pdf_DSCB="/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/plots/applypTcorrplots/CorrFromMC/test/MC2018ggH_m4mu_DSCBfit_beforeaftercorr_AdHocfactors_zerointerc_new100muonsperregion.pdf"

# Apply pT correction factors. Make root file before/after corr.
sed ...
python ${base_dir}/d0_Analyzers/apply_pTcorrfactors_to_H4mu_sample.py

# DSCB plotter: check improvement of sigma(DSCB) before/after corrections.
cd ${base_dir}/Plotters_ROOT/
cp DSCB_fit_m4mu_beforeafterpTcorr_template.C DSCB_fit_m4mu_beforeafterpTcorr.C
sed -i "s/YEAR/${year}/g" DSCB_fit_m4mu_beforeafterpTcorr.C
sed -i "s/INPKL_PTCORRFACTORS/${inpkl_pTcorrfactors}/g" DSCB_fit_m4mu_beforeafterpTcorr.C
sed -i "s/ROOTPATH/${outpath_root_beforeaftercorr}/g" DSCB_fit_m4mu_beforeafterpTcorr.C
sed -i "s/OUTPATH_PDF/${outpath_pdf_DSCB}/g" DSCB_fit_m4mu_beforeafterpTcorr.C
# root -l DSCB_fit_m4mu_beforeafterpTcorr.C

# 
# Make 2D tables.