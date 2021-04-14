"""Derive AdHoc pT Correction Factors

This code runs a series of scripts to derive AdHoc pT correction factors from 
a muon sample (DY, ggH, etc.)

Parameters
----------
infile_path : The root file from which muons will be pulled.
"""
from Utils_Python.Utils_Files import make_dirs #open_pkl, save_to_pkl, check_overwrite
sample_type = "MC"  # "MC", "Data"
year = "2016"
prod_mode = "DY2mu"  # "DY2mu", "DY2e", "H2mu", "H2e", "H4mu", "H4e"
outfile_prefix = "fullanalysis_test01"
infile_path = "/cmsuf/data/store/user/t2/users/ferrico/SingleBS_studies/After/DY_2018_2.root"

outdir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/"


submit_to_slurm = True

inv_m_lim = [60.0, 120.0]
eta_lim = [0.0, 2.4]
pT_lim = [5.0, 1000.0]
d0_lim = [0.0, 1000.0]  # No d0 cuts now. Apply d0 cuts when separating muons into KB2Ds.
dR_max = 100

specified_evt_beg = None # 122000000
specified_evt_end = None

verbose = 1
overwrite = 0
evts_per_batch = 2E6  # The smaller, the faster it will be processed on SLURM.
# SLURM directives.
partition = "hpg2-compute"
mem = (8, 'gb')
nodes = 4
burst = True

make_dirs(os.path.join(outdir, ""))
assert sample_type in ('MC', 'Data')

outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY/outtxt"
outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY/copies"
outpkl_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2018/DY/"

import argparse
parser = argparse.ArgumentParser(description='')
parser.add_argument('-v', action="store_true")
args = parser.parse_args()
if args.v:
    print("Debug info.")

# With `action="store_true"`, just specify the flag;
# don't provide an option, in the true spirit of flags.
parser.add_argument("-v", "--verbose", help="increase output verbosity",
                    action="store_true", default=False)
parser.add_argument('--min', dest='min_Zmass', type=float, help='min for Zmass')
parser.add_argument('--filename', dest='filename', type=str, help='output file')
parser.add_argument('--widthZ', dest='widthZ', type=float, help='Z width in MC or pdg value')
parser.add_argument('--plotBinInfo', dest='binInfo', nargs='+', help='', type=int)
parser.add_argument('--doubleCB_tail',dest='doubleCB_tail', nargs='+', help='', type=float)
args = parser.parse_args()

# Retrieve values from args object.
args.Z_width
massZErr_rel_min = args.min_relM2lErr

# This is the correct way to handle accepting multiple arguments.
# '+' == 1 or more.
# '*' == 0 or more.
# '?' == 0 or 1.
# An int is an explicit number of arguments to accept.
parser.add_argument('--nargs', nargs='+')

# To make the input integers
parser.add_argument('--nargs-int-type', nargs='+', type=int)

#------------------------------#
#--- Step 1: Skim root file ---#
#------------------------------#
print(f"Output dirs:")
print(f"  text:    {outtxt_dir}")
print(f"  copies:  {outcopies_dir}")
print(f"  pickles: {outtxt_dir}")
if submit_to_slurm:
    # Would need to feed all these into this script:
    # year = "2016"
    # prod_mode = "DY2mu"
    # outfile_prefix = "fullanalysis_test01"
    # infile_path = "/cmsuf/data/store/user/t2/users/ferrico/SingleBS_studies/After/DY_2018_2.root"

    # outtxt_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY2mu/output/"
    # outcopies_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY2mu/copies/"
    # outpkl_dir    = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2018DY2mu/pickles/"

    # submit_to_slurm = True

    # inv_m_lim = [60.0, 120.0]
    # eta_lim = [0.0, 2.4]
    # pT_lim = [5.0, 1000.0]
    # d0_lim = [0.0, 1000.0]  # No d0 cuts now. Apply d0 cuts when separating muons into KB2Ds.
    # dR_max = 100

    # specified_evt_beg = None # 122000000
    # specified_evt_end = None

    # verbose = 1
    # overwrite = 0
    # evts_per_batch = 2E6  # The smaller, the faster it will be processed on SLURM.
    # # SLURM directives.
    # partition = "hpg2-compute"
    # mem = (8, 'gb')
    # nodes = 4
    # burst = True
    python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/skim_sample_inbatch_submit.py
    trash = input("[INFO] Press 'Enter' when your SLURM jobs finish")

#-------------------------------------#
#--- Step 2: Sort muons into KB2Ds ---#
#-------------------------------------#
    # Feed in:
    # year = "2016"
    # prod_mode = "DY2mu"
    # eta_ls = equal_entry_bin_edges_eta_mod1_wholenum #eta_bins_geofitvsVXBS
    # pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum #pT_bins_geofitvsVXBS

    # d0_lim = [0.0, 1000.0]  # Only keep muons in this d0 bin.

    # n_bins_dpTOverpT = 100
    # x_lim_dpTOverpT = [-0.4, 0.4]
    # n_bins_qd0 = 100
    # x_lim_qd0 = [-0.01, 0.01]
    # verbose = 1
    # overwrite = 0

    # make_hists = 0
    # indirs = (
    # "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/pickles/",
    # "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/pickles/"
    # )
    # glob_template = "MC2016DY_skim_fullstats_nogenmatching_skim_sample_inbatch_template_copy_*.pkl"
    # outfile = "MC2016DY_skim_fullstats_nogenmatching_0p0_d0_0p01"

    #  the outtxt_dir from last script as input dir:
    # /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr
    # Automatically get 
    python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/merge_muon_collections.py


#------------------------------#
#--- Step 3: Separate KB2Ds ---#
#------------------------------#
    # Feed in:
    # year = "2016"
    # infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/MC2016DY_skim_fullstats_nogenmatching_0p01_d0_1000p0_withGeoFitcorr.pkl"
    # file_prefix = "" #"MC2016DY_fullstats_muoncoll_withkb3dbins_nogenmatching_0p0_d0_0p01"
    # overwrite = 0
    # verbose = 1
    step3_flags = ("-v", "-o", f"--year={year}", f"--infile={infile_path}")
    
    python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/save_kb2ds_separate_dicts.py 

#-----------------------------#
#--- Step 4: Analyze KB3Ds ---#
#-----------------------------#
    # year = "2016"
    # overwrite = 0
    # verbose = 0

    # iters = 2
    # fit_with_zero_interc = True
    # regions = 4
    # min_muons_per_qd0_bin = 100
    # fit_whole_range_first_iter = False  # False gives more consistent fits (with no outlier data).
    # use_data_in_xlim = 1
    # binned_fit = False
    # switch_to_binned_fit = 999999999

    # job_name_base = "MC2016DY_individKB2D_withitergaussfitsonKB3Ds_fitwithzerointerc"  # Also prefix for outfile.
    # # job_name_base = "MC2016DY_individKB2D_withitergaussfitsonKB3Ds_nogenmatching_0p0_d0_0p01"  # Also prefix for outfile.
    # # eta_ls = eta_bins_geofitvsVXBS #equal_entry_bin_edges_eta_mod1_wholenum
    # eta_ls = equal_entry_bin_edges_eta_mod1_wholenum[0:2]
    # # full_pT_ls = pT_bins_geofitvsVXBS  #[20.0, 30.0, 40.0, 50.0] #bin_edges_pT_sevenfifths_to1000GeV_wholenum
    # full_pT_ls = bin_edges_pT_sevenfifths_to1000GeV_wholenum[0:2]

    # template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_singlekb2d_itergaussfitsonkb3ds_template.py"
    # # template_script_main = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py"
    # template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_to_slurm_template.sbatch"
    # # template_script_slurm = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/RochCorrAnalyzers/roch_vs_noroch_slurm.sbatch"

    # # inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/pickles/2016/DY/MC2016DY_skim_fullstats_nogenmatching/kb2d_dictsnofitinfo/MC2016DY_fullstats_muoncoll_withkb3dbins_nogenmatching_0p0_d0_0p01__ETAPART_PTPART.pkl"
    # inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/skim_fullstats_verify/pickles/kb2d_dicts/ETAPART_PTPART.pkl"

    # # /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/output/DeriveCorrections/2016DY/individualKB2Ditergaussfits/MC2016DY_kb2d_ETARANGE_PTRANGE_error.log
    # outdir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu"
    # outdir_copies = os.path.join(outdir, f"{job_name_base}/copies")
    # outdir_pkl    = os.path.join(outdir, f"{job_name_base}/pickles")
    # outdir_txt    = os.path.join(outdir, f"{job_name_base}/output")
    # outdir_pdf    = os.path.join(outdir, f"{job_name_base}/plots")

    # # SLURM directives.
    # partition = "hpg2-compute"
    # mem = (8, 'gb')
    # nodes = 1
    # burst = True
    # time = "01:00:00"  # hh:mm:ss

    step4_flags = ()
    python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/slurm_inbatch_singlekb2d_itergaussfitsonkb3ds.py

#---------------------------#
#--- Step 5: Merge KB2Ds ---#
#---------------------------#
    # Inputs:
    # year = "2016"
    # prod_mode = "DY2mu"
    # inpkl_path_template = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc/pickles/individKB2DwithitergaussfitsonKB3Ds_fitwithzerointerc_*.pkl"
    # overwrite = 0
    # make_pTcorrdict = 1
    # outfilename_base = "muoncoll_itergaussfitsonKB3Ds"
    # outpkl_dir = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY2mu/muoncollwithfitstats"
    python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/merge_individkb2ds_and_make_pTcorrfactordict.py


# IN APPLY CORRECTION SCRIPT:

# Merge KB2Ds but don't make pT corr dict (it was already made).
# Name outfile with prefix like: "muoncoll_withKB2Dfits_beforeaftercorr"
make_pTcorrdict = 0
python merge_individkb2ds_and_make_pTcorrfactordict.py
# Take the muoncoll which gets made and pass it into next script:

# Merge MuonColls.
python /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/merge_muoncoll_withgraphsandwithbeforeaftercorr.py