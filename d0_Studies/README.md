# Ad hoc pT corrections to improve muon momentum resolution

## Analysis Steps

1. Derive pT corrections by studying muons in different regions of phase space
(|eta|, pT, q*d0).
2. Apply the linear pT corrections per muon to see the improvement on dpT/pT
and m4mu.

### Step 1: Derive pT corrections

Do:

```bash
python HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample.py
```

- Make sure all the 'User Parameters' are correct for your purpose.

For very big jobs... This script can take a _very_ long time.
This is due to analyzing something like 12 eta
bins, 13 pT bins, ~10 q*d0 bins and performing ~5 iterated Gaussian fits each
(~7800 unbinned fits!).
To get around this, you can break the job up in several jobs, where each job
will work on a single eta bin:

```bash
# So modify and run this script:
HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py
# which will call this script, so check the parameters beforehand.
HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py
```

Your pT correction factors should be saved in `.pkl` files.

- If you have individual `.pkl` files (say one for each eta bin),
then you can combine them into a single **combined pT corr factor** `.pkl` by:

```bash
python HiggsMassMeasurement/d0_Studies/d0_Analyzers/merge_pT_corr_dcts.py
```

### Step 2: Apply pT corrections

Run the script:

```bash
python HiggsMassMeasurement/d0_Studies/Plotters/plot_m4mu_kinematics_withpTcorr_MC2017ggH.py
```

- This script can make a root file which contains m4mu and m4mu_corr vals.
  - Then a DSCB script can be run: 

---

## Aternative method: Deriving and applying corrections from DY and J/psi samples using Vaex

1. Use muons (from Drell-Yan, J/psi, Upsilon samples, etc.) and bin the muons in eta and pT.
2. For this analysis I used:

`HiggsMassMeasurement/d0_Studies/Plotters_vaex/make_pdf_qd0_hists_multiplesamples_etapTbins.py`

1. Find equal-entry q*d0 bins:
`HiggsMassMeasurement/d0_Studies/d0_Analyzers/find_equalentry_qd0_regions_etapTbins.py`
   - Be sure to check the number of regions (r)!
   - This will produce a `.pkl` with all the equal-entry qd0 region info. 
      * A `.csv` is also produced to validate the info by eye.
      It has the same info as the `.pkl`.

1. Do **iterated Gaussian fits** in dpT/pT dists to retrieve best fit params.
`HiggsMassMeasurement/d0_Studies/d0_Analyzers/findKinBins_doFits_makePDFs.py`
   - This can take a very long time. 
      - Example: 13 minutes to produce 2 PDFs.
   Each PDF had 12 pages. 
   Each page had ~6 plots.

1. Make dpT/pT vs. qd0 graphs.

`python plot_dpToverpT_vs_d0q_graph_binning_eta.py`

- Produces 1 PDF: 
each page shows 1 graph of many lines. 
All lines have constant eta bin and all pT ranges. has 1 graph of 

1. Retrieve pT correction factors and store in pickled dict.

`python /ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plotters/makegraph_dpToverpTvsqd0_getpTcorrfactors.py`

1. Apply pT corrections to muons in Higgs sample. 
Fit a double-sided Crystal Ball (DSCB) function
   
`root -l DSCB_fit.C`

### Using the pT correction factors on dimuon samples

1. First run the script: `d0_Studies/Plotters/applypTcorr_dpToverpT_plots_combinedsamples.py`

   - Produces a pickled dictionary :
   `d0_Studies/pkl_and_csv/Testing/20200709_MC2017_combinesamples_applycorr_fullstats_2p00sigmas__0p0_eta_2p4__5p0_pT_1000p0_GeV.pkl`
   - Works on only 1 year at a time. 