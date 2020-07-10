# How to Improve the Width of the Higgs Boson Mass Resonance using the Correlation between muon pT mismeasurement and d0

## Analysis Steps

1. Use muons from Drell-Yan and J/psi samples (and potentially Higgs samples, Upsilon samples, etc.).
2. Bin muons in eta and pT.
For this analysis I used 
3. 

`HiggsMassMeasurement/d0_Studies/Plotters_vaex/make_pdf_qd0_hists_multiplesamples_etapTbins.py`

1. Find equal-entry q*d0 bins:
`HiggsMassMeasurement/d0_Studies/d0_Analyzers/find_equalentry_qd0_regions_etapTbins.py`
   - Be sure to check the number of regions (r)!
   - This will produce a `.pkl` with all the equal-entry qd0 region info. 
      * A `.csv` is also produced to validate the info by eye.
      It has the same info as the `.pkl`.

1. Do **iterated gaus fit** in DeltapT/pT dists to retrieve best
   - This can take a very long time. 
   Example: 13 minutes to produce 2 PDFs.
   Each PDF had 12 pages. 
   Each page had ~6 plots.
`HiggsMassMeasurement/d0_Studies/d0_Analyzers/findKinBins_doFits_makePDFs.py`

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