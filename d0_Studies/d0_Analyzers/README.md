# How to derive AdHoc pT corr factors

## Main idea

- Derive AdHoc corrections from muon kinematics.
- Apply pT corrections per muon to see the improvement in sigma of dpT/pT
distributions, m(4mu), and m(2mu).

---

## Overview

### Derive Corrections

1. Skim a rootfile (say Drell-Yan) to put muons into a MyMuonCollection.
1. Sort the muons into bins of (eta, pT) - call these KinBin2Ds (KB2Ds).
1. For each KB2D, sort the muons into *n* KB3Ds (eta, pT, q*d0).
1. Analyze the dpT/pT and 1/pT distributions of muons in each KB3D.
1. Graph mean(dpT/pT) vs. q*d0 for each KB2D to extract pT correction factors.

### Apply Corrections

1. Apply pT correction factors to muons in root file or MyMuonCollection.

### Check Results

1. Make DSCB fits of ggH before/after pT corrections per muon
to measure improvement of sigma_DSCB.
1. Measure improvement of muon dpT/pT in (eta, pT) bins.

---

## Detailed Workflow

When working with many muons (>20E6),
it is useful to submit your jobs to SLURM.

**NOTE:**
The workflow below shows processing times for working over a MC 2016 DY sample
which had 122M unskimmed events.
The resulting skimmed `.pkl` file has 28M muons.

### Derive Corrections

1. Skim a root file with:

   ```bash
   python skim_sample_inbatch_submit.py
   ```

   - Which calls: `skim_sample_inbatch_template.py`.
   - This will produce many pickled `MyMuonCollections`.
   - Processing time: from ~`00:10:00 (hh:mm:ss)`
      - (Requesting 2M events per batch.)
      - Output `.pkl` files will be ~1 MB / 18K events.

1. Merge the `.pkl` files with:
   
   ```bash
   python merge_muon_collections.py
   ```

   - Stores inclusive hists of muon kinematics.
   - Sorts muons into KB2Ds based on given eta and pT bins.
   - Overwrites `MyMuonCollection.muon_ls`. Passes MyMuons to KB2Ds.
   - Processing times: from `00:02:42` to `00:21:04`

1. Save KB2Ds into individual pickled dicts with:

   ```bash
   # First start a dev session to access more RAM.
   srun --partition=bigmem --mem=128gb --ntasks=1 --cpus-per-task=8 --time=08:00:00 --pty bash -i
   # Then run the script.
   python save_kb2ds_separate_dicts.py
   ```
   
   - Processing time: `00:03:18`

1. Do recursive Gaussian fits on the dpT/pT dists of each KB3D with:

   ```bash
   python slurm_inbatch_singlekb2d_itergaussfitsonkb3ds.py
   ```

   - Which calls:
      - `submit_singlekb2d_itergaussfitsonkb3ds_template.py`
      - `submit_to_slurm_template.sbatch`
   - For each KB2D, this step produces:
      - dpT/pT vs. q*d0 graph
      - dpT/pT * (1/\<pT\>) vs. q*d0 graph
      - dpT/pT * (\<1/pT\>) vs. q*d0 graph
      - dpT/pT * mu_DSCB(\<1/pT\>) vs. q*d0 graph
   - Processing time (depends on size of KB3Ds): ~`00:24:00` max
   - Jobs start immediately on `hpg2-compute` using, **1 node**, no burst.
   - **TO DO:**
      - Implement SlurmManager class to eliminate `.sbatch` script.
      - Confirm that fits were error-free.

1. Merge KB2Ds (which contain analyzed KB3Ds) into a MyMuonCollection with:

   ```bash
   python merge_individkb2ds_and_make_pTcorrfactordict.py
   ```

   - Saves the following as individual '`.pkl`'s:
      - pT correction factor dict.
      - MyMuonCollection with analyzed KB2Ds and KB3Ds.
   - Processing time: `00:00:19`
      - If you keep each `muon_ls` then ~`00:04:26`.
   - NOTE: These KB2Ds don't yet have before/after pT corrections.

### Apply Corrections

There are a few different scripts to apply pT corrections.

1. Do unbinned IGFs on each KB2D's dpT/pT and dpTcorr/pT distribution:

   ```bash
   python singlekb2ditergaussfit_submit2slurm.py
   ```

   - Which calls:
      - `singlekb2ditergaussfit_slurm_template.py`
   - Processing time for 5 **unbinned** IGFs: ~`00:15:00`
   
1. Merge KB2Ds in MyMuonCollection without creating a new pT corr dict:

   ```bash
   python merge_individkb2ds_and_make_pTcorrfactordict.py
   ```

   - Set `make_pTcorrdict = 0`.

   Then merge this muon_coll with the other muon_coll you made from
   **Derive Corrections**:

   ```bash
   python merge_muoncoll_withgraphsandwithbeforeaftercorr.py
   ```

<!-- 1. Correct muons stored in a MyMuonCollection with:

   ```bash
   python apply_pTcorrfactors_to_H4mu_sample.py
   # or
   python applypTcorrfactors_to_muoncoll.py
   ```

   - Processing time:

1. Then merge KB2DS

   ```bash
   python merge_individkb2ds_and_make_pTcorrfactordict.py
   ```
  -->

### Check Results

1. Check the improvement on sigma of m(4mu) by performing an unbinned DSCB fit
on the m(4mu) dist and plot it with:

   ```bash
   root -l DSCB_fit_m4mu_beforeafterpTcorr.C
   ```

   - Processing time: 

Finally, print and plot the results:

   ```bash
   python make_finalplots_in_muoncoll.py
   ```

---

## Plotting

Make all plots, including KB3D fits and inclusive kinematics plots from
MyMuonCollection with:

```bash
d0_Studies/d0_Analyzers/plot_all_KB3Diterfits.py
```

### 2D Tables

Make 2D tables (eta vs. pT) of sigma(dpT/pT) and improvements on
sigma(dpT/pT):

   ```bash
   d0_Studies/Plotters_Python/make_2D_plot_sigma_table_dpToverpT.py
   ```

---

## Running locally (instead of submitting to SLURM)

1. Same as 1. above, except use `skim_sample.py`.

---

## Locations of scripts

```bash
d0_Studies/d0_Analyzers/
d0_Studies/Plotters_ROOT/
```