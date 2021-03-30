# How to derive AdHoc pT corr factors

## Main idea

1. Skim a rootfile (say Drell-Yan) to put muons into a MyMuonCollection.
1. Sort the muons into bins of (eta, pT) - call these KinBin2Ds (KB2Ds).
1. For each KB2D, sort the muons into *n* KB3Ds (eta, pT, q*d0).
1. Analyze the dpT/pT and 1/pT distributions of muons in each KB3D.
1. Graph mean(dpT/pT) vs. q*d0 for each KB2D to extract pT correction factors.
1. Apply pT correction factors to say ggH sample.
1. Make DSCB fits of ggH before/after pT corrections per muon
to see improvement in sigma_DSCB.

---

## Workflow (submitting jobs to SLURM)

When working with many muons (>20E6),
it is useful to submit your jobs to SLURM.

**NOTE:**
The workflow below shows processing times for working over a MC 2016 DY sample
which had 122M unskimmed events.
The resulting skimmed `.pkl` file has 28M muons.

1. Skim a root file with:

   ```bash
   python skim_sample_inbatch_submit.py
   ```

   - Which calls: `skim_sample_inbatch_template.py`.
   - This will produce many pickled `MyMuonCollections`.
   - Processing time: from ~`00:10:00 (hh:mm:ss)`
      - (Requesting 2M events per batch.)

1. Merge the `.pkl` files with:
   
   ```bash
   python merge_muon_collections.py
   ```

   - Makes inclusive hists of muon kinematics.
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

1. Do recursive Gaussian fits on the dpT/pT and 1/pT dists of each KB3D with:

   ```bash
   python slurm_inbatch_singlekb2d_itergaussfits.py
   ```

   - Which calls:
      - `submit_singlekb2d_itergaussfits_template.py`
      - `submit_to_slurm_template.sbatch`
   - For each KB2D, this step produces:
      - dpT/pT vs. q*d0 graph
      - dpT/pT * (1/\<pT\>) vs. q*d0 graph
      - dpT/pT * (\<1/pT\>) vs. q*d0 graph
      - dpT/pT * mu_DSCB(\<1/pT\>) vs. q*d0 graph
   - **TO DO:** Implement SlurmManager class to eliminate `.sbatch` script.
   - Processing time (depends on size of KB3Ds): ~`00:24:00` max

1. Merge KB2Ds (which contain analyzed KB3Ds) into a MyMuonCollection with:

   ```bash
   python merge_individkb2ds_and_make_pTcorrfactordict.py
   ```

   - Saves the following as individual '`.pkl`'s:
      - pT correction factor dict.
      - MyMuonCollection with analyzed KB2Ds and KB3Ds.
   - Processing time: `00:00:19`
      - If you keep each `muon_ls` then ~`00:04:26`.

1. Apply the pT correction factors to muons stored in a MyMuonCollection with:

   ```bash
   python apply_pTcorrfactors_to_H4mu_sample.py
   ```

   - Processing time: 

1. Do an unbinned DSCB fit of the m(4mu) dist and plot it with:

   ```bash
   root DSCB_fit_m4mu_beforeafterpTcorr.C
   ```

   - Processing time: 

---

## Apply pT Corr Factors


```bash
python apply_pTcorrfactors_to_muoncoll_submit2slurm.py
```

   - Which calls:
      - `apply_pTcorrfactors_to_muoncoll_slurm_template.py`

      
Split and save kb2ds

Submit KB2Ds to SLURM.

Then merge KB2DS.
python merge_individkb2ds_and_make_pTcorrfactordict.py
---

## Running locally (instead of submitting to SLURM)

1. Same as 1. above, except use `skim_sample.py`.

---

## Locations of scripts

```bash
d0_Studies/d0_Analyzers/
d0_Studies/Plotters_ROOT/
```

---

## Plotting

Make all plots, including KB3D fits and inclusive kinematics plots from
MyMuonCollection with:

```d0_Studies/d0_Analyzers/plot_all_KB3Diterfits.py```