# Rochester Correction (RC) Study

## Comparing RCed muons to non-RCed muons

1. Create KinBin2Ds (eta, pT) and kinematic distributions of muons with and without RC:

   * You need 2 root files: one _with_ RC and the other _without_ RC.
   * Both files must contain exactly the same events.
   * Check the parameters inside of this script and then run it:

    ```python
    d0_Studies/RochCorrAnalyzers/roch_vs_noroch_kb2dmaker_inclusivehistplotter.py
    ```

2. Perform iterated Gaussian fits on each KB2D:

   * Get the `.pkl` path from step 1.
   * If you are doing some testing (a few iterations, a couple of eta ranges),
   you can run the script locally using (a).
   * Otherwise, submit the script to SLURM using (b).

   ```python
   # (a)
   d0_Studies/RochCorrAnalyzers/python roch_vs_noroch_itergausfit_local.py

   # (b)
   # First modify: d0_Studies/RochCorrAnalyzers/roch_vs_noroch_slurm.sbatch
   d0_Studies/d0_Analyzers/submit_to_slurm_inbatch.py
   ```

   **NOTE:** If you ran your jobs on SLURM, then the `.pkl` files will all be
   separate from one another.
   Use this script to merge the pickled dicts together:

   ```python
   d0_Studies/d0_Analyzers/merge_pkld_dcts.py
   ```

3. Make 2D tables of results:

    ```python
    d0_Studies/Plotters_Python/plotter_2D.py
    ```