"""SLURM Script Duplicator and Submitter

This code makes a copy of a script that you wish to run on SLURM.
It is useful because it makes many copies of the main script and
of the SLURM submissions script, each copy only differing
from the next by the eta bin used. 

User can put in a whole list of eta values,
and each bin will be considered:

Suppose: full_eta_ls = [0.2, 0.4, 0.6, 0.8]
Then this code will produce a copy of the main script
and of the SLURM submission script using eta_bin = [0.2, 0.4].
The next copy will use eta_bin = [0.4, 0.6], etc.
"""
import subprocess
import shutil
import os

#-----------------------#
#--- User Parameters ---#
#-----------------------#
template_main_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/derive_pTcorrfactors_from_ggH_sample_template.py"
template_slurm_script = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/submit_derive_bigmem_template.sbatch"

full_eta_ls = [0.00, 0.20, 0.40, 0.60, 0.80, 1.00, 1.25, 1.50, 1.75, 2.00, 2.10, 2.20, 2.30, 2.40]
# full_eta_ls = [1.5, 1.75, 2.0]

delete_copies = False
#------------------------#
#--- Script functions ---#
#------------------------#
def make_eta_name(eta_ls):
    """Parse a 2-element list of eta values into a title string."""
    eta_str = f"{eta_ls[0]}eta{eta_ls[1]}"
    eta_str = eta_str.replace(".", "p")
    return eta_str

def replace_value(old, new, script):
    """Use the `sed` command to replace `old` with `new` in `script`"""
    cmd = ["sed", "-i", f"s|{old}|{new}|g", script]
    output = subprocess.run(cmd)

def make_new_files(eta_ls, template_main_script, template_slurm_script):
    """Use eta_ls in string form to name files, make copies, and return the new names."""
    eta_name = make_eta_name(eta_ls)
    suffix = f"_{eta_name}.py"
    new_main = template_main_script.replace(".py", f"_copy_{eta_name}.py")
    new_slurm = template_slurm_script.replace(".sbatch", f"_copy_{eta_name}.sbatch")

    # Make copies.
    shutil.copyfile(template_main_script, new_main)
    shutil.copyfile(template_slurm_script, new_slurm)

    return (new_main, new_slurm)

def submit_slurm_script(main, slurm, eta_ls):
    """Replace values in scripts corresponding to eta_la and submit new SLURM script.

    NOTE: Works for a 2-element eta_ls: [eta_min, eta_max]
    """
    eta_name = make_eta_name(eta_ls)
    # Replace phrases in scripts.
    replace_value("REPLACE_ETA_LS", eta_ls, main)
    replace_value("REPLACE_ETA_NAME", eta_name, slurm)
    replace_value("REPLACE_NEW_FILE", main, slurm)
    # Submit the code to SLURM.
    output = subprocess.run(["sbatch", slurm])

if __name__ == "__main__":
    for eta_min,eta_max in zip(full_eta_ls[:-1], full_eta_ls[1:]):
        eta_ls = [eta_min, eta_max]
        print(f"[INFO] Making copies of scripts and submitting SLURM script for eta range: {eta_ls}")
        new_main, new_slurm = make_new_files(eta_ls, template_main_script, template_slurm_script)
        submit_slurm_script(new_main, new_slurm, eta_ls)
        if delete_copies:
            print(f"[INFO] Removing copies of scripts.")
            os.remove(new_main)
            os.remove(new_slurm)
