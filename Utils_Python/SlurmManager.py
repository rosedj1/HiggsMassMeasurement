from Utils_Python.Utils_Files import check_overwrite

class SlurmManager:
    """A class to handle the details of working with SLURM.
    FIXME: Not implemented yet.
    """
    def __init__(self): #, path_main, path_sbatch):
        self.file_copies_ls = []
        pass

    def copy_file(self, src, dst):
        """Copy `src` file to `dst`. Store the paths."""
        try:
            shutil.copyfile(src, dst)
            self.file_copies_ls.append(dst)
        except:
            print(f"[WARNING] Could not copy\n{src}to\n{dst}")

    def write_copies_to_file(self, outf):
        """Save the file paths of file_copies_ls to file `outf`."""
        # with open(outf) as f:
        #     f.
        pass

class SLURMSubmitter:
    """Generate a new SLURM script."""

    def __init__(self, verbose=False):
        self.directive_dct = {}
        self.prescript_text = """
            pwd; hostname; date
            source ~/.bash_profile
            cd /blue/avery/rosedj1/HiggsMassMeasurement/
            conda activate my_root_env
            source setup.sh
            echo 'Packages loaded.'
            cd /blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/d0_Analyzers/
            echo 'Starting script...'
        """
        self.postscript_text = """
            echo 'Script finished!'
        """
        self.verbose = verbose

    def prep_directives(self, job_name, output_txt, email="rosedj1@ufl.edu", time="08:00:00", acct="avery", burst=False, mem=1, partition="hpg2-compute", nodes=1):
        """Store SLURM directives.
        
        Parameters
        ----------
        job_name : str
            The name of the job displayed when doing `squeue -u <username>`.
        output_txt : str
            The full path to the stdout and stderr files.
        email : str
            Email address to which SLURM sends alerts.
        time : str
            Max length of job, using this format: hh:mm:ss
        acct : str
            Use resources on this group account.
        burst : bool
            If True, use 9x more resources from your group, if available.
        mem : str
            RAM to use (GB).
        partition : str
            Send jobs to this partition on HPG.
            Options: 'hpg2-compute', 'bigmem'
        nodes : str
            Number of nodes (servers, computers) to use.
        """
        self.directive_dct["job-name"] = job_name
        self.directive_dct["output"] = output_txt
        self.directive_dct["error"] = f"{output_txt.split('.')[0]}_error.log"
        self.directive_dct["mail-type"] = "ALL"
        self.directive_dct["mail-user"] = email
        self.directive_dct["time"] = time
        self.directive_dct["account"] = acct
        self.directive_dct["qos"] = "avery-b" if burst else "avery"
        self.directive_dct["mem"] = mem
        self.directive_dct["partition"] = partition
        self.directive_dct["nodes"] = nodes

    def write_directives(self, f):
        """Write stored directives to file object `f`."""
        f.write("#!/bin/bash\n")
        for dctv, val in self.directive_dct.items():
            f.write(f"#SLURM --{dctv}={val}\n")
        f.write("\n")

    def write_cmds(self, cmdtup, f):
        """Write commands in `cmdtup` to file object `f`.

        cmdtup : tuple of str
            Each element corresponds to one line of commands.
        """
        if isinstance(cmdtup, str):
            # Put into a tuple so the string won't be split up.
            cmdtup = (cmdtup, "")
        for cmd in cmdtup:
            self.write_text(cmd, f, newline=True)

    def write_text(self, text, f, newline=True):
        """Write `text` to file object `f`, followed by newline."""
        nl = "\n" if newline else ""
        f.write(f"{text}{nl}")

    def make_slurm_script(self, slurm_outpath, cmdtup, overwrite=False):
        """
        Write SLURM script to `slurm_outpath` that executes commands in `cmdtup`.
        Return 0 on success.
        
        Parameters
        ----------
        slurm_outpath : str
            The absolute path of the output SLURM script.
        cmdtup : tuple of strings
            A tuple of commands to execute.
            NOTE: Use triple double-quotes: "
            Each element corresponds to a command to be executed on a single
            line by the interpreter.
            Example:
                say you want to do `ls -l` followed by `python code.py`
            Then do: 
                cmdtup = ("ls -l", "python code.py")
        """
        try:
            assert all(d is not None for d in self.directive_dct.values())
        except AssertionError:
            # User hasn't specified all directives yet.
            print("[FAIL] You need to specify all directives.")
            return 1
        try:
            assert len(cmdtup) > 0
        except AssertionError:
            print("[FAIL] You need to specify some commands.")
            return 1
        check_overwrite(slurm_outpath, overwrite)
        # output = shell_cmd(f"touch {slurm_outpath}")  # Doesn't delete file if file exists.
        with open(slurm_outpath, "w") as f:
            if self.verbose:
                print(f"Writing directives to SLURM script.")
            self.write_directives(f)
            if self.verbose:
                print(f"Writing pre-script instructions to SLURM script.")
            self.write_text(self.prescript_text, f)
            if self.verbose:
                print(f"Writing commands to SLURM script.")
            self.write_cmds(cmdtup, f)
            if self.verbose:
                print(f"Writing post-script instructions to SLURM script.")
            self.write_text(self.postscript_text, f)
        if self.verbose:
            print(f"SLURM script successfully written:\n{slurm_outpath}")
        return 0
            
    def submit_script(self):
        """"""
        pass