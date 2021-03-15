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
    """FIXME: Class not complete yet.

    The idea is to have this class create a SLURM script
    and submit it for the user.
    """
    def __init__(self):
        self.directive_dct = {}

    def prep_directives(self, job_name, output_txt, output_err, email, time="08:00:00", acct="avery", burst=False, mem=1, partition="hpg2-compute", nodes=1):
        """Store SLURM directives."""
        self.directive_dct["job-name"] = job_name
        self.directive_dct["output"] = output_txt if ".log" in output_txt else f"{output_txt}.log"
        self.directive_dct["error"] = job_%j_error.log
        self.directive_dct["mail-type"] = ALL
        self.directive_dct["mail-user"] = rosedj1@ufl.edu
        self.directive_dct["time"] = 12:00:00
        self.directive_dct["account"] = avery
        self.directive_dct["qos"] = avery
        self.directive_dct["mem"] = 4gb
        self.directive_dct["partition"] = partition
        self.directive_dct["nodes"] = nodes

    def write_directives(self, f):
        """Write stored directives to file `f`."""
        f.write("#!/bin/bash\n")
        for dctv, val in self.directive_dct.items():
            f.write(f"#SLURM --{dctv}={val}\n")
        f.write("\n")

    def make_slurm_script(self, outpath, cmdtup, overwrite=False):
        """Write SLURM script to outpath that executes commands in `cmdtup`.
        
        Parameters
        ----------
        outpath : str
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
        assert all(d is not None for d in self.directive_dct.values())
        check_overwrite(outpath, overwrite)
        output = shell_cmd(f"touch {outpath}")  # Doesn't delete file if file exists.
        with open() as f:
            f.write()
            for cmd in cmdtup:
                f.write(f" ")
            

    import subprocess
    cmd = ["sed", "-i", f"s|{old}|{new}|g", script]
        touch out_script_path

    def submit_script(self):
        """"""
        pass