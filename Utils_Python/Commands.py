import subprocess

def shell_cmd(cmd_str, get_stdout=True, outfile=None, verbose=False):
    """Execute `cmd_str` from your shell and return the CompletedProcess.
    
    Say you want to execute: `ls -l some_dir/`
    Then do: shell_cmd("ls -l some_dir/")

    Parameters
    ----------
    cmd_str : str
        String of commands to execute in shell.
    get_stdout : bool
        If True, will store stdout.
        To print stdout as a str:
        `results = subprocess.run(<cmds>)`
        `results.stdout`

        NOTE: If `outfile` is provided, then stdout will not be printed.
              Instead it will be written to `outfile`.
    outfile : str or None
        File path or name of file. Write stdout to this file.
    """
    if verbose:
        print(f"Executing command: {cmd_str}")
    cmd_parts_ls = cmd_str.split(" ")
    if outfile is not None:
        with open(outfile, "w") as f:
            return subprocess.run(cmd_parts_ls, stdout=f)
    if get_stdout:
        return subprocess.run(cmd_parts_ls, stdout=subprocess.PIPE, universal_newlines=True)
    return subprocess.run(cmd_parts_ls)