import subprocess

def shell_cmd(cmd_str):
    """Execute `cmd_str` from your shell and return the CompletedProcess.
    
    Say you want to execute: `ls -l some_dir/`
    Then do: shell_cmd("ls -l some_dir/")
    """
    cmd_parts_ls = cmd_str.split(" ")
    return subprocess.run(cmd_parts_ls)