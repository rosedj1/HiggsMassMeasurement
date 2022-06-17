import os
import json
from shutil import copy2
import pickle
import subprocess

def make_dirs(*d, verbose=True):                       
    """If the directory (d) does not exist, then make it. 
    
    NOTE: Makes directories recursively.
    """
    for thisdir in d:
        if not os.path.exists(thisdir):
            if (verbose):
                print(f"[INFO] Directory not found. Creating it:\n{thisdir}")
            os.makedirs(os.path.abspath(thisdir))

def check_overwrite(*outfiles, overwrite=False):
    """Raise an error if outfile exists and overwrite == False."""
    for f in outfiles:
        if os.path.exists(f) and not overwrite:
            emsg = (
                f"Not allowed to overwrite file since it already exists"
                f"\n{f}\n"
                f"To write over this file, set overwrite = True."
                )
            raise RuntimeError(emsg)
    
def make_str_title_friendly(str_, keep_whitespace=False):
    title = str_.replace('-13','muPlus')
    title = title.replace('+13','muMinus')
    title = title.replace('&&','_')
    title = title.replace('&','_')
    title = title.replace('<','_lt_')
    title = title.replace('>','_gt_')
    if not (keep_whitespace): 
        title = title.replace(' ','')
    title = title.replace('.','p')
    title = title.replace('(','')
    title = title.replace(')','')
    title = title.replace('-','neg')
    title = title.replace('+','pos')
    title = title.replace('||','or')
    title = title.replace('|','or')
    return title

def open_pkl(inpath):
    """Open a pickled object stored in `inpath`."""
    with open(inpath, "rb") as p: 
        return pickle.load(p)

def open_json(inpath):
    """Open a pickled object stored in `inpath`."""
    with open(inpath, "rb") as p: 
        return json.load(p)

def save_to_json(obj, outpath, overwrite=True, sort_keys=False):
    """Write one obj to 'outpath.json'."""
    check_overwrite(outpath, overwrite=overwrite)
    with open(outpath, 'w') as f:
        json.dump(obj, f, indent=4, sort_keys=sort_keys)
    print(f"[INFO] JSON file written:\n{outpath}\n")

def save_to_pkl(obj, outpkl_path, overwrite=True):
    """Write one obj to pickle."""
    check_overwrite(outpkl_path, overwrite=overwrite)
    with open(outpkl_path, 'wb') as output:
        pickle.dump(obj, output, protocol=2)
    print(f"[INFO] Pickle file written:\n{outpkl_path}\n")

def replace_value(old, new, script):
    """Use the `sed` command to replace `old` with `new` in `script`."""
    cmd = ["sed", "-i", f"s|{old}|{new}|g", script]
    output = subprocess.run(cmd)

def root2feather(infile_root, outfile_fullpath):
    """
    Convert a .root file to DataFrame (DF) and store the DF as a binary .feather file.
    Use this function once and then you only have to load the .feather each time you want
    to analyze the data. 
    
    Parameters
    ----------
    infile_root : str
        The full path to the .root file. (/path/to/root/file)
    outfile_fullpath : str
        The full path to the outfile (.feather). (/path/to/feather/file)
        You should make sure it ends with '.feather'.
    """
    arr = root_numpy.root2array(infile_root)
    pd.DataFrame(arr).to_feather(outfile_fullpath)

class Py3toPy2Converter:
    """
    # FIXME: Not complete.
    Class to insert the proper encoding line in Python2 files
    but take them away in Python3 files.
    """
    def __init__(self, file_ls):
        import subprocess
        self.file_ls = file_ls

    def insert_fstring_coding(self, file_ls):
        """Insert f-string encoding at the top of each file."""
        fstring_code_str = "# -*- coding: future_fstrings -*-"
        # NOTE: The code below is under development.
        # Not sure if the `sed` cmd will interfere with regex.
        # sed '1 s/^/# -*- coding: future_fstrings -*-\n/' file1

# I don't like this function.
# Doesn't provide much flexibility.
# def copy_file(inDir, fileName, outDir):
#     '''
#     Make directory (outDir) if it does not exist, then copy file (inDir + filename) to outDir.
#     Also checks to see if file already exists in outDir. If so, do nothing.

#     inDir       = source dir which contains fileName
#     outDir      = destination dir to put fileName
#     fileName    = file name
#     '''
#     ## Make directory if it doesn't exist.
#     make_dirs(outDir)
#     ## Check if file already exists in outDir.
#     if not os.path.exists(outDir + fileName):
#         copy2(inDir + fileName, outDir + fileName)