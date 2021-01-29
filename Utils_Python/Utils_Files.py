import os
from shutil import copy2

##__________________________________________________________________||
def copyFile(inDir,fileName,outDir):
    '''
    Make directory (outDir) if it does not exist, then copy file (inDir + filename) to outDir.
    Also checks to see if file already exists in outDir. If so, do nothing.

    inDir       = source dir which contains fileName
    outDir      = destination dir to put fileName
    fileName    = file name
    '''
    ## Make directory if it doesn't exist
    makeDirs(outDir)
    ## Check if file already exists in outDir
    if not os.path.exists(outDir+fileName):
        copy2(inDir+fileName, outDir+fileName)

##______________________________________________
def makeDirs(dir_, verbose=True):                       
    """
    If the directory (dir_) does not exist, then make it. 
    Makes directories recursively.
    """
    if not os.path.exists(dir_):
        if (verbose):
            print(f"[INFO] Directory not found. Creating {dir_}.")    
        os.makedirs(os.path.abspath(dir_))

##__________________________________________________________________||
#def mkdir_p(path):
#    # http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
#    try:
#        os.makedirs(path)
#    except OSError as exc: # Python >2.5
#        if exc.errno == errno.EEXIST and os.path.isdir(path):
#            pass
#        else: raise

def check_overwrite(outfile, overwrite=False):
    """
    Raises an error if outfile exists and overwrite == False. 
    """
    if os.path.exists(outfile) and not (overwrite):
        err_msg = (
            f"Not allowed to overwrite file since it already exists"
            f"\n{outfile}\n"
            f"To write over this file, set overwrite = True.\n"
            )
        raise RuntimeError(err_msg)
    
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
    Class to insert the proper encoding line in Python2 files
    but take them away in Python3 files.
    # FIXME: Not complete.
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