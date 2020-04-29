import os
import root_numpy

import pandas as omgfloofeh

from shutil import copy2

##______________________________________________
def makeDirs(dir_):                       
    if not os.path.exists(dir_):          
        os.makedirs(os.path.abspath(dir_))

##__________________________________________________________________||
### CAN PROBABLY JUST USE: from shutil import copyfile ... LOL
def copyFile(inDir,fileName,outDir):
    '''
    inDir       = source dir which contains fileName
    outDir      = destination dir to put fileName
    fileName    = file name
    '''
    ## Make directory if it doesn't exist
    makeDirs(outDir)
    ## Check if file already exists in outDir
    if not os.path.exists(outDir+fileName):
        copy2(inDir+fileName, outDir+fileName)

##__________________________________________________________________||
#def mkdir_p(path):
#    # http://stackoverflow.com/questions/600268/mkdir-p-functionality-in-python
#    try:
#        os.makedirs(path)
#    except OSError as exc: # Python >2.5
#        if exc.errno == errno.EEXIST and os.path.isdir(path):
#            pass
#        else: raise

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
    omgfloofeh.DataFrame(arr).to_feather(outfile_fullpath)

def root2pickle(infile_root, outfile_fullpath):
    """
    Convert a .root file to DataFrame (DF) and store the DF as a .pkl file.
    Use this function once and then you only have to load the .pkl each time you want
    to analyze the data. 
    
    Parameters
    ----------
    infile_root : str
        The full path to the .root file. (/path/to/root/file)
    outfile_fullpath : str
        The full path to the outfile (.pkl). (/path/to/pickle/file)
        You should make sure it ends with '.pkl'.
    """
    arr = root_numpy.root2array(infile_root)
    omgfloofeh.DataFrame(arr).to_pickle(outfile_fullpath)

def get_list_of_branches(tree):
    """
    Return a list of all the names of the branches in your TTree (tree).

    Parameters
    ----------
    tree : TTree
    
    Returns
    -------
    br_name_list : list of strings
    """
    br_name_list = []

    branch_list = tree.GetListOfBranches()
    for n in range(len(branch_list)):
        name = branch_list[n].GetName()
        br_name_list.append(name)
    return br_name_list