import os
from shutil import copy2

##__________________________________________________________________||
### CAN PROBABLY JUST USE: from shutil import copyfile ... LOL
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
def makeDirs(dir_):                       
    """
    If the directory (dir_) does not exist, then make it. 
    """
    if not os.path.exists(dir_):          
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
