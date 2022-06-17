#!/usr/bin/env ipython -i --
"""ROOT File Opener

Purpose:
    From the command line, open a root file in iPython, set the TFile object to
    `f`. If a TTree path was provided then open the TTree as `t`.

Syntax:
    In your terminal, you can do:
    * openrootfile /home/work/myfile.root
    * openrootfile /home/work/myfile.root treepath
    * openrootfile /home/work/myfile.root treepath -v (verbose)

    This supposes you have a root file at: /home/work/myfile.root
    Inside the root file, there is a TTree called "treepath".

    You will be taken into the iPython interpreter in which the root file has
    been attached as `f` and treepath is set as `t`.

Notes:
    Put this function in your ~/.bash_profile:

    openrootfile() {
        if [ $# -eq 1 ]; then
            # Infile.
            ipython -i /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/Utils_Python/openrootfile.py -- -f ${1}
        elif [ $# -eq 2 ]; then
            # Infile, treepath.
            ipython -i /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/Utils_Python/openrootfile.py -- -f ${1} -t ${2}
        elif [ $# -eq 3 ]; then
            # Infile, treepath, and verbose.
            ipython -i /cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/Utils_Python/openrootfile.py -- -f ${1} -t ${2} ${3}
        fi
    }
"""
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument(
    '-f', '--file', dest = 'infile', type = str,
    help = "Input root file.")
parser.add_argument(
    '-t', '--treepath', dest = 'treepath', type = str, default="",
    help = "Path to TTree.")
parser.add_argument(
    '-v', '--verbose', dest="verbose", action="store_true",
    help="If TTree is found, print info. Default is False.")
        
args = parser.parse_args()
infile = args.infile
treepath = args.treepath
verbose = args.verbose

import ROOT as rt
print("Imported ROOT as rt.")

# Open file.
f = rt.TFile.Open(infile)
print(
    f"Root file attached as TFile object `f`.\n"
    f"#=== f.ls() ===#"
    )
f.ls()

# Open TTree if provided.
try:
    t = f.Get(treepath)
    print(
        f"TTree ('{treepath}') opened as `t`.\n"
        f"TTree.GetEntries() = {t.GetEntries()}"
        )
    if verbose:
        print(f"#=== TTree.Show(0) ===#")
        t.Show(0)
except (IndexError, AttributeError):
    pass
