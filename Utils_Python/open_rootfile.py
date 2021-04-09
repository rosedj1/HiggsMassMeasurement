"""ROOT File Opener

Purpose:
    From the command line, open a root file in iPython, set the TFile object to
    `f`. If a TTree name was provided then open the TTree as `t`.

Syntax:
    In your terminal, do:
    openrootfile /home/work/myfile.root mytree

    This supposes you have a root file at: /home/work/myfile.root
    Inside the root file, there is a TTree called "mytree".

    You will be taken into the iPython interpreter in which the root file has
    been attached as `f` and mytree is set as `t`.

Note:
    Put the following in your .bash_profile.
    Be sure to replace the angled brackets with the actual path:

    openrootfile() {
        ipython -i <fullpath/to/thisscript.py> ${1} ${2}
    }
"""
import sys
print("Imported sys.")
import ROOT as rt
print("Imported ROOT as rt.")

# Open file.
infile = sys.argv[1]
f = rt.TFile(infile)
print("Root file attached as TFile object `f`.")
f.ls()

# Open TTree if provided.
try:
    intree = str(sys.argv[2])
    t = f.Get(intree)
    print(f"TTree ('{intree}') opened as `t`.")
except IndexError:
    pass
