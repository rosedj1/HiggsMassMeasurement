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
    #print(f"[WARNING] Failed to open TTree.")
