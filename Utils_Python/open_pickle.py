import pickle
import sys

inpkl = sys.argv[1]

with open(inpkl, "rb") as p:
    d = pickle.load(p)
