### Reading root files with vaex
I have not been able to convert root files directly into 
hdf5 files that are readable by vaex. 
The closest I came was using the shell command `root2hdf5`
(after doing cmsenv). Vaex was able to read it in, but the
resulting vaex DataFrame essentially had a single column 
called "passedEvents". I wasn't able to parse any of it.

Instead, I have found two options to convert a root file
into a vaex DF (VDF). Let's call them (A) and (B):

*Option A:*
root -> array -> dataframe -> csv -> hdf5
     (1)      (2)          (3)    (4)

(1) from root_numpy import root2array
    arr = root2array("path/to/file.root", "passedEvents")

(2) import numpy
    numpy.savez_compressed("compressed_arr.npz", arr)

(3) Move npz to local. 
    import pandas
    df = pandas.DataFrame(numpy.load("compressed_arr.npz")['arr_0'])

(4) import vaex
    vdf = vaex.from_pandas(df)
    vdf.export_hdf5("sexy_file.hdf5")

Pros: 
  - The npz is relatively light-weight. 
Cons: 
  - You will have to store a huge DF in local 
    memory during steps (3) and (4).
----------

*Option B:*
root -> array -> dataframe -> hdf5
     (1)      (2)          (3)

(1) from root_numpy import root2array
    arr = root2array("path/to/file.root", "passedEvents")

(2) import pandas 
    df = pandas.DataFrame(arr)

(3) df.to_csv("bigass_file.csv")

(4) import vaex
    vdf = vaex.from_csv("bigass_file.csv")
    vdf.export_hdf5("sexy_file.hdf5")

Pros: 
  - Hypothetically could be done 100% remotely, 
    on melrose say, if we could only get vaex to work!
  - More straightforward method. 
Cons: 
  - The csv file is usually enormous (~50 GB).
----------
To read the vdf, do: 
    vdf = vaex.open("sexy_file.hdf5")
This file will be ~instantly read.
It will not be stored in memory!
Values are retrieved only when needed.
<3 vaex
