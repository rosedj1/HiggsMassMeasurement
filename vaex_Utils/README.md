## Reading root files with vaex
I have not been able to convert root files directly into 
hdf5 files that are readable by vaex. 
The closest I came was using the shell command `root2hdf5`
(after doing `cmsenv`). Vaex was able to read the **.hdf5** file in, but the
resulting vaex DataFrame essentially had a single column 
called "passedEvents". I wasn't able to parse any of it.

Instead, I have discovered two options to convert a **.root** file
into a vaex-readable **.hdf5** file. Let's call the Options (A) and (B):

### Option A
    root -> array -> pandas DataFrame -> vaex DataFrame -> hdf5
        (1)      (2)                 (3)               (4)
        
**Pros:** The npz is relatively light-weight. 
  
**Cons:** You will have to store a huge DF in local 
memory during steps (3) and (4).

(1) Convert root file to array (yikes).
```python
from root_numpy import root2array
arr = root2array("path/to/file.root", "passedEvents")
```
(2) Convert array to pandas DataFrame (DF).
```python
import pandas
df = pandas.DataFrame(arr)
```
(3) Convert pandas DF to vaex DF.
```python
import vaex
vdf = vaex.from_pandas(df)
```
(4) Export the vaex DF to a vaex-readable **.hdf5** file.
```python
vdf.export_hdf5("sexy_file.hdf5")
```

##### Recommendation: 
Are people yelling at you because you are using too much RAM 
on the remote machines? Then after step (1), save the array 
as an `.npz` file:

(1b)
```python
import numpy
numpy.savez_compressed("compressed_arr.npz", arr)
```

Then move it to another machine that won't get you in trouble.
Then replace step (2) above with (2'):

(2') Uncompress the **.npz** file and load it into a pandas DF:
```python
import pandas
df = pandas.DataFrame(numpy.load("compressed_arr.npz")['arr_0'])
```

---

### Option B
    root -> array -> dataframe -> csv -> hdf5
        (1)      (2)          (3)    (4)

**Pros:**
- Hypothetically could be done 100% remotely, 
  if you can import `vaex` on a remote machine
- More straightforward method. 
  
**Cons:**
- The **.csv** file it makes can be enormous (~50 GB).

Follow steps (1), (2) above, then:

(3) Convert the pandas DF to a **.csv** file. 
```python
df.to_csv("bigass_file.csv")
```

(4) Load the **.csv** file into a vaex DF:
```python
import vaex
vdf = vaex.from_csv("bigass_file.csv")
```

Finally, step (5) above.

----------

To read the vdf, do: 
    vdf = vaex.open("sexy_file.hdf5")
This file will be ~instantly read.
It will not be stored in memory!
Values are retrieved only when needed.
<3 vaex
