import vaex  # The holy grail.
from ROOT import TFile
from Utils_vaex.vaex_fns import prepare_vaex_df

infile_path_data_dict = {
    # Key must be: "{data_type}_{year}_{name}_{extension}"
    "MC_2016_DY_hdf5"      : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2016/MC_2016_DY.hdf5",
    "MC_2017_DY_hdf5"      : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_DY.hdf5",

    "MC_2017_Jpsi_hdf5"    : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Jpsi.hdf5",
    "MC_2017_Jpsi_root"    : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Jpsi.root",
    "MC_2017_Upsilon_root" : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2017/MC_2017_Upsilon.root",

    "MC_2018_DY_hdf5"      : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2018/MC_2018_DY.hdf5",
    "MC_2018_Jpsi_hdf5"    : "/ufrc/avery/rosedj1/HiggsMassMeasurement/Samples/MC_2018/MC_2018_Jpsi.hdf5",
}

sample_info_dict = {
    "Jpsi"    : {"LaTeX_name" : r"J/\psi", 
                 "inv_mass_cut_lim" : [2.9, 3.3],
                 "dR_cut" : 0.005},
    "DY"      : {"LaTeX_name" : r"\mathrm{Z}", 
                 "inv_mass_cut_lim" : [60.0, 120.0],
                 "dR_cut" : 0.002},
    "Upsilon" : {"LaTeX_name" : r"\Upsilon", 
                 "inv_mass_cut_lim" : [8.5, 11.0]},
    "Higgs"   : {"LaTeX_name" : r"H", 
                 "inv_mass_cut_lim" : [105.0, 140.0]},
}

class Sample:
    def __init__(self, data_type, year, name, extension):
        """
        A data sample. 

        Parameters
        ----------
        data_type : str
            E.g., "MC" or "Data"
        year : str
            E.g., "2016", "2017", "2018"
        name : str
            E.g., "Jpsi", "Upsilon", "DY", "Higgs"
        extension : str
            File extension of data, like: "hdf5", "root"
        """
        self.data_type = data_type
        self.year = year
        self.name = name
        self.extension = extension
        
        self.get_sample_info()
        self.get_data_path()
        self.get_data()
        self.prep_data()
    
    def get_sample_info(self):
        self.LaTeX_name       = sample_info_dict[self.name]["LaTeX_name"]
        self.inv_mass_cut_lim = sample_info_dict[self.name]["inv_mass_cut_lim"]
        self.dR_cut           = sample_info_dict[self.name]["dR_cut"]

    def get_data_path(self):
        """
        Based on data_type, year, name, and extension, 
        go get file path where data are stored.
        """
        key = f"{self.data_type}_{self.year}_{self.name}_{self.extension}"
        self.infile_path = infile_path_data_dict[key]
        
    def get_data(self):
        if self.extension in "hdf5":
            # Assume a vaex DataFrame.
            self.vdf_original = vaex.open(self.infile_path)
        elif self.extension in "root":
            # FIXME: I think this breaks ROOT.
            self.tfile = TFile.Open(self.infile_path, "read") 
            
    def prep_data(self):
        """
        Adds more columns to vaex DataFrame, like kinematic differences.
        Also treats muons independently by making a 'concatenated' VDF.
        """
        if self.extension in "hdf5":
            self.vdf_prepped = prepare_vaex_df(self.vdf_original)