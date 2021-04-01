import ROOT as r

ID_charge_dct = {
     13 : -1,
    -13 :  1,
}

class MyParticle:
    """A class to organize generic particle info."""

    def __init__(self):
        self.isFermion = None
        self.isLepton = None
        self.flavor = None
       
        self.pT = None
        self.eta = None
        self.phi = None
        self.mass = None
       
        self.gen_pT = None
        self.gen_eta = None
        self.gen_phi = None
        self.gen_mass = None

        self.pT_withFSR = None
        self.eta_withFSR = None
        self.phi_withFSR = None
        self.mass_withFSR = None

        self.charge = None
        self.ID = None
        self.d0 = None

    def charge_from_ID(self, ID):
        return ID_charge_dct[ID]
    
    def set_PtEtaPhiMass(self, pT, eta, phi, mass):
        self.pT = pT
        self.eta = eta
        self.phi = phi
        self.mass = mass

    def set_GENPtEtaPhiMass(self, gpT, geta, gphi, gmass):
        self.gen_pT = gpT
        self.gen_eta = geta
        self.gen_phi = gphi
        self.gen_mass = gmass

    def set_PtEtaPhiMass_withFSR(self, pTFSR, etaFSR, phiFSR, massFSR):
        self.pT_withFSR = pTFSR
        self.eta_withFSR = etaFSR
        self.phi_withFSR = phiFSR
        self.mass_withFSR = massFSR

    def get_PtEtaPhiMass(self):
        """Return a 4-tup of muon's RECO pT, eta, phi, mass."""
        return (self.pT, self.eta, self.phi, self.mass)

    def get_GENPtEtaPhiMass(self):
        """Return a 4-tup of muon's GEN pT, eta, phi, mass."""
        return (self.gen_pT, self.gen_eta, self.gen_phi, self.gen_mass)
 
    def get_PtEtaPhiMass_withFSR(self):
        """Return a 4-tup of muon's pT, eta, phi, mass, corrected for FSR."""
        return (self.pT_withFSR, self.eta_withFSR, self.phi_withFSR, self.mass_withFSR)

    def verify_charge(self):
        """Return True if charge of particle is compatible with its ID."""
        q_ = ID_charge_dct[self.ID]
        assert abs(q_) == abs(self.charge)

    def get_LorentzVector(self, kind=""):
        """Return a ROOT.Math.LorentzVector object of this particle.
        
        kind : str
            Can be: "gen", "reco", or "withFSR".
        """
        # k = kind.lower()
        # print(f"k is {k}")
        if kind in "gen":
            tup = self.get_GENPtEtaPhiMass()
        elif kind in "reco":
            tup = self.get_PtEtaPhiMass()
        elif kind in "withFSR":
            tup = self.get_PtEtaPhiMass_withFSR()
        else:
            raise ValueError(f"Parameter `kind` ({kind}) not understood.")
        return r.Math.PtEtaPhiMVector(*tup)

#----------------------------------------#

class MyMuon(MyParticle):
    """A class to organize muon info."""

    def __init__(self, q):
        """Initialize a muon using it's charge (q)."""
        self.isFermion = True
        self.isLepton = True
        self.flavor = "muon"
        self.gen_mass = 0.1057  # GeV
        self.charge = q

        self.ID = self.ID_from_charge(q)

    def ID_from_charge(self, q):
        if q == -1:
            return 13
        elif q == 1: 
            return -13
        else:
            raise TypeError(f"Based on the charge ({q}), this is not a muon!")