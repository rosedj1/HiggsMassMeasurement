import ROOT

charge_dict = {
    '+':'-13', 
    '-':'13',
}

color_dict = {
    1:'black',
    2:'red',
    3:'green',
    4:'blue',
    5:'m',
    6:'c',
    7:'y',
    8:'limegreen',
    9:'darkviolet',
    10:'orange',
    11:'yellow',
    12:'pink',
    13:'turquoise',
}

color_dict_RooFit = {
    1:600,#,'kBlue',
    2:632,#'kRed',
    3:416,#'kGreen',
    4:880,#'kViolet',
    5:800,#'kOrange',
    # Now cycle the colors, but go darker:
    6:ROOT.kBlue+2,#'kCyan',
    7:ROOT.kRed+2,#'kMagenta', # looks like kViolet
    8:ROOT.kGreen+2,#'kTeal', # looks like kCyan
    9:ROOT.kViolet+2,#'kPink', # looks like kRed
    10:ROOT.kOrange+2,#'kSpring', # looks like kGreen
    11:920,#'kGray',
    # 6:432,#'kCyan',
    # 7:616,#'kMagenta', # looks like kViolet
    # 8:840,#'kTeal', # looks like kCyan
    # 9:900,#'kPink', # looks like kRed
    # 10:820,#'kSpring', # looks like kGreen
    # From https://root.cern.ch/doc/master/classTColor.html
    # kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416,
    # kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
    # kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
}

process_dct = {
    "DY2mu" : {"process" : r"q#bar{q} #rightarrow Z #rightarrow #mu#mu",
               "m_inv_label" : r"m_{2#mu}",
               "m_inv_lim" : [60, 120],
               "m_inv_nbins" : 120},
    "DY2e"  : {"process" : r"q#bar{q} #rightarrow Z #rightarrow ee",
               "m_inv_label" : r"m_{2e}",
               "m_inv_lim" : [105, 140],
               "m_inv_nbins" : 120},
    "H4mu"  : {"process" : r"H #rightarrow ZZ #rightarrow 4#mu",
               "m_inv_label" : r"m_{4#mu}",
               "m_inv_lim" : [105, 140],
               "m_inv_nbins" : 140},
    "H4e"  : {"process" : r"H #rightarrow ZZ #rightarrow 4e",
               "m_inv_label" : r"m_{4e}",
               "m_inv_lim" : [105, 140],
               "m_inv_nbins" : 140},
    "H2mu"  : {"process" : r"H rightarrow #mu#mu",
               "m_inv_label" : r"m_{2#mu}",
               "m_inv_lim" : [105, 140],
               "m_inv_nbins" : 140},
    "H2e"  : {"process" : r"H rightarrow ee",
               "m_inv_label" : r"m_{2e}",
               "m_inv_lim" : [105, 140],
               "m_inv_nbins" : 140},

}

label_LaTeX_dict = {
    "pT1"  : {"label":r"$p_{T1}^{\mathrm{REC}}$", 
              "independent_label":r"$p_{T}^{\mathrm{REC}}$", 
              "units":"GeV"},
    
    "pT2"  : {"label":r"$p_{T2}^{\mathrm{REC}}$", 
              "independent_label":r"$p_{T}^{\mathrm{REC}}$", 
              "units":"GeV"},
    
    'eta1' : {"label":r"$\eta_{1}^{\mathrm{REC}}$", 
              "independent_label":r"$\eta^{\mathrm{REC}}$", 
              "units":""},
    
    'eta2' : {"label":r"$\eta_{2}^{\mathrm{REC}}$", 
              "independent_label":r"$\eta^{\mathrm{REC}}$", 
              "units":""},
    
    'theta1' : {"label":r"$\theta_{1}^{\mathrm{REC}}$", 
              "independent_label":r"$\theta^{\mathrm{REC}}$", 
                "units":""},
    
    'theta2' : {"label":r"$\theta_{2}^{\mathrm{REC}}$", 
              "independent_label":r"$\theta^{\mathrm{REC}}$", 
                "units":""},
    
    'phi1' : {"label":r"$\phi_{1}^{\mathrm{REC}}$", 
              "independent_label":r"$\phi^{\mathrm{REC}}$",
              "units":""},
    'phi2' : {"label":r"$\phi_{2}^{\mathrm{REC}}$", 
              "independent_label":r"$\phi^{\mathrm{REC}}$",
              "units":""},
    
    'genLep_pt1' : {"label":r"$p_{T1}^{\mathrm{GEN}}$", 
                    "independent_label":r"$p_{T}^{\mathrm{GEN}}$",
                    "units":"GeV"},
    'genLep_pt2' : {"label":r"$p_{T2}^{\mathrm{GEN}}$", 
                    "independent_label":r"$p_{T}^{\mathrm{GEN}}$",
                    "units":"GeV"},
    
    'genLep_eta1' : {"label":r"$\eta_{1}^{\mathrm{GEN}}$", 
                     "independent_label":r"$\eta^{\mathrm{GEN}}$",
                     "units":""},
    'genLep_eta2' : {"label":r"$\eta_{2}^{\mathrm{GEN}}$", 
                     "independent_label":r"$\eta^{\mathrm{GEN}}$",
                     "units":""},
    
    'genLep_phi1' : {"label":r"$\phi_{1}^{\mathrm{GEN}}$", 
                     "independent_label":r"$\phi^{\mathrm{GEN}}$",
                     "units":""},
    'genLep_phi2' : {"label":r"$\phi_{2}^{\mathrm{GEN}}$", 
                     "independent_label":r"$\phi^{\mathrm{GEN}}$",
                     "units":""},
    
    'genLep_theta1' : {"label":r"$\theta_{1}^{\mathrm{GEN}}$", 
                       "independent_label":r"$\theta^{\mathrm{GEN}}$",
                       "units":""},
    'genLep_theta2' : {"label":r"$\theta_{2}^{\mathrm{GEN}}$", 
                       "independent_label":r"$\theta^{\mathrm{GEN}}$",
                       "units":""},
    
    "delta_pT1" : {"label":r"$\Delta p_{T1} \equiv p_{T1}^{\mathrm{REC}} - p_{T1}^{\mathrm{GEN}} $", 
                   "independent_label":r"$\Delta p_{T} \equiv p_{T}^{\mathrm{REC}} - p_{T}^{\mathrm{GEN}} $",
                   "units":"GeV"},
    "delta_pT2" : {"label":r"$\Delta p_{T2} \equiv p_{T2}^{\mathrm{REC}} - p_{T2}^{\mathrm{GEN}} $", 
                   "independent_label":r"$\Delta p_{T} \equiv p_{T}^{\mathrm{REC}} - p_{T}^{\mathrm{GEN}} $",
                   "units":"GeV"},
    
    "delta_pToverGenpT"    : {"label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{GEN}}$", 
                            "independent_label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{GEN}}$", 
                            "units":""},
    "delta_pToverGenpT1"    : {"label":r"$ \Delta p_{T1} \ / p_{T1}^{\mathrm{GEN}}$", 
                            "independent_label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{GEN}}$", 
                            "units":""},
    "delta_pToverGenpT2"    : {"label":r"$ \Delta p_{T2} \ / p_{T2}^{\mathrm{GEN}}$", 
                            "independent_label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{GEN}}$", 
                            "units":""},
    "delta_pToverRecpT" : {"label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{REC}}$", 
                            "independent_label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{REC}}$", 
                            "units":""},
    "delta_pToverRecpT1" : {"label":r"$ \Delta p_{T1} \ / p_{T1}^{\mathrm{REC}}$", 
                            "independent_label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{REC}}$", 
                            "units":""},
    "delta_pToverRecpT2" : {"label":r"$ \Delta p_{T2} \ / p_{T2}^{\mathrm{REC}}$", 
                            "independent_label":r"$ \Delta p_{T} \ / p_{T}^{\mathrm{REC}}$", 
                            "units":""},
    "delta_pToverGenpTsqred" : {"label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{GEN}})^2$", 
                        "independent_label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{GEN}})^2$", 
                        "units":r"$\mathrm{GeV}^{-1}$"},
    "delta_pToverGenpTsqred1" : {"label":r"$ \Delta p_{T1} \ / (p_{T1}^{\mathrm{GEN}})^2$", 
                        "independent_label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{GEN}})^2$", 
                        "units":r"$\mathrm{GeV}^{-1}$"},
    "delta_pToverGenpTsqred2" : {"label":r"$ \Delta p_{T2} \ / (p_{T2}^{\mathrm{GEN}})^2$", 
                        "independent_label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{GEN}})^2$", 
                        "units":r"$\mathrm{GeV}^{-1}$"},
    "delta_pToverRecpTsqred" : {"label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{REC}})^2$", 
                        "independent_label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{REC}})^2$", 
                        "units":r"$\mathrm{GeV}^{-1}$"},
    "delta_pToverRecpTsqred1" : {"label":r"$ \Delta p_{T1} \ / (p_{T1}^{\mathrm{REC}})^2$", 
                        "independent_label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{REC}})^2$", 
                        "units":r"$\mathrm{GeV}^{-1}$"},
    "delta_pToverRecpTsqred2" : {"label":r"$ \Delta p_{T2} \ / (p_{T2}^{\mathrm{REC}})^2$", 
                        "independent_label":r"$ \Delta p_{T} \ / (p_{T}^{\mathrm{REC}})^2$", 
                        "units":r"$\mathrm{GeV}^{-1}$"},
    
    
#     "delta_pToverpT2_squared" : {"label":r"$( p_{T}^{\mathrm{REC}} - p_{T}^{\mathrm{GEN}} ) / (p_{T}^{\mathrm{GEN}})^2 $ GeV$^{-1}$", "units":r"GeV$^{-1}$"},
    
    "delta_R1" : {"label":r"$\Delta R_{1} = \sqrt{  (\Delta \eta_{1})^2 + (\Delta \phi_{1})^2   }$", 
                  "independent_label":r"$\Delta R = \sqrt{  (\Delta \eta)^2 + (\Delta \phi)^2   }$",
                  "units":""},
    "delta_R2" : {"label":r"$\Delta R_{2} = \sqrt{  (\Delta \eta_{2})^2 + (\Delta \phi_{2})^2   }$", 
                  "independent_label":r"$\Delta R = \sqrt{  (\Delta \eta)^2 + (\Delta \phi)^2   }$",
                  "units":""},
    
    "delta_eta1" : {"label":r"$\Delta \eta_{1} = \eta_{1}^{\mathrm{REC}} - \eta_{1}^{\mathrm{GEN}}$", 
                    "independent_label":r"$\Delta \eta = \eta^{\mathrm{REC}} - \eta^{\mathrm{GEN}}$",
                    "units":""},
    "delta_eta2" : {"label":r"$\Delta \eta_{2} = \eta_{2}^{\mathrm{REC}} - \eta_{2}^{\mathrm{GEN}}$", 
                    "independent_label":r"$\Delta \eta = \eta^{\mathrm{REC}} - \eta^{\mathrm{GEN}}$",
                    "units":""},
    
    "delta_theta1" : {"label":r"$\Delta \theta_{1} = \theta_{1}^{\mathrm{REC}} - \theta_{1}^{\mathrm{GEN}}$", 
                      "independent_label":r"$\Delta \theta = \theta^{\mathrm{REC}} - \theta^{\mathrm{GEN}}$",
                      "units":""},
    "delta_theta2" : {"label":r"$\Delta \theta_{2} = \theta_{2}^{\mathrm{REC}} - \theta_{2}^{\mathrm{GEN}}$", 
                      "independent_label":r"$\Delta \theta = \theta^{\mathrm{REC}} - \theta^{\mathrm{GEN}}$",
                      "units":""},
    
    "delta_phi1" : {"label":r"$\Delta \phi_{1} = \phi_{1}^{\mathrm{REC}} - \phi_{1}^{\mathrm{GEN}}$", 
                    "independent_label":r"$\Delta \phi = \phi^{\mathrm{REC}} - \phi^{\mathrm{GEN}}$",
                    "units":""},
    "delta_phi2" : {"label":r"$\Delta \phi_{2} = \phi_{2}^{\mathrm{REC}} - \phi_{2}^{\mathrm{GEN}}$", 
                    "independent_label":r"$\Delta \phi = \phi^{\mathrm{REC}} - \phi^{\mathrm{GEN}}$",
                    "units":""},
    
    "d0BS1" : {"label":r"$d_{0, \mathrm{lep1}}^{ \mathrm{BS} }$", 
              "independent_label":r"$d_{0}^{ \mathrm{BS} }$",
              "independent_label_ROOT":r"d_{0}^{\mathrm{BS}}", 
               "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    "d0BS2" : {"label":r"$d_{0, \mathrm{lep2}}^{ \mathrm{BS} }$", 
              "independent_label":r"$d_{0}^{ \mathrm{BS} }$",
              "independent_label_ROOT":r"d_{0}^{\mathrm{BS}}",
               "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    "d0PV1" : {"label":r"$d_{0, \mathrm{lep1}}^{ \mathrm{PV} }$", 
              "independent_label":r"$d_{0}^{ \mathrm{PV} }$", 
               "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    "d0PV2" : {"label":r"$d_{0, \mathrm{lep2}}^{ \mathrm{PV} }$", 
              "independent_label":r"$d_{0}^{ \mathrm{PV} }$", 
               "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    
    "qd0BS1" : {"label":r"$q(\mu_{1}^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{BS} }$", 
              "independent_label":r"$q(\mu^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{BS} }$",
              "independent_label_ROOT":r"qd_{0}^{\mathrm{BS}}",
                 "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    "qd0BS2" : {"label":r"$q(\mu_{2}^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{BS} }$", 
              "independent_label":r"$$q(\mu^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{BS} }$",
              "independent_label_ROOT":r"qd_{0}^{\mathrm{BS}}",
                 "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    "qd0PV1" : {"label":r"$q(\mu_{1}^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{PV} }$", 
              "independent_label":r"$q(\mu^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{PV} }$", 
                 "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    "qd0PV2" : {"label":r"$q(\mu_{2}^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{PV} }$", 
              "independent_label":r"$q(\mu^{\pm, \mathrm{REC} }) * d_{0}^{ \mathrm{PV} }$", 
                 "units":"cm", "default_bin_limits":[-0.01, 0.01, 0.0002], "default_x_limits":[-0.012, 0.012]},
    
    "massZ"    : {"label":r"$m_{\mu^{+}\mu^{-}}$",        "independent_label":r"$m_{\mu^{+}\mu^{-}}$", "units":"GeV", "default_bin_limits":[60,120,0.4], "default_x_limits":[50, 130]},
    "massZErr" : {"label":r"$\delta m_{\mu^{+}\mu^{-}}$", "independent_label":r"$\delta m_{\mu^{+}\mu^{-}}$", "units":"GeV", "default_bin_limits":[0,2,0.01],   "default_x_limits":[-0.1, 2.1]},

    
# Unused branches: 
# ', 'massZErr', 'massZ_vtx', 'massZ_vtx_FSR', 'massErrZ_vtx',
#        'massErrZ_vtx_FSR', 'massZ_vtxChi2', 'massZ_vtx_BS',
#        'm1', 'm2', 'Id1', 'Id2', 'Tight1', 'Tight2', 'pterr1', 'pterr2', 'weight','GENmass2l', 
#         'nFSRPhotons',  'genLep_p1', 'genLep_p2', 'p1', 'p2', 
#         'delta_Rtheta1', 'delta_Rtheta2', 'delta_pToverGenpT2'
}