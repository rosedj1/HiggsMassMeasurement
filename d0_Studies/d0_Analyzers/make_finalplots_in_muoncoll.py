# from ROOT import TCanvas
import os
from Utils_Python.Utils_Files import open_pkl, make_dirs, check_overwrite
from Utils_ROOT.Printer import CanvasPrinter

outpdf_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/kb2ditergaussfits_beforeafterGeoFitcorr/test/test02.pdf"
# outtxt_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/output/"
# mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p0_d0_0p01.pkl")
mu_coll = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/muoncoll_itergaussfitsonKB2Ds_0p0_d0_0p01_1000bins_neg0p2xlim0p2.pkl")

overwrite = 0

# kb2d = open_pkl("/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_nogenmatching_0p0_d0_0p01/MC2016DY_individKB2D_withitergaussfitsonKB3Ds_nogenmatching_0p0_d0_0p01_0p0eta2p4_10p0pT15p0.pkl")
make_dirs(os.path.dirname(outpdf_path))
check_overwrite(outpdf_path, overwrite=overwrite)
printer = CanvasPrinter(show_plots=0)
printer.make_plots_pretty()
# printer.make_pdf_of_plots()
# kb3d.frame_dpTOverpT.Draw()

# Let's make some plots.

# mu_coll.make_plots_beforeafterpTcorr()
# plot_ls = mu_coll.hist_inclusive_ls
# c = TCanvas()
c = printer.canv
c.Print(outpdf_path + "[")
# with open(outtxt_path, "w") as f:
for ct, kb2d in enumerate(mu_coll.KinBin2D_dict.values()):
    if ct == 1:
        print(vars(kb2d))
    # kb2d.overwrite_muon_info(delete_all=True)
    kb2d.make_beforeafterpTcorr_frames()
    kb2d.draw_beforeafterpTcorr()
    c.Print(outpdf_path)
    # plot_ls.append(kb2d.frame_dpTOverpT)
    # plot_ls.append(kb2d.frame_dpTOverpT_corr)
    print(f"{kb2d.eta_range}, {kb2d.pT_range}, {kb2d.}, {kb2d.}, {kb2d.}, {kb2d.}, {kb2d.}")
c.Print(outpdf_path + "]")
# printer.make_pdf_of_plots(plot_ls, outpdf_path)
# save_to_pkl(mu_coll, outpkl_path)