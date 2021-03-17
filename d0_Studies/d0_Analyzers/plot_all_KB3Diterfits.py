"""Make a PDF of all iterated Gaussian fits from KB3Ds."""

from Utils_Python.Utils_Files import open_pkl
from Utils_ROOT.Printer import CanvasPrinter

infile = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/pickles/final_muon_coll/MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds.pkl"
outpath_pdf = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/DeriveCorr/MC2016DY/plots/test/MC2016DY_finalmuoncoll_allitergaussfitsonKB2DsandKB3Ds.pdf"

mu_coll = open_pkl(infile)
printer = CanvasPrinter(show_plots=1)
printer.make_plots_pretty()
# plots = [kb3d.frame_dpTOverpT for kb2d in mu_coll.KinBin2D_dict.values() for kb3d in kb2d.KinBin3D_dict.values()]
# printer.make_pdf_of_plots(plots, outpath_pdf)

#--- ReferenceError: attempt to access a null-pointer
#--- gr_ls = list(mg.GetListOfGraphs())
# printer.canv.Print(outpdf_path + "[")
# mu_coll.make_all_multigraphs(eta_range, scale_by_1divpT=False)
# print("Made all MultiGraphs with scale_by_1divpT=False.")
# mu_coll.make_all_multigraphs(eta_range, scale_by_1divpT=True)
# print("Made all MultiGraphs with scale_by_1divpT=True.")
# print(f"multigraph_ls = {mu_coll.multigraph_ls}")
# print(f"multigraph_1divpT_ls = {mu_coll.multigraph_1divpT_ls}")
# mu_coll.draw_all_multigraphs(outpdf_path, printer)  # Needs: self.multigraph_ls, look into: draw_mg_and_fits
# print("Drew all MultiGraphs.")
# all_plots = mu_coll.get_all_kb3d_plots("dpT/pT iterfit", "qd0 hists") + mu_coll.get_all_kb2d_plots("dpT/pT hists", "qd0 hists")
# printer.draw_plots(all_plots, outpdf_path)
# printer.canv.Print(outpdf_path + "]")