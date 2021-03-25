/**
This code produces a PDF which shows an m4mu distribution,
before and after pT corrections have been applied per muon.
This is part of the GeoFit study (also called "ad hoc pT correction").

Both the uncorr and corr m4mu dists are then fit with DSCB
functions. These are unbinned fits.
The fit parameters are put on the plot.

A root file is used as an input.
Inside the root file, you should have a TTree with 2 branches:
- m4mu and m4mu_corr

Author: Jake Rosenzweig
OG Date: 2020-07ish
Updated: 2021-01-20
**/
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <stdio.h>
#include <map>
#include <utility>
#include <iterator>

#include "TROOT.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TLegend.h"
#include "TMultiGraph.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TMath.h"
#include "TTree.h"
#include "TTreeIndex.h"
#include "TH2F.h"
#include "TLatex.h"
#include "TLine.h"
#include "TGraphAsymmErrors.h"
#include "Math/QuantFuncMathCore.h"

#include "TSystem.h"
#include "TStyle.h"
#include "TPaveText.h"

#include "TPaveLabel.h"
#include "TLegend.h"

#include "TLorentzRotation.h"
#include "TVector3.h"
#include "TLorentzVector.h"
//
//#include <vector>
//#include <fstream>
////
//#include "TRandom3.h"
//
//#include "RooRealVar.h"
//#include "RooArgSet.h"
//#include "RooGaussian.h"
//#include "RooBreitWigner.h"
//#include "RooProdPdf.h"
//#include "RooDataSet.h"
//#include "RooGlobalFunc.h"
//#include "RooDataHist.h"
//#include "RooHistPdf.h"
//#include "RooCBShape.h"
//#include "RooMinuit.h"
//#include "RooFormulaVar.h"
//#include "RooAddPdf.h"
//#include "RooGenericPdf.h"
//
//#include "RooPlot.h"
//
//using namespace std;
//
//void MioSkim(){
//
////   std::vector<float>* Z_mass
#include "RooMyPDF_DSCB.h"
#include "RooRealVar.h"

//using namespace RooFit;

void DSCB_fit_m4mu_beforeafterpTcorr(Bool_t draw = false) {

//----- User Params -----//
Double_t m4mu_min = 105.0;  // 121.2, 128.4
Double_t m4mu_max = 140.0;
Int_t n_bins = 100;
string year = "2018";
string fs = "ggH";
string derive_from_sample = "q#bar{q} #rightarrow Z #rightarrow 2#mu"; // Sample from which pT corr factors were derived.
// string derive_from_sample = "DY+J/#psi"; // Sample from which pT corr factors were derived.
bool plot_residuals = true;  // If false, then ratio of hists will be plotted.
TString infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromAdHocfactors_fullstats_zerointerc_new100muonsperregion.root";
// TString infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2016_m4mu_m4mucorrfromAdHocfactors_fullstats_zerointerc_new.root";
// TString infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018ggH_m4muvals_fullstats_usingXunwuHmumucorrfactorsfrommacro_withoutFSR.root";
// TString infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018_m4mu_m4mucorrfromDYJpsifactors_fullstats_noFSR.root";
TString outfile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/plots/applypTcorrplots/CorrFromMC/MC2018ggH_m4mu_DSCBfit_beforeaftercorr_AdHocfactors_zerointerc_new100muonsperregion.pdf";
/* TString infile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/rootfiles/ggH_skimmed/MC2018ggH_m4muvals_fullstats_usingXunwuHmumucorrfactorsfrommacro_withoutFSR.root";
TString outfile_path = "/cmsuf/data/store/user/t2/users/rosedj1/HiggsMassMeasurement/d0_studies/Plots/applypTcorrplots/CorrFromMC/MC2018_m4mu_ggH_DSCBfit_corrfromMCDY2018Xunwusmacro_withoutFSR.pdf"; */
// This is the color of the fits. 
// The data points will be darker (+2).
Int_t color_line = kBlue;
Int_t color_line_corr = kRed;
Int_t color_increm = 2;

Double_t size_text = 0.020;
//----- Automatons -----//
Int_t color_marker = color_line + color_increm;
Int_t color_marker_corr = color_line_corr + color_increm;

//----- Main -----//
if (!draw) gROOT->SetBatch();  // Do not draw plots to screen.

// Open the file and set the TTree. 
TFile* infile = TFile::Open(infile_path);

TTree* tree;
if (infile) {
    std::cout << "File opened successfully." << std::endl;
    tree = (TTree*)gDirectory->Get("tree");
}
else{
    std::cout << "[ERROR] File not opened." << std::endl; 
    return -1;
}
if (!tree) {
    cout << "[ERROR] Tree not opened." << endl;
    return -1;
}

// Bin the data in the TTree to show the integral.
Int_t n_tot = tree->GetEntries();
cout << "tree->GetEntries() = " << n_tot << endl;
tree->SetBranchStatus("m4mu",1);
tree->SetBranchStatus("m4mu_corr",1);
Float_t ptr_m4mu;
Float_t ptr_m4mu_corr;
tree->SetBranchAddress("m4mu", &ptr_m4mu);
tree->SetBranchAddress("m4mu_corr", &ptr_m4mu_corr);
TH1D* h_m4mu = new TH1D("h_m4mu", "h_m4mu", n_bins, m4mu_min, m4mu_max);
TH1D* h_m4mu_corr = new TH1D("h_m4mu_corr", "h_m4mu_corr", n_bins, m4mu_min, m4mu_max);
h_m4mu->Sumw2();
h_m4mu_corr->Sumw2();
for (int i = 0; i < n_tot; i++) {
// for (int i = 0; i < 5000; i++) {
    tree->GetEntry(i);
    h_m4mu->Fill(ptr_m4mu);
    h_m4mu_corr->Fill(ptr_m4mu_corr);
}
Int_t n_hist = h_m4mu->GetEntries();
cout << "number of entries in hist after filling: " << n_hist << endl;
// if (n_hist != n_tot) {
//     cout << "n_hist doesn't have the same entries as tree" << endl;
//     cout << "n_hist=" << n_hist << ", n_tree=" << n_tot << endl;
//     return -1;
// }
// if (n_hist != h_m4mu->Integral()) {
//     cout << "h_m4mu doesn't have the same entries as h->Integral()" << endl;
//     cout << "n_hist=" << n_hist << ", integral=" << h_m4mu->Integral() << endl;
//     return -1; 
// }

// DSCB fit.
RooRealVar m4mu("m4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
RooRealVar Mean("#mu", "#mu", 125, 120, 130);
RooRealVar Sigma("#sigma", "#sigma", 1, 0, 10);//sigma[decay]);
RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
RooRealVar ExpL("n_{L}", "n_{L}", 1, 0, 30);//expL[decay]);
RooRealVar ExpR("n_{R}", "n_{R}", 1, 1, 50);//expR[decay]);

RooRealVar m4mu_corr("m4mu_corr", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
RooRealVar Mean_corr("#mu^{corr}", "#mu^{corr}", 125, 120, 130);
RooRealVar Sigma_corr("#sigma^{corr}", "#sigma^{corr}", 1, 0, 10);//sigma[decay]);
RooRealVar AlphaL_corr("#alpha_{L}^{corr}", "#alpha_{L}^{corr}", 1, 0, 30);//alphaL[decay]);
RooRealVar ExpL_corr("n_{L}^{corr}", "n_{L}^{corr}", 1, 0, 30);//expL[decay]);
RooRealVar AlphaR_corr("#alpha_{R}^{corr}", "#alpha_{R}^{corr}", 1, 0, 30);//alphaR[decay]);
RooRealVar ExpR_corr("n_{R}^{corr}", "n_{R}^{corr}", 1, 1, 50);//expR[decay]);

RooDataSet rds("rds", "dataset with m4mu", tree, m4mu);
RooDataSet rds_corr("rds_corr", "dataset with m4mu_corr", tree, m4mu_corr);

RooMyPDF_DSCB DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
RooMyPDF_DSCB DSCB_corr("DSCB_corr", "DSCB_corr", m4mu_corr, Mean_corr, Sigma_corr, AlphaL_corr, ExpL_corr, AlphaR_corr, ExpR_corr);

TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 800, 600);
TPad* ptop = new TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0);
TPad* pbot = new TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25);
// c_MC->SetFrameFillColor(0);
// c_MC->cd(1)->SetBottomMargin(0.2);  // Moves bottom margin up to y_low = 0.2.
// c_MC->SetLogy();
// c_MC->Draw();
TString title1 = Form("%s MC %s, Unbinned Double-Sided CB Fit \n", fs.c_str(), year.c_str());
TString title2 = Form("(using p_{T} corr. factors derived from %s)", derive_from_sample.c_str());
TString title = title1 + title2;
RooPlot* xframe = m4mu.frame(RooFit::Title(title));
RooPlot* xframe_corr = m4mu_corr.frame(RooFit::Title(title));
RooPlot* framePull = m4mu.frame(RooFit::Title("")); // Ratio plot 
RooPlot* framePull_corr = m4mu_corr.frame(RooFit::Title("")); // Ratio plot

// Reisize y-axis.
Double_t max_val = max(h_m4mu->GetMaximum(), h_m4mu_corr->GetMaximum());
xframe->SetMaximum(max_val * 1.1);
xframe_corr->SetMaximum(max_val * 1.1);
xframe_corr->GetYaxis()->SetTitle(Form("Events / %.2f", h_m4mu->GetBinWidth(1)));

rds.plotOn(xframe, RooFit::MarkerColor(color_marker), RooFit::Binning(n_bins));
rds_corr.plotOn(xframe_corr, RooFit::MarkerColor(color_marker_corr), RooFit::Binning(n_bins));

DSCB.fitTo(rds, RooFit::Range(m4mu_min, m4mu_max), RooFit::PrintLevel(-1));
// DSCB.plotOn(xframe, RooFit::LineColor(color_line));
DSCB.plotOn(xframe, RooFit::LineColor(color_line), RooFit::Binning(n_bins));
DSCB.paramOn(xframe, RooFit::Layout(0.14, 0.34, 0.88));
xframe->getAttText()->SetTextSize(size_text);
xframe->getAttText()->SetTextColor(color_line);

DSCB_corr.fitTo(rds_corr, RooFit::Range(m4mu_min, m4mu_max), RooFit::PrintLevel(-1));
DSCB_corr.plotOn(xframe_corr, RooFit::LineColor(color_line_corr));
DSCB_corr.paramOn(xframe_corr, RooFit::Layout(0.66, 0.86, 0.88));
xframe_corr->getAttText()->SetTextSize(size_text);
xframe_corr->getAttText()->SetTextColor(color_line_corr);

// Get number of unbinned entries in TTree, making sure it passes m4l cut. 
// string cuts = Form("m4mu > %.1f && m4mu < %.1f", m4mu_min, m4mu_max);
// Double_t integral = tree->GetEntries(cuts.c_str());
Double_t integral = h_m4mu->Integral();
TString integ = Form("Integral = %.0f", integral);
TLatex* tex = new TLatex(0.12, 0.45, integ);
tex->SetNDC();
tex->SetTextColor(color_marker);

ptop->SetBottomMargin(0);
pbot->SetTopMargin(0);
pbot->SetBottomMargin(0.25);
ptop->Draw();
pbot->Draw();
ptop->cd();
ptop->SetTicks(1,1);
xframe->Draw();
tex->Draw("same");

// string cuts_corr = Form("m4mu_corr > %.1f && m4mu_corr < %.1f", m4mu_min, m4mu_max);
// Double_t integral_corr = tree->GetEntries(cuts_corr.c_str());
Double_t integral_corr = h_m4mu_corr->Integral();
TString integ_corr = Form("Integral = %.0f", integral_corr);
TLatex* tex_corr = new TLatex(0.65, 0.45, integ_corr);
tex_corr->SetNDC();
tex_corr->SetTextColor(color_marker_corr);
tex_corr->Draw("same");
xframe_corr->Draw("same"); 

// Add to ratio plots.
pbot->cd();
pbot->SetTicks(1,1);

// On the lower pad, either plot the residuals or make a hist ratio.
Float_t text_size_ratioplot = 0.10;
if (plot_residuals) {
    framePull->addObject((TObject*)xframe->pullHist(), "p");
    framePull_corr->addObject((TObject*)xframe_corr->pullHist(), "p");

    // TPad* padPull = new TPad("padPull", "padPull", 0., 0., 1., 0.2);
    // padPull->cd(2);
    // padPull->Draw();
    framePull->GetXaxis()->SetLabelSize(text_size_ratioplot);
    framePull->GetYaxis()->SetLabelSize(text_size_ratioplot);
    framePull->GetXaxis()->SetTitleSize(text_size_ratioplot);
    framePull->GetYaxis()->SetTitleSize(text_size_ratioplot);
    framePull->SetTitle("");
    framePull->GetYaxis()->SetTitle("Residuals");
    framePull->SetMinimum(-5.);
    framePull->SetMaximum(5.);
    framePull->SetNdivisions(207, "Y");
    framePull->SetTickLength(0.04, "XY");
    framePull->getAttMarker()->SetMarkerColor(color_marker);
    framePull->Draw();

    framePull_corr->GetXaxis()->SetLabelSize(text_size_ratioplot);
    framePull_corr->GetYaxis()->SetLabelSize(text_size_ratioplot);
    framePull_corr->GetXaxis()->SetTitleSize(text_size_ratioplot);
    framePull_corr->GetYaxis()->SetTitleSize(text_size_ratioplot);
    framePull_corr->SetTitle("");
    // framePull_corr->SetYTitle("Residuals");
    framePull_corr->GetYaxis()->SetTitle("Residuals");
    framePull_corr->SetMinimum(-5.);
    framePull_corr->SetMaximum(5.);
    framePull_corr->SetNdivisions(207, "Y");
    framePull_corr->SetTickLength(0.04, "XY");
    framePull_corr->getAttMarker()->SetMarkerColor(color_marker_corr);
    framePull_corr->Draw("same");
} 
// else {
//     ratio->GetXaxis()->SetLabelSize(text_size_ratioplot);
//     ratio->GetYaxis()->SetLabelSize(text_size_ratioplot);
//     ratio->GetXaxis()->SetTitleSize(text_size_ratioplot);
//     ratio->GetYaxis()->SetTitleSize(text_size_ratioplot);
//     ratio->GetYaxis()->SetTitleOffset(0.3);  // Default is 0.005.
//     ratio->SetXTitle("m_{4#mu} (GeV)");
//     ratio->SetYTitle("corrected / uncorr.");
//     ratio->SetMinimum(0.);
//     ratio->SetMaximum(2.);
//     ratio->SetNdivisions(207, "Y");
//     ratio->GetXaxis()->SetTickLength(0.12);
//     ratio->SetLineColor(kGreen+2);
//     gStyle->SetOptStat(0);
//     ratio->Draw();
// }

// Add a line at x = 1 onto ratio plot.
TLine* lineRef = new TLine(m4mu_min, 0., m4mu_max, 0.);
lineRef->Draw("same");

c_MC->Print(outfile_path);// + ".pdf");

return 0;
}
