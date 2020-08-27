#include "RooMyPDF_DSCB.h"
#include "RooRealVar.h"

// using namespace RooFit;

void DSCB_fit_synchwithFilippo(Bool_t draw = false) {

//----- User Params -----//
Double_t m4mu_min = 105;
Double_t m4mu_max = 140;
string year = "2018";

// This is the color of the fits.
// The data points will be darker (+2).
Int_t color_line = kBlack;
Int_t color_marker = color_line; // 15;

Int_t color_line_corr = kGreen+2;
Int_t color_marker_corr = color_line_corr; // kGreen-3;

//----- Main -----//
if (!draw) gROOT->SetBatch();  // Do not draw plots to screen.

// Open the file and set the TTree. 
TString infile_path = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/MC2018_ggF_synchwithFilippo_basiccuts_usingFSR.root";
TString outfile_path = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plots/m4mu_DSCB_fits/20200713_MC2018_ggF_synchwithFilippo_basiccuts_usingFSR_DSCBfits_final.pdf";
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
TH1D* h_m4mu = new TH1D("h_m4mu", "h_m4mu", 100, m4mu_min, m4mu_max);
TH1D* h_m4mu_corr = new TH1D("h_m4mu_corr", "h_m4mu_corr", 100, m4mu_min, m4mu_max);
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
RooRealVar m4mu_corr("m4mu_corr", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
RooDataSet rds("rds", "dataset with m4mu", tree, m4mu);
RooDataSet rds_corr("rds_corr", "dataset with m4mu_corr", tree, m4mu_corr);

RooRealVar Mean("#mu", "#mu", 125, 120, 130);
RooRealVar Sigma("#sigma", "#sigma", 1, 0, 10);//sigma[decay]);
RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
RooRealVar ExpL("n_{L}", "n_{L}", 1, 0, 30);//expL[decay]);
RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
RooRealVar ExpR("n_{R}", "n_{R}", 1, 1, 50);//expR[decay]);

RooRealVar Mean_corr("#mu^{corr}", "#mu^{corr}", 125, 120, 130);
RooRealVar Sigma_corr("#sigma^{corr}", "#sigma^{corr}", 1, 0, 10);//sigma[decay]);
RooRealVar AlphaL_corr("#alpha_{L}^{corr}", "#alpha_{L}^{corr}", 1, 0, 30);//alphaL[decay]);
RooRealVar ExpL_corr("n_{L}^{corr}", "n_{L}^{corr}", 1, 0, 30);//expL[decay]);
RooRealVar AlphaR_corr("#alpha_{R}^{corr}", "#alpha_{R}^{corr}", 1, 0, 30);//alphaR[decay]);
RooRealVar ExpR_corr("n_{R}^{corr}", "n_{R}^{corr}", 1, 1, 50);//expR[decay]);
RooMyPDF_DSCB DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
RooMyPDF_DSCB DSCB_corr("DSCB_corr", "DSCB_corr", m4mu_corr, Mean_corr, Sigma_corr, AlphaL_corr, ExpL_corr, AlphaR_corr, ExpR_corr);

TCanvas *c_MC = new TCanvas("c_MC", "c_MC", 800, 600);
TPad* ptop = new TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0);
TPad* pbot = new TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25);
// c_MC->SetFrameFillColor(0);
// c_MC->cd(1)->SetBottomMargin(0.2);  // Moves bottom margin up to y_low = 0.2.
// c_MC->SetLogy();

TString title = Form("MC %s, Double-Sided CB Fit (Unbinned)", year.c_str());
RooPlot* xframe = m4mu.frame(RooFit::Title(title));
RooPlot* xframe_corr = m4mu_corr.frame(RooFit::Title(title));

rds.plotOn(xframe, RooFit::MarkerColor(color_marker));
rds_corr.plotOn(xframe_corr, RooFit::MarkerColor(color_marker_corr));
// Reisize y-axis.
Double_t max_val = max(h_m4mu->GetMaximum(), h_m4mu_corr->GetMaximum());
xframe->SetMaximum(max_val * 1.1);

Double_t size_text = 0.020;
DSCB.fitTo(rds, RooFit::Range(m4mu_min, m4mu_max), RooFit::PrintLevel(-1));
DSCB.plotOn(xframe, RooFit::LineColor(color_line));
DSCB.paramOn(xframe, RooFit::Layout(0.14, 0.34, 0.7));
xframe->getAttText()->SetTextSize(size_text);
xframe->getAttText()->SetTextColor(color_line);

DSCB_corr.fitTo(rds_corr, RooFit::Range(m4mu_min, m4mu_max), RooFit::PrintLevel(-1));
DSCB_corr.plotOn(xframe_corr, RooFit::LineColor(color_line_corr));
DSCB_corr.paramOn(xframe_corr, RooFit::Layout(0.66, 0.86, 0.7));
xframe_corr->getAttText()->SetTextSize(size_text);
xframe_corr->getAttText()->SetTextColor(color_line_corr);

// Get number of unbinned entries in TTree, making sure it passes m4l cut. 
// string cuts = Form("m4mu > %.1f && m4mu < %.1f", m4mu_min, m4mu_max);
// Double_t integral = tree->GetEntries(cuts.c_str());
// Double_t integral = h_m4mu->Integral();
// TString integ = Form("Integral = %.1f", integral);
// TLatex* tex = new TLatex(0.12, 0.45, integ);
// tex->SetNDC();
// tex->SetTextColor(color_marker);

ptop->SetBottomMargin(0);
pbot->SetTopMargin(0);
pbot->SetBottomMargin(0.25);
ptop->Draw();
pbot->Draw();
ptop->cd();
ptop->SetTicks(1,1);
xframe->Draw();
// tex->Draw("same");

// string cuts_corr = Form("m4mu_corr > %.1f && m4mu_corr < %.1f", m4mu_min, m4mu_max);
// Double_t integral_corr = tree->GetEntries(cuts_corr.c_str());

// Double_t integral_corr = h_m4mu_corr->Integral();
// TString integ_corr = Form("Integral = %.1f", integral_corr);
// TLatex* tex_corr = new TLatex(0.65, 0.45, integ_corr);
// tex_corr->SetNDC();
// tex_corr->SetTextColor(color_marker_corr);
// tex_corr->Draw("same");
xframe_corr->Draw("same"); 

// Add a legend.
// FIXME: Do not know how to get correct DSCB line/point in TLegend.
// ptop->cd();
// TLegend* leg = new TLegend(0.14, 0.3, 0.29, 0.45);
// leg->AddEntry(DSCB, "Uncorrected", "lep");
// leg->AddEntry(DSCB_corr, "Corrected", "lep");
// leg->AddEntry("DSCB", "Uncorrected"); //, "lep");
// leg->AddEntry("DSCB_corr", "Corrected"); //, "lep");
// leg->AddEntry(xframe, "Uncorrected", "lep");
// leg->AddEntry(xframe_corr, "Corrected", "lep");
// leg->Draw();
// ptop->Update();


// Add ratio plots.
pbot->cd();
pbot->SetTicks(1,1);
TH1F* ratio = new TH1F("h_ratio", "", 100, m4mu_min, m4mu_max);
ratio->Divide(h_m4mu_corr, h_m4mu);
// ratio->addObject((TObject*)xframe->pullHist(), "p");

// TPad* padPull = new TPad("padPull", "padPull", 0., 0., 1., 0.2);
// padPull->cd(2);
// padPull->Draw();
Float_t text_size_ratioplot = 0.10;
ratio->GetXaxis()->SetLabelSize(text_size_ratioplot);
ratio->GetYaxis()->SetLabelSize(text_size_ratioplot);
ratio->GetXaxis()->SetTitleSize(text_size_ratioplot);
ratio->GetYaxis()->SetTitleSize(text_size_ratioplot);
ratio->GetYaxis()->SetTitleOffset(0.3);  // Default is 0.005.
ratio->SetXTitle("m_{4#mu} (GeV)");
ratio->SetYTitle("corrected / uncorr.");
ratio->SetMinimum(0.);
ratio->SetMaximum(2.);
ratio->SetNdivisions(207, "Y");
ratio->GetXaxis()->SetTickLength(0.12);
ratio->SetLineColor(kGreen+2);
gStyle->SetOptStat(0);
ratio->Draw();

// Add a line at x = 1 onto ratio plot.
TLine* lineRef = new TLine(105.,1.,140.,1.);
lineRef->Draw("same");

c_MC->Print(outfile_path);// + ".pdf");

return 0;
}