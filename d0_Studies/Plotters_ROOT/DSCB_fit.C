#include "RooMyPDF_DSCB.h"
#include "RooRealVar.h"

// using namespace RooFit;

void DSCB_fit(Bool_t draw = false) {

if (!draw) gROOT->SetBatch();  // Do not draw plots to screen.

// Open the file and set the TTree. 
TString infilename = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/root_files/m4mu_withcorr_MC_2017_ggF_fullstats_pTlt200.root";
TString outpath_file = "/ufrc/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Plots/m4mu_DSCB_fits/test_04.pdf";
TFile* infile = TFile::Open(infilename);

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

// tree->SetBranchStatus("m4mu",1);
// tree->SetBranchStatus("m4mu_corr",1);

// Float_t m4mu; 
// Float_t m4mu_corr; 
// tree->SetBranchAddress("m4mu", &m4mu);
// tree->SetBranchAddress("m4mu_corr", &m4mu_corr);             

// Make m4mu hist.
// TH1F * h_m4mu = new TH1F("h_m4mu", "m_{4#mu} Distribution (Higgs Production Mode: ggF)", 140, 105, 140);

// Decide on number of events to run over.
// Int_t n_total = tree->GetEntries();
// Int_t n_loop;
// if (n == -1) {
//     n_loop = tree->GetEntries();
// }
// else {
//     n_loop = n;
// }

// cout << "...Looping over events..." << endl;
// for (int i = 0; i < n_loop; i++) {
//     tree->GetEntry(i);

//     h_m4mu->Fill(m4mu);
//     cout << "Event " << i << endl;
//     cout << "The value of m4mu is " << m4mu << endl;
//     cout << "The value of m4mu_corr is " << m4mu_corr << endl;
// }

// DSCB fit.
RooRealVar m4mu("m4mu", "m_{4#mu}", 105, 140, "GeV");
RooRealVar m4mu_corr("m4mu_corr", "m_{4#mu}", 105, 140, "GeV");

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
c_MC->SetFrameFillColor(0);
// c_MC->cd(1)->SetBottomMargin(0.2);
// 	c_MC->SetLogy();

TString title = "MC 2017, Double-Sided CB Fit (Unbinned)";
RooPlot* xframe = m4mu.frame(RooFit::Title(title));
RooPlot* xframe_corr = m4mu_corr.frame(RooFit::Title(title));
rds.plotOn(xframe);
rds_corr.plotOn(xframe_corr);

Int_t color = kBlue;
Double_t size_text = 0.020;
DSCB.fitTo(rds, RooFit::Range(105, 140));
DSCB.plotOn(xframe, RooFit::LineColor(color));
DSCB.paramOn(xframe, RooFit::Layout(0.15, 0.35, 0.90));
xframe->getAttText()->SetTextSize(size_text);
xframe->getAttText()->SetTextColor(color);

Int_t color_corr = kRed;
DSCB_corr.fitTo(rds_corr, RooFit::Range(105, 140));
DSCB_corr.plotOn(xframe_corr, RooFit::LineColor(color_corr));
// RooPlot *framePull_DATA = x.frame("");
// framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
// DSCB.plotOn(xframe, Components("bkg"), LineColor(kBlue), LineStyle(kDashed));
DSCB_corr.paramOn(xframe_corr, RooFit::Layout(0.7, 0.90, 0.90));
xframe_corr->getAttText()->SetTextSize(size_text);
xframe_corr->getAttText()->SetTextColor(color_corr);

xframe->Draw();
xframe_corr->Draw("same"); 



// TPad* padPull_DATA =  new  TPad("padPull","padPull",0.,0.,1.,0.2);
// padPull_DATA->Draw();
// padPull_DATA->cd(0);
// framePull_DATA->GetYaxis()->SetLabelSize(0.1);
// framePull_DATA->GetXaxis()->SetLabelSize(0.1);
// framePull_DATA->SetMinimum(-5.);
// framePull_DATA->SetMaximum(5.);
// framePull_DATA->Draw();
// TLine* lineRef = new TLine(105,0,140,0.);
// lineRef->Draw("same");

c_MC->Print(outpath_file);// + ".pdf");

// param.push_back(Mean.getVal());
// param.push_back(Sigma.getVal());

// Draw the plots.
// TCanvas * c1 = new TCanvas("c1", "c1", 600, 700);
// h_m4mu->Draw();
// c1->SaveAs("/ufrc/avery/rosedj1/HiggsMassMeasurement/Utils_ROOT/TestingROOT/test_draw.pdf");

return 0;
}