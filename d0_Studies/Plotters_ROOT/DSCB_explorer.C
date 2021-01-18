#include "RooMyPDF_DSCB.h"
#include "RooArgSet.h"
#include <vector>

// Functions.
void skim_ggH_file(TString infile_path, TString outfile_path);
// vector<double> linspace(double start, double stop, int num); 
// RooFitResult*  DSCB_unbinnedfit(TTree* tree, Double_t m4mu_min, Double_t m4mu_max);
// void           skim_ggH_file(TString infile_path, TString outfile_path);
// TTree*         sort_tree(TTree* tree);

vector<double> linspace(double start, double stop, int num) {
    int steps = num - 1;
    vector<double> pts;
    // Width of each segment.
    double width = (stop - start) / steps;
    pts.push_back(start);
    // Create the intermediate points, between start and stop.
    for (int k = 1; k < steps; k++) {
        pts.push_back(start + k*width);
    }
    // Finally add the end point.
    pts.push_back(stop);
    return pts;
}

RooMyPDF_DSCB* prep_dscb(RooRealVar& mass4mu) {
    /* Return a DSCB pdf of the m4mu values stored in tree. */
    //    RooRealVar mass4mu("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
    //    RooRealVar passedFullSelection("passedFullSelection","passedFullSelection", 0, 1);
    //    RooRealVar finalState("finalState","finalState", 0, 1);

    //    RooRealVar& x = mass4mu;
    //    cout << "x is a reference to mass4mu" << endl;
    //    cout << "Address of x=      " << &x << endl;
    //    cout << "Address of mass4mu=" << &mass4mu << endl;
    //    cout << "Value of x=      " << x << endl;
    //    cout << "Value of mass4mu=" << mass4mu << endl;
    RooRealVar Mean("#mu", "#mu", 125, 120, 130);
    RooRealVar Sigma("#sigma", "#sigma", 1, 0, 10);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 1, 0, 30);//expL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
    RooRealVar ExpR("n_{R}", "n_{R}", 1, 1, 50);//expR[decay]);
    RooMyPDF_DSCB* dscb = new RooMyPDF_DSCB("DSCB", "DSCB", mass4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    cout << "Test 1" << endl;
    return dscb;
}

RooDataSet prep_roodataset(TTree* tree, RooRealVar& mass4mu,
                           Double_t m4mu_min, Double_t m4mu_max, RooMyPDF_DSCB& DSCB) {
    /* Return a RooDataSet of the mass4mu values in the specified [min,max]. */
    TString cuts = Form("(%f < mass4mu) && (mass4mu < %f)", m4mu_min, m4mu_max);
    cout << "Applying cuts: " << cuts << endl;
    //    const char* cuts = "passedFullSelection==1 && finalState==1 && (mass4mu";
    //    RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu, passedFullSelection, finalState), cuts);
    RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu), cuts);
    //    RooDataSet* rds = new RooDataSet("rds", "roodataset m4mu", tree, RooArgSet(mass4mu), cuts);
    return rds;
}

RooFitResult* DSCB_unbinnedfit(
    TTree* tree,
    Float_t m4mu_min, Float_t m4mu_max,
    Float_t dm4mu_min, Float_t dm4mu_max,
    TCanvas* canv, TString outfile_path) {
    // Double_t m4mu_min, Double_t m4mu_max,
    // Double_t m4muErr_min, Double_t m4muErr_max
    /* 
     * Do an unbinned DSCB fit on the mass4mu data stored in `tree`.
     * Apply certain selections in the RooDataSet.
     * Return the fit result.
     */
    if (!tree) {
        cout << "[ERROR] Tree not opened." << endl;
        // Throw exception here.
    }
    // Prepare the DSCB variables.
    RooRealVar mass4mu("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
    RooRealVar Mean("mean", "#mu", 125, 120, 130);
    RooRealVar Sigma("sigma", "#sigma", 1, 0, 10);//sigma[decay]);
    RooRealVar AlphaL("alphaL", "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
    RooRealVar ExpL("nL", "n_{L}", 1, 0, 30);//expL[decay]);
    RooRealVar AlphaR("alphaR", "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
    RooRealVar ExpR("nR", "n_{R}", 1, 1, 50);//expR[decay]);
    RooMyPDF_DSCB DSCB("DSCB", "DSCB", mass4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // Prepare the cuts.
    RooRealVar passedFullSelection("passedFullSelection","passedFullSelection", 0, 1);
    RooRealVar finalState("finalState","finalState", 0, 1);
    RooRealVar mass4lErr("mass4lErr","mass4lErr", 0, 1000, "GeV");
    TString cuts1 = "passedFullSelection==1 && finalState==1";
    TString cuts2 = Form("%f < mass4mu && mass4mu < %f", m4mu_min, m4mu_max);
    TString cuts3 = Form("%f < mass4lErr && mass4lErr < %f", m4muErr_min, m4muErr_max);
    TString cuts = cuts1 + " && " + cuts2 + " && " + cuts3
    cout << "Applying the following cuts:" << endl;
    cout << "  " << cuts << endl;
    RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu, mass4lErr, passedFullSelection, finalState), cuts);
    // RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu), cuts);
    // Do the DSCB fit.
        // Now 
    RooFitResult* result = DSCB.fitTo(rds,
                                 RooFit::Range(m4mu_min, m4mu_max),
                                 RooFit::PrintLevel(-1),
                                 RooFit::Save()
                                 );

    TString fs = "ggH";
    TString year = "2018";
    Float_t size_text = 0.025;
    // New code.
    // TPad* ptop = new TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0);
    // TPad* pbot = new TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25);
    TString title = Form("%s MC %s, Unbinned DSCB Fit \n", fs.c_str(), year.c_str());

    RooPlot* xframe = mass4mu.frame(RooFit::Title(title));
    RooPlot* framePull = mass4mu.frame(RooFit::Title("")); // Ratio plot
    // Make binned version of data.
    TH1F* h_m4mu = rds.createHistogram("rds", 70);  // Name of data, n_bins.
    Double_t max_val = h_m4mu->GetMaximum();
    xframe->SetMaximum(max_val * 1.1);
    xframe->GetYaxis()->SetTitle(Form("Events / %.2f", h_m4mu->GetBinWidth(1)));
    // xframe->addTH1(h_m4mu)
    // When the RooDataSet plots on the frame, it will look binned.
    rds.plotOn(xframe, RooFit::MarkerColor(1), RooFit::Binning(n_bins));
    // DSCB.plotOn(xframe, RooFit::LineColor(color_line));
    DSCB.plotOn(xframe, RooFit::LineColor(2), RooFit::Binning(n_bins));
    DSCB.paramOn(xframe, RooFit::Layout(0.14, 0.34, 0.88));
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(2);

    Double_t integral = h_m4mu->Integral();
    TString integ = Form("Integral = %.0f", integral);
    TLatex* tex = new TLatex(0.12, 0.45, integ);
    tex->SetNDC();
    tex->SetTextColor(1);

    // ptop->SetBottomMargin(0);
    // pbot->SetTopMargin(0);
    // pbot->SetBottomMargin(0.25);
    // ptop->Draw();
    // pbot->Draw();
    // ptop->cd();
    // ptop->SetTicks(1,1);
    xframe->Draw();
    tex->Draw("same");

    // Add to ratio plots.
    // pbot->cd();
    // pbot->SetTicks(1,1);

    // Float_t text_size_ratioplot = 0.10;
    // framePull->addObject((TObject*)xframe->pullHist(), "p");
    // framePull_corr->addObject((TObject*)xframe_corr->pullHist(), "p");

    // TPad* padPull = new TPad("padPull", "padPull", 0., 0., 1., 0.2);
    // padPull->cd(2);
    // padPull->Draw();

    // framePull->GetXaxis()->SetLabelSize(text_size_ratioplot);
    // framePull->GetYaxis()->SetLabelSize(text_size_ratioplot);
    // framePull->GetXaxis()->SetTitleSize(text_size_ratioplot);
    // framePull->GetYaxis()->SetTitleSize(text_size_ratioplot);
    // framePull->SetTitle("");
    // framePull->GetYaxis()->SetTitle("Residuals");
    // framePull->SetMinimum(-5.);
    // framePull->SetMaximum(5.);
    // framePull->SetNdivisions(207, "Y");
    // framePull->SetTickLength(0.04, "XY");
    // framePull->getAttMarker()->SetMarkerColor(color_marker);
    // framePull->Draw();

    // framePull_corr->GetXaxis()->SetLabelSize(text_size_ratioplot);
    // framePull_corr->GetYaxis()->SetLabelSize(text_size_ratioplot);
    // framePull_corr->GetXaxis()->SetTitleSize(text_size_ratioplot);
    // framePull_corr->GetYaxis()->SetTitleSize(text_size_ratioplot);
    // framePull_corr->SetTitle("");
    // // framePull_corr->SetYTitle("Residuals");
    // framePull_corr->GetYaxis()->SetTitle("Residuals");
    // framePull_corr->SetMinimum(-5.);
    // framePull_corr->SetMaximum(5.);
    // framePull_corr->SetNdivisions(207, "Y");
    // framePull_corr->SetTickLength(0.04, "XY");
    // framePull_corr->getAttMarker()->SetMarkerColor(color_marker_corr);
    // framePull_corr->Draw("same");

    // Add a line at x = 1 onto ratio plot.
    // TLine* lineRef = new TLine(m4mu_min, 0., m4mu_max, 0.);
    // lineRef->Draw("same");

    canv->Print(outfile_path);// + ".pdf");

    return result;
}

void skim_ggH_file(TString infile_path, TString outfile_path) {
    /* Opens a ggH TTree stored inside a root file (infile_path).
     * The TTree is cloned. 
     * The following event selections are made:
     * * 105 < mass4mu < 140
     * * passedFullSelection == 1
     * * finalState == 1
     * Saves a '.root' file with the skimmed TTree inside.
     */
    TFile* f = TFile::Open(infile_path);
    TTree* t = (TTree*)f->Get("passedEvents");
    t->SetBranchStatus("*", 1);
    // Get mass4mu values which passFullSelection, finalState, and largest
    // mass4mu cut.
    Float_t mass4mu;
    Int_t finalState;
    Bool_t passedFullSelection;
    t->SetBranchAddress("mass4mu", &mass4mu);
    t->SetBranchAddress("finalState", &finalState);
    t->SetBranchAddress("passedFullSelection", &passedFullSelection);

    TFile newfile(outfile_path, "recreate");
    TTree* skimtree = t->CloneTree(0);

    bool pass_m4mu;
    bool pass_fullselect;
    bool pass_fs;
    Int_t n_evts = t->GetEntries();
    // Fill the new skimtree with only the events which pass selection.
    for (auto k : ROOT::TSeqI(n_evts)) {
        t->GetEntry(k);
        pass_m4mu = (105 < mass4mu) && (mass4mu < 140);
        pass_fullselect = (passedFullSelection == 1);
        pass_fs = (finalState == 1);
        if (pass_m4mu && pass_fullselect && pass_fs)
            skimtree->Fill();
    }
    cout << "skim tree has: " << skimtree->GetEntries() << endl;

    skimtree->Print();
    newfile.Write();

    return;
}

TTree* sort_tree(TTree* tree) {
    //     /*
    //     Return a TTree pointer with sorted mass4mu values
    //     originally found in tree.
    //     NOTES: Only returns one branch of tree!
    //     */
    //     vector<Float_t> vec;

    //     // Retrieve values from TTree.
    //     tree->SetBranchStatus("*", 0);
    //     tree->SetBranchStatus("mass4mu", 1);
    //     tree->SetBranchStatus("mass4lErr", 1);
    //     tree->SetBranchStatus("finalState", 1);
    //     tree->SetBranchStatus("passedFullSelection", 1);
    //     Float_t mass4mu;
    //     Float_t mass4lErr;
    //     Int_t finalState;
    //     Bool_t passedFullSelection;
    //     tree->SetBranchAddress("mass4mu", &mass4mu);
    //     tree->SetBranchAddress("mass4lErr", &mass4lErr);
    //     tree->SetBranchAddress("finalState", &finalState);
    //     tree->SetBranchAddress("passedFullSelection", &passedFullSelection);
        
    //     // Fill vector.
    //     Int_t n = tree->GetEntries();
    //     for (auto k : ROOT::TSeqI(n)) {
    //         tree->GetEntry(k);
    //         vec.push_back(m4mu);
    //     }
    //     cout << "...Sorting the vector..." << endl;
    //     sort(vec.begin(), vec.end());
    //     int ctr = 0;
    //     cout << "...Making sorted TTree..." << endl;
    // //    TFile newfile("deleteme.root", "recreate"); // For memory purposes.
    //     TTree* sortedtree = tree->CloneTree(0);
    //     Float_t m4mu_sorted;
    //     sortedtree->SetBranchAddress("mass4mu", &m4mu_sorted);
    //     // Fill it with the sorted vector elements.
    //     for (const auto& value: vec) { 
    //         m4mu_sorted = value;
    //         sortedtree->Fill();
    //     }
    //     return sortedtree;
}

TTree* get_ttree_from_liteskim(
    TString infile_path,
    Float_t &mass4mu,
    Float_t &mass4lErr,
    Int_t &finalState,
    Bool_t &passedFullSelection) {
    /* Open a TTree in infile and activate a few branches. */
    //   Double_t m4mu_min = 121.2;
    //   Double_t m4mu_max = 128.4;
    TFile* f = TFile::Open(infile_path);
    TTree* tree = (TTree*)f->Get("passedEvents");

    // Turn on appropriate branches.
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("mass4mu", 1);
    tree->SetBranchStatus("mass4lErr", 1);
    tree->SetBranchStatus("finalState", 1);
    tree->SetBranchStatus("passedFullSelection", 1);
    // Link variable address to branch:
    tree->SetBranchAddress("mass4mu", &mass4mu);
    tree->SetBranchAddress("mass4lErr", &mass4lErr);
    tree->SetBranchAddress("finalState", &finalState);
    tree->SetBranchAddress("passedFullSelection", &passedFullSelection);
    return tree;
}

void draw_ttree(TTree* tree, Int_t n_bins, Float_t m4mu_min, Float_t m4mu_max) {
        // Make 1 DSCB plot.
        TH1D* h_m4mu = new TH1D("h_m4mu", "h_m4mu", n_bins, m4mu_min, m4mu_max);
        h_m4mu->Sumw2();
        TString cuts1 = Form("(%f < mass4mu) && (mass4mu < %f)", m4mu_min, m4mu_max);
        TString cuts2 = "passedFullSelection==1 && finalState==1";
        TString cuts = cuts1 + " && " + cuts2;
        cout << "Applying cuts: " << cuts << endl;
        gROOT->SetBatch(kTRUE);
        TCanvas* c = new TCanvas("canv", "Test canv", 600, 620);
        TString outfile_path = "/afs/cern.ch/work/d/drosenzw/Higgs/HiggsMassMeasurement/d0_Studies/output/DSCBscanoutput/test02.pdf";
        c->Print(outfile_path + "[");
        // h_m4mu->Draw();
        tree->Draw("mass4mu", cuts1);
        c->Print(outfile_path);
        tree->Draw("mass4mu", cuts2);
        c->Print(outfile_path);
        tree->Draw("mass4mu", cuts);
        c->Print(outfile_path + "]");
}

void DSCB_explorer() {
//    skim_ggH_file()
    TString infile_path = "/eos/home-d/drosenzw/Samples/LiteSkim/MC2018ggH_liteskim_fullstats.root";
    Double_t m4mu_min = 105.0;
    Double_t m4mu_max = 140.0;
    Int_t n_bins = 100;

    // Branch variables.
    Float_t mass4mu;
    Float_t mass4lErr;
    Int_t finalState;
    Bool_t passedFullSelection;

    TTree* tree = get_ttree_from_liteskim(infile_path, mass4mu, mass4lErr, finalState, passedFullSelection);
    draw_ttree(tree, n_bins, m4mu_min, m4mu_max);
    return;

//     RooRealVar mass4mu("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
// //    RooRealVar Mean("#mu", "#mu", 125, 120, 130);
// //    RooRealVar Sigma("#sigma", "#sigma", 1, 0, 10);//sigma[decay]);
// //    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
// //    RooRealVar ExpL("n_{L}", "n_{L}", 1, 0, 30);//expL[decay]);
// //    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
// //    RooRealVar ExpR("n_{R}", "n_{R}", 1, 1, 50);//expR[decay]);
// //    RooMyPDF_DSCB DSCB("DSCB", "DSCB", mass4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
//     cout << "Test 1" << endl;

//     TCanvas *canv = new TCanvas("canv", "canv", 600, 620);


// //    RooRealVar* mass4mu = new RooRealVar("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
//     RooMyPDF_DSCB* pDSCB = prep_dscb(mass4mu);
// //    RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu, passedFullSelection, finalState), cuts);
// //    RooDataSet* rds = new RooDataSet("rds", "roodataset m4mu", tree, RooArgSet(mass4mu), cuts);
//     RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu), cuts);
//     cout << "Test 2" << endl;
// //    RooDataSet rds = prep_roodataset(tree, mass4mu, m4mu_min, m4mu_max, DSCB);
//     RooFitResult* result = DSCB_unbinnedfit(*pDSCB, rds);
// //    RooFitResult* result = DSCB.fitTo(rds,
// //                                 RooFit::PrintLevel(-1),
// //                                 RooFit::Save()
// //                                 );
//     cout << "Test 3" << endl;
//     result->Print();

//     return;
}

// void DSCB_explorer_BROKEN() {
//     /* HAVING TROUBLE WITH MEMORY ALLOCATION */
// //    skim_ggH_file()
//     TString infile_path = "skimtree_MC2018ggH_sorted.root";
//     Double_t m4mu_min = 105.0;
//     Double_t m4mu_max = 140.0;
//  //   Double_t m4mu_min = 121.2;
//  //   Double_t m4mu_max = 128.4;
// 
//     TFile* f = TFile::Open(infile_path);
//     TTree* tree = (TTree*)f->Get("passedEvents");
//     tree->SetBranchStatus("*", 1);
// 
//     RooRealVar mass4mu("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
//     // Make the corresponding pointer, which stores address of mass4mu.
// //    RooRealVar* m4mu = &mass4mu;
// //    RooRealVar* mass4mu = new RooRealVar("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
//     RooMyPDF_DSCB DSCB = prep_dscb(mass4mu);
//     cout << "Test 2" << endl;
//     RooDataSet rds = prep_roodataset(tree, mass4mu, m4mu_min, m4mu_max, DSCB);
//     cout << "Test 3" << endl;
// //    RooFitResult* result = DSCB_unbinnedfit(DSCB, rds);
//     RooFitResult* result = DSCB.fitTo(rds,
// //                                 RooFit::Range(m4mu_min, m4mu_max),
//                                  RooFit::PrintLevel(-1),
//                                  RooFit::Save()
//                                  );
//     cout << "Test 4" << endl;
//     result->Print();
// //    TTree* sortedtree = skimtree->CloneTree(0);
// //    sortedtree = sort_tree(skimtree);
// //    sortedtree->Print();
// //    newfile.Write();
// 
//     return;
// }
