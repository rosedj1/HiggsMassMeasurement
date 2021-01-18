RooFitResult* fit_and_draw_DSCB(
    TTree* tree,
    Float_t m4mu_min, Float_t m4mu_max,
    Float_t m4muErr_min, Float_t m4muErr_max,  // Passed in as a percent.
    vector<Double_t> &integral_vec, TCanvas* canv, TString outfile_path, Int_t n_bins=100) {
    // Double_t m4mu_min, Double_t m4mu_max,
    // Double_t m4muErr_min, Double_t m4muErr_max
    /** 
     * Do an unbinned DSCB fit on the mass4mu data stored in `tree`.
     * Apply certain selections in the RooDataSet.
     * Return the fit result.

     Parameters
     ----------
     n_bins : Number of bins to be places across [105, 140] GeV
     */
    if (!tree) {
        cout << "[ERROR] Tree not opened." << endl;
        // Throw exception here.
    }
    // Prepare the DSCB variables.
    // RooRealVar mass4mu("mass4mu", "m_{4#mu}", m4mu_min, m4mu_max, "GeV");
    RooRealVar mass4mu("mass4mu", "m_{4#mu}", 105, 140, "GeV");
    RooRealVar Mean("#mu", "#mu", 125, 120, 130);
    RooRealVar Sigma("#sigma", "#sigma", 1, 0, 10);//sigma[decay]);
    RooRealVar AlphaL("#alpha_{L}", "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
    RooRealVar AlphaR("#alpha_{R}", "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
    RooRealVar ExpL("n_{L}", "n_{L}", 1, 0, 50);//expL[decay]);
    /**
     * The right tail of the m(4mu) resonance decays very quickly.
     * This means that n_R (related to the exponential decay of the tail)
     * will typically be O(100), but capping it at 100 doesn't change
     * the other fit params much based on preliminary testing.
     * So leave it at 100 for faster fitting.
     */ 
    RooRealVar ExpR("n_{R}", "n_{R}", 1, 0, 100);//expR[decay]); (1,1,50)
    RooMyPDF_DSCB DSCB("DSCB", "DSCB", mass4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // Prepare the cuts.
    RooRealVar passedFullSelection("passedFullSelection","passedFullSelection", 0, 1);
    RooRealVar finalState("finalState","finalState", 0, 1);
    RooRealVar mass4lErr("mass4lErr","mass4lErr", 0, 1000, "GeV");
    TString sel_passfulls_and_fs = "passedFullSelection==1 && finalState==1";
    TString sel_mass4mu = Form("%f < mass4mu && mass4mu < %f", m4mu_min, m4mu_max);
    TString sel_relmasserr = Form("%f < mass4lErr/mass4mu*100.0 && mass4lErr/mass4mu*100.0 < %f", m4muErr_min, m4muErr_max);
    // TString cuts = cuts1 + " && " + cuts2 + " && " + cuts3;
    TString cuts_wo_mass4mu = sel_passfulls_and_fs + " && " + sel_relmasserr;
    TString all_cuts = cuts_wo_mass4mu + " && " + sel_mass4mu;
    cout << "***** Fitting a DSCB to mass4mu with the following cuts: *****" << endl;
    cout << "  " << cuts_wo_mass4mu << endl;
    RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(mass4mu, mass4lErr, passedFullSelection, finalState), cuts_wo_mass4mu);
    // Make a RooDataSet to plot black points along full mass4mu range:
    RooDataSet rds_withmass4mucuts("rds_withmass4mucuts", "roodataset m4mu", tree, RooArgSet(mass4mu, mass4lErr, passedFullSelection, finalState), all_cuts);
    // Do the DSCB fit.
    gStyle->SetOptStat("iouRMe");
    RooFitResult* result = DSCB.fitTo(rds,
                                 RooFit::Range(m4mu_min, m4mu_max),
                                 RooFit::PrintLevel(-1),
                                 RooFit::Save()
                                 );

    const char* fs = "ggH";
    const char* year = "2018";
    Float_t size_text = 0.02;
    // New code.
    // TPad* ptop = new TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0);
    // TPad* pbot = new TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25);
    TString title_p1 = Form("%s MC %s DSCB Unbinned Fit, ", fs, year);
    // Need extra spaces in title all of a sudden...
    TString title_p2 = Form("%.2f <   #deltam_{4#mu}/m_{4#mu} < %.2f%%", m4muErr_min, m4muErr_max);
    TString title = title_p1 + title_p2;
    RooPlot* xframe = mass4mu.frame(RooFit::Range(105, 140), RooFit::Title(title));
    // RooPlot* framePull = mass4mu.frame(RooFit::Title("")); // Ratio plot
    // Make binned version of data.
    TString hist_name = Form("h_m4mu_%f_%f", m4mu_min, m4mu_max);
    TH1* h_m4mu = rds.createHistogram(hist_name, mass4mu, RooFit::Binning(n_bins, 105, 140));
    h_m4mu->SetLineWidth(0);
    // Double_t max_val = h_m4mu->GetMaximum();
    // xframe->SetMaximum(max_val * 1.1);
    // xframe->addTH1(h_m4mu)
    // When the RooDataSet plots on the frame, it will look binned.
    Double_t bin_width = h_m4mu->GetBinWidth(1);
    xframe->GetYaxis()->SetTitle(Form("Events / (%.2f GeV)", bin_width));
    // Int_t n_bins_in_fitrange = static_cast<Int_t>(round((m4mu_max - m4mu_min) / bin_width));
    // rds.plotOn(xframe, RooFit::MarkerColor(1), RooFit::Binning(n_bins_in_fitrange));
    rds.plotOn(xframe, RooFit::MarkerColor(1), RooFit::Binning(n_bins));
    // DSCB.plotOn(xframe, RooFit::LineColor(color_line));
    DSCB.plotOn(xframe, RooFit::LineColor(2), RooFit::LineWidth(3), RooFit::Range(m4mu_min, m4mu_max));  //RooFit::Binning(n_bins), 
    DSCB.paramOn(xframe, RooFit::Layout(0.20, 0.40, 0.80));
    xframe->getAttText()->SetTextSize(size_text);
    xframe->getAttText()->SetTextColor(2);
    // h_m4mu->SetTitle(title);
    // h_m4mu->Draw("same");  // For the stats box.
    xframe->Draw();

    TLatex* tex_fitparams = new TLatex(0.20, 0.82, "DSCB Fit Parameters:");
    tex_fitparams->SetNDC();
    tex_fitparams->SetTextColor(2);
    tex_fitparams->SetTextSize(size_text);
    tex_fitparams->Draw("same");
    
    // Evaluate the integral over the fit range.
    Double_t integral = rds_withmass4mucuts.sumEntries();
    // cout << "Fitting over range:\n" << sel_mass4mu << endl;
    // cout << "Using cuts: " << all_cuts << endl;
    integral_vec.push_back(integral);
    TString fitrange_txt = Form("#bf{#splitline{Fit range:}{%.1f < m_{4#mu} < %.1f GeV}}", m4mu_min, m4mu_max);
    TString integ = Form("#bf{Integral = %.0f}", integral);
    TPaveText* pave = new TPaveText(0.68, 0.65, 0.905, 0.80, "NDC");  //  # NDC = normalized coord.
    pave->SetFillColor(0);  //
    pave->SetFillStyle(1001);  //  # Solid fill.
    pave->SetBorderSize(1);  //) # Use 0 for no border.
    pave->SetTextAlign(11);  // # 11 is against left side, 22 is centered vert and horiz.
    pave->SetTextSize(size_text);  //
    pave->AddText(fitrange_txt);  //  # (x, y, "text")
    pave->AddText(integ);
    pave->Draw("same");
    
    // TLatex* tex_integ = new TLatex(0.20, 0.45, integ);
    // tex_integ->SetNDC();
    // tex_integ->SetTextColor(2);
    // tex_integ->SetTextSize(size_text);
    // tex_integ->Draw("same");

    
    // TLatex* tex_fitrange = new TLatex(0.20, 0.38, fitrange_txt);
    // tex_fitrange->SetNDC();
    // tex_fitrange->SetTextColor(2);
    // tex_fitrange->SetTextSize(size_text);
    // tex_fitrange->Draw("same");

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