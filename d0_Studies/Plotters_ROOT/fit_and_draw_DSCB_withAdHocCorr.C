#include "RooFitResult.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TString.h"
#include "RooRealVar.h"
#include "RooMyPDF_DSCB.h"
#include "RooDataSet.h"
#include "TROOT.h"
#include "RooPlot.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TStyle.h"
#include "TAxis.h"
#include "RooArgList.h"

using namespace std;

RooFitResult* fit_and_draw_DSCB(
    TTree* tree,
    Float_t m4mu_min, Float_t m4mu_max,
    Float_t m4muErr_min, Float_t m4muErr_max,  // Passed in as a percent.
    vector<Double_t> &integral_vec,
    TString method, const char* year, 
    TCanvas* canv, TString outfile_path,
    Int_t n_bins=100,
    Int_t color_line=kRed,
    Bool_t show_after_corr=false,
    Double_t mean_before_corr=0.0,
    Double_t sigma_before_corr=0.0,
    Bool_t zoom=false) {
    // Double_t m4mu_min, Double_t m4mu_max,
    // Double_t m4muErr_min, Double_t m4muErr_max
    /** 
     * Do an unbinned DSCB fit on the mass4mu data stored in `tree`.
     * Apply certain selections in the RooDataSet.
     * Return the fit result.

     Parameters
     ----------
     n_bins : Number of bins to be placed across [105, 140] GeV
     */

    const char* fs = "ggH";
    // const char* year = "2017";
    // const char* corr_factor_sample = "DY+J/#psi";

    TString title_p1 = Form("MC %s %s Unbinned DSCB Fit", fs, year);
    TString title_p2 = Form(" (using %s p_{T} corr. factors)", method.Data());
    TString title = title_p1 + title_p2;
    // TString title_p2 = Form("%.2f < #deltam_{4#mu}/m_{4#mu} < %.2f%%", m4muErr_min, m4muErr_max);
    // TString title_p2 = Form(" (using p_{T} corr. factors derived from %s muons)", corr_factor_sample);
    // if (method == "AdHoc") {
    //     title_p2 = Form(" (using AdHoc p_{T} corr. factors)");
    // } else if (method == "GeoFit") {
    //     title_p2 = Form(" (using GeoFit p_{T} corr. factors)");
    // } else {
    //     throw std::invalid_argument("`method` must be either 'AdHoc' or 'GeoFit'.");
    // }

    Float_t gbl_x_min = 105.0;
    Float_t gbl_x_max = 140.0;
    Float_t size_text  = 0.017;  // 0.02
    Double_t tex_horiz = 0.20;
    Double_t tex_vert  = 0.82;
    Int_t color_marker = color_line + 2;

    if (!tree) {
        cout << "[ERROR] Tree not opened." << endl;
        // Throw exception here.
    }
    
    TString name_mu = "#mu";
    TString name_sig = "#sigma";
    TString name_aL = "#alpha_{L}";
    TString name_aR = "#alpha_{R}";
    TString name_expL = "n_{L}";
    TString name_expR = "n_{R}";

    TString corr_suffix = "^{corr.}";
    TString fit_param_title = "DSCB Fit Params:";

    if (show_after_corr) {
        tex_horiz = 0.70;
        name_mu   += corr_suffix;
        name_sig  += corr_suffix;
        name_aL   += corr_suffix;
        name_aR   += corr_suffix;
        name_expL += corr_suffix;
        name_expR += corr_suffix;
        fit_param_title = "DSCB Fit Params after corr:";
    }

    // Prepare the DSCB variables.
    /**
     * The right tail of the m(4mu) resonance decays very quickly.
     * This means that n_R (related to the exponential decay of the tail)
     * will typically be O(100), but capping it at 100 doesn't change
     * the other fit params much based on preliminary testing.
     * So leave it at 100 for faster fitting.
     */
    RooRealVar m4mu("m4mu", "m_{4#mu}", gbl_x_min, gbl_x_max, "GeV");
    RooRealVar m4mu_corr("m4mu_corr", "m_{4#mu}", gbl_x_min, gbl_x_max, "GeV");
    RooRealVar Mean(name_mu, "#mu", 125, 120, 130);
    RooRealVar Sigma(name_sig, "#sigma", 1, 0, 10);//sigma[decay]);
    RooRealVar AlphaL(name_aL, "#alpha_{L}", 1, 0, 30);//alphaL[decay]);
    RooRealVar AlphaR(name_aR, "#alpha_{R}", 1, 0, 30);//alphaR[decay]);
    RooRealVar ExpL(name_expL, "n_{L}", 1, 0, 50);//expL[decay]);
    RooRealVar ExpR(name_expR, "n_{R}", 1, 0, 100);//expR[decay]); (1,1,50)

    // RooMyPDF_DSCB DSCB;
    // if (show_after_corr) {
    //     DSCB = RooMyPDF_DSCB("DSCB", "DSCB", m4mu_corr, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // } else {
    //     DSCB = RooMyPDF_DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    // }
    RooMyPDF_DSCB DSCB = RooMyPDF_DSCB("DSCB", "DSCB", m4mu, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);
    RooMyPDF_DSCB DSCB_corr = RooMyPDF_DSCB("DSCB_corr", "DSCB", m4mu_corr, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);

    // Prepare the cuts.
    // TString cut1 = Form("%f < ", m4mu_min) + name_m4mu + " && "
    // TString cut1 = Form("%f > ", m4mu_min) + name_m4mu + " && "
    // TString sel_mass4mu = Form("%f < %s && %s < %f", m4mu_min, name_m4mu.Data(), name_m4mu.Data(), m4mu_max);
    cout << "***** Fitting a DSCB to m4mu with the following cuts: *****" << endl;
    TString all_cuts = Form("%f < m4mu && m4mu < %f", m4mu_min, m4mu_max);
    cout << "  * " << all_cuts << endl;

    // Make a RooDataSet to plot black points along full mass4mu range:
    RooDataSet rds("rds", "roodataset m4mu", tree, RooArgSet(m4mu), "");
    RooDataSet rds_corr("rds_corr", "roodataset m4mu_corr", tree, RooArgSet(m4mu_corr), "");
    RooDataSet rds_inmasswindow("rds_inmasswindow", "rds_inmasswindow m4mu", tree, RooArgSet(m4mu), all_cuts);
    cout << "rds.numEntries() =" << rds.numEntries() << endl;
    cout << "rds.sumEntries() =" << rds.sumEntries() << endl;
    cout << "rds_corr.numEntries() =" << rds_corr.numEntries() << endl;
    cout << "rds_corr.sumEntries() =" << rds_corr.sumEntries() << endl;
    cout << "rds_inmasswindow.numEntries() =" << rds_inmasswindow.numEntries() << endl;
    cout << "rds_inmasswindow.sumEntries() =" << rds_inmasswindow.sumEntries() << endl;

    // Do the DSCB fit.
    gStyle->SetOptStat("iouRMe");
    // RooFitResult* result = DSCB.fitTo(rds,
    //                              RooFit::Range(m4mu_min, m4mu_max),
    //                              RooFit::PrintLevel(-1),
    //                              RooFit::Save()
    //                              );
    RooFitResult* result;
    if (show_after_corr) {
        result = DSCB_corr.fitTo(rds_corr,
                                 RooFit::Range(m4mu_min, m4mu_max),
                                 RooFit::PrintLevel(-1),
                                 RooFit::Save()
                                 );
    } else {
        result = DSCB.fitTo(rds,
                                 RooFit::Range(m4mu_min, m4mu_max),
                                 RooFit::PrintLevel(-1),
                                 RooFit::Save()
                                 );
    }

    // TPad* ptop = new TPad("ptop", "pad main", 0.0, 0.25, 1.0, 1.0);
    // TPad* pbot = new TPad("pbot", "pad ratio", 0.0, 0.0, 1.0, 0.25);
    RooPlot* xframe;
    RooPlot* xframe_corr;
    if (zoom) {
        xframe = m4mu.frame(RooFit::Range(m4mu_min, m4mu_max), RooFit::Title(title));
        xframe_corr = m4mu_corr.frame(RooFit::Range(m4mu_min, m4mu_max), RooFit::Title(title));
    } else {
        xframe = m4mu.frame(RooFit::Range(gbl_x_min, gbl_x_max), RooFit::Title(title));
        xframe_corr = m4mu_corr.frame(RooFit::Range(gbl_x_min, gbl_x_max), RooFit::Title(title));
    }
    
    // Make binned version of data.
    // Should be no need to make a separate hist.
    // When the RooDataSet plots on the frame, it will look binned.
    TString hist_name = Form("h_m4mu_%fm4mu%f_%fm4muErr%f_%s", m4mu_min, m4mu_max, m4muErr_min, m4muErr_max, name_mu.Data());
    // h_m4mu->SetLineWidth(0);
    // xframe->addTH1(h_m4mu)
    TH1* h_m4mu = rds.createHistogram(hist_name, m4mu, RooFit::Binning(n_bins, gbl_x_min, gbl_x_max));
    Double_t max_val = h_m4mu->GetMaximum();
    xframe->SetMaximum(max_val * 1.1);
    // Double_t bin_width = h_m4mu->GetBinWidth(1);
    Double_t bin_width = (gbl_x_max - gbl_x_min) / n_bins;
    
    xframe->GetYaxis()->SetTitle(Form("Events / (%.2f GeV)", bin_width));
    // Int_t n_bins_in_fitrange = static_cast<Int_t>(round((m4mu_max - m4mu_min) / bin_width));
    // rds.plotOn(xframe, RooFit::MarkerColor(1), RooFit::Binning(n_bins_in_fitrange));
    rds.plotOn(xframe, RooFit::MarkerColor(color_marker), RooFit::Binning(n_bins));
    rds_corr.plotOn(xframe_corr, RooFit::MarkerColor(color_marker), RooFit::Binning(n_bins));
    // DSCB.paramOn(xframe, RooFit::Layout(tex_horiz, tex_horiz + 0.22, tex_vert - 0.02));
    // xframe->getAttText()->SetTextSize(size_text);
    // xframe->getAttText()->SetTextColor(color_line);
    // h_m4mu->SetTitle(title);
    // h_m4mu->Draw("same");  // For the stats box.
    if (show_after_corr) {
        DSCB_corr.plotOn(xframe_corr, RooFit::LineColor(color_line), RooFit::LineWidth(3), RooFit::Range(m4mu_min, m4mu_max));  //RooFit::Binning(n_bins), 
        xframe_corr->Draw("same");  // A plot was drawn before, so draw on it.
    } else {
        DSCB.plotOn(xframe, RooFit::LineColor(color_line), RooFit::LineWidth(3), RooFit::Range(m4mu_min, m4mu_max));  //RooFit::Binning(n_bins), 
        xframe->Draw();
    }

    // TLatex* tex = new TLatex(tex_horiz, tex_vert, fit_param_title);
    // tex->SetNDC();
    // tex->SetTextColor(color_line);
    // tex->SetTextSize(size_text);
    // tex->Draw("same");

    // Draw fit parameters.
    Double_t y_drop = 0.03;
    TLatex tex;
    tex.SetNDC();
    tex.SetTextColor(color_line);
    tex.SetTextSize(size_text);
    tex.DrawLatex(tex_horiz, tex_vert, fit_param_title);
    tex.DrawLatex(tex_horiz, tex_vert - 1.0*y_drop, Form("%s = %.3f #pm %.3f", name_aL.Data(), AlphaL.getVal(), AlphaL.getError()));
    tex.DrawLatex(tex_horiz, tex_vert - 2.0*y_drop, Form("%s = %.3f #pm %.3f", name_aR.Data(), AlphaR.getVal(), AlphaR.getError()));
    tex.DrawLatex(tex_horiz, tex_vert - 3.0*y_drop, Form("%s = %.3f #pm %.3f", name_mu.Data(), Mean.getVal(), Mean.getError()));
    tex.DrawLatex(tex_horiz, tex_vert - 4.0*y_drop, Form("%s = %.3f #pm %.3f", name_sig.Data(), Sigma.getVal(), Sigma.getError()));
    tex.DrawLatex(tex_horiz, tex_vert - 5.0*y_drop, Form("%s = %.3f #pm %.3f", name_expL.Data(), ExpL.getVal(), ExpL.getError()));
    tex.DrawLatex(tex_horiz, tex_vert - 6.0*y_drop, Form("%s = %.3f #pm %.3f", name_expR.Data(), ExpR.getVal(), ExpR.getError()));


    
    // Double_t y_max_fitparams = tex_vert - 0.02;
    // Double_t y_min_fitparams = y_max_fitparams - 0.02;
    // TPaveText* pave_fitparams = new TPaveText(tex_horiz, 0.185, tex_horiz + 0.22, 0.319, "NDC");  //  # NDC = normalized coord.
    // pave->SetFillColor(0);  // 0 = white
    // pave->SetFillStyle(1001);  //  # Solid fill.
    // pave->SetBorderSize(0);  //) # Use 0 for no border.
    // pave->SetTextAlign(11);  // # 11 is against left side, 22 is centered vert and horiz.
    // pave->SetTextSize(size_text);  //
    // pave->AddText(fitrange_txt);  //  # (x, y, "text")
    // pave->AddText(integ_txt);
    // // pave->SetTextColor(color_marker);
    // pave->Draw("same");


    Double_t bigger_txt_size = size_text + 0.01;
    if (show_after_corr) {
        // Draw the mass shift and the improvement in sigma.
        Double_t sigma_corr = Sigma.getVal();
        Double_t mass_shift = (Mean.getVal() - mean_before_corr) * 1000.0;  // Convert to MeV.
        Double_t improve = (sigma_corr - sigma_before_corr) / sigma_before_corr * 100.0;  // A percentage!
        TString imp_txt = Form("#bf{#left|#frac{#Delta#sigma}{#sigma}#right| = %.1f %%}", abs(improve));
        TString shift_txt = Form("#bf{#Delta#mu = %.0f MeV}", mass_shift);
        tex.SetNDC();
        tex.SetTextColor(1);
        tex.SetTextSize(bigger_txt_size);
        tex.DrawLatex(tex_horiz, 0.55, imp_txt);
        tex.DrawLatex(tex_horiz, 0.50, shift_txt);
        // TLatex* tex_imp = new TLatex(tex_horiz, 0.40, imp_txt);  // y-val = 0.25
        // tex_imp->SetNDC();
        // tex_imp->SetTextColor(2);
        // tex_imp->SetTextSize(size_text + 0.01);
        // tex_imp->Draw("same");
    }

    Double_t integral = rds_inmasswindow.sumEntries();
    // cout << "Fitting over range:\n" << sel_mass4mu << endl;
    // cout << "Using cuts: " << all_cuts << endl;
    integral_vec.push_back(integral);
    TString fitrange_txt = Form("#splitline{Fit range:}{%.1f < m_{4#mu} < %.1f GeV}", m4mu_min, m4mu_max);
    tex.SetTextColor(1);
    tex.SetTextSize(size_text);
    tex.DrawLatex(0.20, 0.55, fitrange_txt);

    // TString integ_txt = Form("#bf{Integral = %.0f}", integral);
    // TLatex* tex_integ = new TLatex(tex_horiz, 0.25, integ);
    // tex_integ->SetNDC();
    // tex_integ->SetTextColor(color_line);
    // tex_integ->SetTextSize(size_text);
    // tex_integ->Draw("same");
    // TPaveText* pave = new TPaveText(0.20, 0.185, 0.42, 0.319, "NDC");  //  # NDC = normalized coord.
    // pave->SetFillColor(0);  // 0 = white
    // pave->SetFillStyle(1001);  //  # Solid fill.
    // pave->SetBorderSize(0);  //) # Use 0 for no border.
    // pave->SetTextAlign(11);  // # 11 is against left side, 22 is centered vert and horiz.
    // pave->SetTextSize(size_text);  //
    // pave->AddText(fitrange_txt);  //  # (x, y, "text")
    // pave->AddText(integ_txt);
    // pave->SetTextColor(color_marker);
    // pave->Draw("same");


    
    
    // TLatex* tex_fitrange = new TLatex(0.20, 0.38, fitrange_txt);
    // tex_fitrange->SetNDC();
    // tex_fitrange->SetTextColor(2);
    // tex_fitrange->SetTextSize(size_text);
    // tex_fitrange->Draw("same");

// RooPlot* framePull = mass4mu.frame(RooFit::Title("")); // Ratio plot
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

    // canv->Print(outfile_path);// + ".pdf");

    return result;
}