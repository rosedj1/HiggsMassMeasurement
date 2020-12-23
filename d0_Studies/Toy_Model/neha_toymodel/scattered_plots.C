 #include "TRandom.h"
#include "TAxis.h"
// for 1000 pts 
   void scattered_plots() {

   gROOT->SetBatch(kTRUE);
   Double_t sigma_st = .0001;
   Double_t sigma1 = .01; 
   Double_t a;
   TRandom3 *random_number = new TRandom3();


// defining pT values

   Double_t par[7] = {.00114,.00057,.000285,.0001425,.000114,.000076,.000057}; //correct
  // Double_t par[1];
 //  par[0] = {.0001425};
   Double_t pT[7] = {5,10,20,40,50,75,100};
 //  Double_t pT[1];
 //  pT[0] = {40};
/*
   Double_t par[7] = {.000568,.000142,.00114,.000114,.000285,.000057,.000076};
   Double_t pT[7] = {10,40,5,50,20,100,75};*/

// define graphs for scattered plot AND for gain and resolution with pT

   TGraph *gr1[7], *gr2[7], *gr_corrected[7], *gr_corrected_a[7];
   //TGraph *gain_plot;
   TGraph *sigma_pt, *sigma_pt_fit;
  // TGraph *gain_plot_a;
   TGraph *sigma_a, *sigma_a_fit;
   TGraphErrors *gain_plot, *gain_plot_a;
 
// defining canvas for plots
   TCanvas *pT_plots[7], *pT_plots_fit[7];
   TF1 *line_fit[7];
   TCanvas *a_plots[7], *a_plots_fit[7], *c_rand_dist[7];
   TF1 *line_fit_a[7];
   TH1F* h_rand[7];
   Int_t number_loop;

// To save pdf for scattered plots
   TString scattered_plot[7]={"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_5_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_10_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_20_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_40_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_50_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_75_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_100_second.pdf"};
  
   TString scattered_plot_fit[7]={"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_5_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_10_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_20_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_40_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_50_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_75_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/pt_100_second_fit.pdf"};

   Double_t slope_fit[7],intercept_fit[7];

   // Double_t pt_range[7] = {.007,.013,.028,0.065,0.065,.08,.12};
   // Double_t number[7] = {100,100,100,100,100,100,100};
   // Double_t pt_range_before_below[7] = {.01,.025,.038,0.08,0.095,.15,.18};
   // Double_t pt_range_before_above[7] = {.01,.025,.038,0.08,0.09,.15,.18};

  TString scattered_plot_a[7]={"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_5_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_10_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_20_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_40_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_50_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_75_second.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_100_second.pdf"};
  
   TString scattered_plot_a_fit[7]={"/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_5_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_10_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_20_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_40_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_50_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_75_second_fit.pdf",
   "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/a_100_second_fit.pdf"};
  
   TString outpath_pdf = "/blue/avery/rosedj1/HiggsMassMeasurement/d0_Studies/Toy_Model/output/plots/";
   TString fullpath_canv_rand_dist[7]={
      outpath_pdf + "randnum_pT5.pdf",
      outpath_pdf + "randnum_pT10.pdf",
      outpath_pdf + "randnum_pT20.pdf",
      outpath_pdf + "randnum_pT40.pdf",
      outpath_pdf + "randnum_pT50.pdf",
      outpath_pdf + "randnum_pT75.pdf",
      outpath_pdf + "randnum_pT100.pdf"
      };

   Double_t slope_fit_a[7],intercept_fit_a[7];
   TCanvas *parabolic_graphs[1002];
//   TString plots_parabola[5]={"first.pdf","second.pdf","third.pdf","fourth.pdf","fifth.pdf"};

  // for pT
/*
    Double_t pt_range[7] = {.013,.065,.007,0.065,0.028,.12,.08};
    Double_t number[7] = {100,100,100,100,100,100,100};
    Double_t pt_range_before_below[7] = {.025,0.08,.01,.095,.038,.18,.15};
    Double_t pt_range_before_above[7] = {.025,.08,.01,0.09,0.038,.18,.15};
*/
// defining pt ranges for projection plots
  // Double_t pt_range[7]={.015,.025,.04,.07,.100,.130,.160};
  // Double_t number[7]={150,250,400,700,1000,1300,1600};
   //Double_t name,name2;

// defining histograms for projection plots
/*   TH1F *h[7],*pt_projection[7],*h_fit[7],*pt_projection_fit[7];
   for(Int_t t=0;t<7;t++){
    TString name = TString::Format("h[%d]",t);
    TString name2 = TString::Format("pt_projection[%d]",t);
    TString name_fit = TString::Format("h_fit[%d]",t);
    TString name2_fit = TString::Format("pt_projection_fit[%d]",t);

    h[t] = new TH1F(name,"projection before correction",150,-.025,.025);
    h[t]->GetXaxis()->SetTitle("d0 (cm)");
    h[t]-> GetXaxis() -> CenterTitle();
    h[t]->GetYaxis()->SetTitle("no. of events");


    pt_projection[t] = new TH1F(name2,"projection before correction",number[t],-pt_range_before_below[t],pt_range_before_above[t]);
    pt_projection[t]->GetXaxis()->SetTitle("#frac{#Delta p_{T}}{p_{T}}");
    pt_projection[t]->GetXaxis()->CenterTitle();
    pt_projection[t]->GetYaxis()->SetTitle("no. of events");
   // pt_projection[t]->GetXaxis()->SetRangeUser(-.18, .18);


    h_fit[t] = new TH1F(name_fit,"projection after correction",150,-.025,.025);
    h_fit[t]->GetXaxis()->SetTitle("d0 (cm)");
    h_fit[t]-> GetXaxis() -> CenterTitle();
    h_fit[t]->GetYaxis()->SetTitle("no. of events");

    pt_projection_fit[t] = new TH1F(name2_fit,"projection after correction",number[t],-pt_range[t],pt_range[t]);
    pt_projection_fit[t]->GetXaxis()->SetTitle("#frac{#Delta p_{T}}{p_{T}}");
    pt_projection_fit[t]-> GetXaxis() -> CenterTitle();
    pt_projection_fit[t]->GetYaxis()->SetTitle("no. of events");

   }

   //for a 

      TH1F *h_a[7],*a_projection[7],*h_a_fit[7],*a_projection_fit[7];
   for(Int_t t=0;t<7;t++){
    TString name_a = TString::Format("h_a[%d]",t);
    TString name2_a = TString::Format("a_projection[%d]",t);
    TString name_a_fit = TString::Format("h_a_fit[%d]",t);
    TString name2_a_fit = TString::Format("a_projection_fit[%d]",t);

    h_a[t] = new TH1F(name_a,"projection before correction",150,-.025,.025);
    h_a[t]->GetXaxis()->SetTitle("d0 (cm)");
    h_a[t]-> GetXaxis() -> CenterTitle();
    h_a[t]->GetYaxis()->SetTitle("no. of events");


    a_projection[t] = new TH1F(name2_a,"projection before correction",number[t],-pt_range_before_below[t],pt_range_before_above[t]);
    a_projection[t]->GetXaxis()->SetTitle("#frac{#Delta a}{a}");
    a_projection[t]->GetXaxis()->CenterTitle();
    a_projection[t]->GetYaxis()->SetTitle("no. of events");
   // pt_projection[t]->GetXaxis()->SetRangeUser(-.18, .18);


    h_a_fit[t] = new TH1F(name_a_fit,"projection after correction",150,-.025,.025);
    h_a_fit[t]->GetXaxis()->SetTitle("d0 (cm)");
    h_a_fit[t]-> GetXaxis() -> CenterTitle();
    h_a_fit[t]->GetYaxis()->SetTitle("no. of events");

    a_projection_fit[t] = new TH1F(name2_a_fit,"projection after correction",number[t],-pt_range[t],pt_range[t]);
    a_projection_fit[t]->GetXaxis()->SetTitle("#frac{#Delta a}{a}");
    a_projection_fit[t]-> GetXaxis() -> CenterTitle();
    a_projection_fit[t]->GetYaxis()->SetTitle("no. of events");

   }*/
//--- MAIN ---//
// first loop for 7 different pT values
  for(Int_t pt_loop=0; pt_loop<7; pt_loop++){
      a = par[pt_loop];
      Double_t slopevalue[1002], slopeerr[1002], interceptvalue[1002],intercepterr[1002],d0value[1002], d0err[1002],diff_slope[1002];
      Double_t pT_graph[1002],pT_slope[1002],pT_corrected[1002],new_dpT[1002],new_da[1002];
      TString h_title = Form("Dist. of values from pdf: G(0, 0.01) for pT=%f", pT[pt_loop]);
      TString h_x_label = "arbitrary independent var";
      TString h_y_label = "frequency";
      TString all_titles = Form("%s;%s;%s", h_title.Data(), h_x_label.Data(), h_y_label.Data());
      h_rand[pt_loop] = new TH1F(Form("h_rand_%f", pT[pt_loop]), all_titles, 500, -0.05, 0.05);
      cout << "hist number " << pt_loop << " successfully made." << endl;
      number_loop =0;
   // for tossing toys for a single pT
      for(Int_t k =0; k<1002;k++){
         
      
         Double_t y[14];
         Double_t x[14] = {2.9,6.8,10.9,16.0,25.5,33.9,41.9,49.8,60.8,69.2,78.0,86.8,96.5,108.0}; //correct.
         for (Int_t i=0;i<14;i++) {
            // getting random points      
            Double_t rand = random_number->Gaus(0,sigma1);
            y[i] = a*x[i]*x[i] + rand;
            h_rand[pt_loop]->Fill(rand);
         }
         TString parabolic_name = TString::Format("parabolic_graphs_pT%f_%d",pT[pt_loop],k);
      // fitting random points in a graph and obtaining fit parameters
         parabolic_graphs[k] = new TCanvas(parabolic_name,"fit"); //,200,10,700,500);
         parabolic_graphs[k]->cd();
      
         TGraph *gr = new TGraph(14,x,y);
      // TGraph *gr = new TGraph(n,x,y);
         gr->SetLineColor(2);
         gr->SetLineWidth(4);
         gr->SetMarkerColor(4);
         gr->SetMarkerStyle(21);
         gr->SetTitle("parabolic_fit");
         gr->GetXaxis()->SetTitle("x");
         gr->GetYaxis()->SetTitle("y");
         gr->Draw("AP*");

         TF1 *line0 = new TF1("line0","[2]*x*x+[1]*x+[0]",0,150);
         line0->SetParameter(2,a);
         line0->SetParameter(1,0.0);
         line0->SetParameter(0,0.0);
         gr->Fit("line0");
         gr->Paint("|>");

         // TString parabola_name = plots_parabola[k];

         gStyle->SetOptFit();
         // parabolic_graphs[k]->SaveAs(parabola_name);

         slopevalue[k] = line0->GetParameter(2); // a_fit
         slopeerr[k] = line0->GetParError(2);    // a_fit_err
         interceptvalue[k] = line0->GetParameter(1); // NOT USED.
         intercepterr[k] = line0->GetParError(1);    // NOT USED.
         d0value[k] = line0->GetParameter(0);
         d0err[k] = line0->GetParError(0);
         diff_slope[k] = (slopevalue[k] - a)/a;
      
         pT_graph[k] = (0.0057/slopevalue[k]);
         pT_slope[k] = -((pT_graph[k] - pT[pt_loop])/pT[pt_loop]);
         number_loop++;
         
      }
      cout<<"number times loop  "<<number_loop<<endl;
      cout<<"pT  "<<"a  "<<"a from graph   "<<" pt from graph "<<"#Delta a/a   "<<"#Delta pT/pT  "<<" loop number  "<<endl; 
      Int_t i;
   // for(i = 0;i<1002;i++){
   //  cout<<pT[pt_loop]<<"  "<<par[pt_loop]<<"  "<<slopevalue[i]<<"  "<<pT_graph[i]<<"  "<<diff_slope[i]<<"  "<<pT_slope[i]<<"  "<<pt_loop<<endl;
   // }
      cout<<"number of values printed"<<i<<endl;
      /*
      for(Int_t fl=0; fl<1002; fl++){
      h[pt_loop]->Fill(d0value[fl]);
      pt_projection[pt_loop]->Fill(pT_slope[fl]);
      h_a[pt_loop]->Fill(d0value[fl]);
      a_projection[pt_loop]->Fill(diff_slope[fl]);

      }*/
   // scattered plots gr1 and gr2 
      gr1[pt_loop]= new TGraph(1002,d0value,diff_slope); // diff_slope = (a_fit - a)/a
      gr2[pt_loop]= new TGraph(1002,d0value,pT_slope);

      // for pT

      TString name_plots = TString::Format("pT_plots[%d]",pt_loop);
      TString name_plots_fit = TString::Format("pT_plots_fit[%d]",pt_loop);

      TString relation_plot = scattered_plot[pt_loop];
      TString relation_plot_fit = scattered_plot_fit[pt_loop];

      pT_plots[pt_loop] = new TCanvas(name_plots,"scattered plot pt",1000,1000);

      pT_plots[pt_loop]->cd();
      gr2[pt_loop]->SetMarkerColor(kBlue);
      gr2[pt_loop]->SetMarkerSize(.2);
      gr2[pt_loop]->SetMarkerStyle(20);
      gr2[pt_loop]->SetLineColor(kBlue);
   // gr2[0]->SetName("g1");
      gr2[pt_loop]->GetXaxis()->SetRangeUser(-.02, .02);
      gr2[pt_loop]->GetYaxis()->SetRangeUser(-.15, .15);

      gr2[pt_loop]->GetXaxis()->SetTitle("d0 (cm)");
      gr2[pt_loop]->GetYaxis()->SetTitle(" #frac{#Delta p_{T}}{p_{T}}");
      //gr->GetYaxis()->SetRangeUser(0,.3);
      gr2[pt_loop]->SetTitle(" #frac{#Delta p_{T}}{p_{T}} vs d0");
      gr2[pt_loop]->Draw("AP");
      pT_plots[pt_loop]->SaveAs(relation_plot);
   
   pT_plots_fit[pt_loop]  = new TCanvas(name_plots_fit);

   pT_plots_fit[pt_loop]->cd();

   TString fit_line_name = TString::Format("line_fit[%d]",pt_loop);
   //TString line_name2 = TString::Format("line_fit[%d]",pt_loop);

   line_fit[pt_loop] = new TF1(fit_line_name,"[1]*x+[0]",-10,10);
   line_fit[pt_loop]->SetParameter(0,0);
   line_fit[pt_loop]->SetParameter(1,1.0);
   line_fit[pt_loop]->SetParName(0,"intercept");  
   line_fit[pt_loop]->SetParName(1,"slope");  
   gr2[pt_loop]->Fit(fit_line_name);
   gr2[pt_loop]->Draw("AP");
   pT_plots_fit[pt_loop]->SaveAs(relation_plot_fit);

      // a plots

      TString name_plots_a = TString::Format("a_plots[%d]",pt_loop);
      TString name_plots_a_fit = TString::Format("a_plots_fit[%d]",pt_loop);
      TString name_c_rand_dist = TString::Format("c_rand_dist[%d]",pt_loop);

      TString relation_plot_a = scattered_plot_a[pt_loop];
      TString relation_plot_a_fit = scattered_plot_a_fit[pt_loop];
      TString this_canv_fullpath = fullpath_canv_rand_dist[pt_loop];

      a_plots[pt_loop] = new TCanvas(name_plots_a,"scattered plot with a ",1000,1000);

      a_plots[pt_loop]->cd();
      gr1[pt_loop]->SetMarkerColor(kBlue);
      gr1[pt_loop]->SetMarkerSize(.2);
      gr1[pt_loop]->SetMarkerStyle(20);
      gr1[pt_loop]->SetLineColor(kBlue);
   // gr2[0]->SetName("g1");
      gr1[pt_loop]->GetXaxis()->SetRangeUser(-.02, .02);
      gr1[pt_loop]->GetYaxis()->SetRangeUser(-.15, .15);

      gr1[pt_loop]->GetXaxis()->SetTitle("d0 (cm)");
      gr1[pt_loop]->GetYaxis()->SetTitle(" #frac{#Delta a}{a}");
      //gr->GetYaxis()->SetRangeUser(0,.3);
      gr1[pt_loop]->SetTitle(" #frac{#Delta a}{a} vs d0");
      gr1[pt_loop]->Draw("AP");
      a_plots[pt_loop]->SaveAs(relation_plot_a);
   
   a_plots_fit[pt_loop]  = new TCanvas(name_plots_a_fit);

   a_plots_fit[pt_loop]->cd();

   TString fit_line_name_a = TString::Format("line_fit_a[%d]",pt_loop);
   //TString line_name2 = TString::Format("line_fit[%d]",pt_loop);

   line_fit_a[pt_loop] = new TF1(fit_line_name_a,"[1]*x+[0]",-10,150);
   line_fit_a[pt_loop]->SetParameter(0,0);
   line_fit_a[pt_loop]->SetParameter(1,1.0);
   line_fit_a[pt_loop]->SetParName(0,"intercept");  
   line_fit_a[pt_loop]->SetParName(1,"slope");  
   gr1[pt_loop]->Fit(fit_line_name_a);
   gr1[pt_loop]->Draw("AP");
   a_plots_fit[pt_loop]->SaveAs(relation_plot_a_fit);

   c_rand_dist[pt_loop] = new TCanvas(name_c_rand_dist);
   h_rand[pt_loop]->Draw("hist");
   gStyle->SetOptStat("iouRMe");
   c_rand_dist[pt_loop]->SaveAs(this_canv_fullpath);
   /*
   // correction applied here
   slope_fit[pt_loop] = line_fit[pt_loop]->GetParameter(1);
   intercept_fit[pt_loop] = line_fit[pt_loop]->GetParameter(0);

   slope_fit_a[pt_loop] = line_fit_a[pt_loop]->GetParameter(1);
   intercept_fit_a[pt_loop] = line_fit_a[pt_loop]->GetParameter(0);

      for(Int_t i=0; i<1002;i++){
      new_dpT[i] = pT_slope[i] - (intercept_fit[pt_loop] + slope_fit[pt_loop] * d0value[i]);
      }

      for(Int_t i=0; i<1002;i++){
      new_da[i] = diff_slope[i] - (intercept_fit_a[pt_loop] + slope_fit_a[pt_loop] * d0value[i]);
      }

      for(Int_t fl=0; fl<1002; fl++){
      h_fit[pt_loop]->Fill(d0value[fl]);
      pt_projection_fit[pt_loop]->Fill(new_dpT[fl]);
      h_a_fit[pt_loop]->Fill(d0value[fl]);
      a_projection_fit[pt_loop]->Fill(new_da[fl]);
   }

   // corrected scattered plot
   gr_corrected[pt_loop] = new TGraph(1002,d0value,new_dpT);
   gr_corrected_a[pt_loop] = new TGraph(1002,d0value,new_da);*/
 }
// for loop ending here

// projection after correction
/* Int_t i;
 TCanvas *fit_final_pt = new TCanvas("fit_final_pt","projection after correction",1200,1200);
 fit_final_pt->Divide(3,3);
 for(Int_t l=1;l<8;l++){
  i = l-1;
  fit_final_pt->cd(l);
  gr_corrected[i]->SetMarkerColor(kBlue);
   gr_corrected[i]->SetMarkerSize(.5);
   gr_corrected[i]->SetMarkerStyle(20);
   gr_corrected[i]->SetLineColor(kBlue);
  // gr2[0]->SetName("g1");
   gr_corrected[i]->GetXaxis()->SetRangeUser(-.02, .02);
   gr_corrected[i]->GetYaxis()->SetRangeUser(-.15, .15);
//   gr2[0]->GetYaxis()->SetRangeUser(-.01, .02);

   gr_corrected[i]->GetXaxis()->SetTitle("d0 (cm)");
   gr_corrected[i]->GetYaxis()->SetTitle(" #frac{#Delta p_{T}}{p_{T}}");
   //gr->GetYaxis()->SetRangeUser(0,.3);
   gr_corrected[i]->SetTitle(" #frac{#Delta p_{T}}{p_{T}} vs d0");
   gr_corrected[i]->Draw("AP");

   }
  fit_final_pt->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/after_refitting_pt.pdf");


 TCanvas *fit_final_a = new TCanvas("fit_final_a","projection after correction",1200,1200);
 fit_final_a->Divide(3,3);
 for(Int_t l=1;l<8;l++){
  i = l-1;
  fit_final_a->cd(l);
  gr_corrected_a[i]->SetMarkerColor(kBlue);
   gr_corrected_a[i]->SetMarkerSize(.5);
   gr_corrected_a[i]->SetMarkerStyle(20);
   gr_corrected_a[i]->SetLineColor(kBlue);
  // gr2[0]->SetName("g1");
   gr_corrected_a[i]->GetXaxis()->SetRangeUser(-.02, .02);
   gr_corrected_a[i]->GetYaxis()->SetRangeUser(-.15, .15);
//   gr2[0]->GetYaxis()->SetRangeUser(-.01, .02);

   gr_corrected_a[i]->GetXaxis()->SetTitle("d0 (cm)");
   gr_corrected_a[i]->GetYaxis()->SetTitle(" #frac{#Delta a}{a}");
   //gr->GetYaxis()->SetRangeUser(0,.3);
   gr_corrected_a[i]->SetTitle(" #frac{#Delta a}{a} vs d0");
   gr_corrected_a[i]->Draw("AP");

   }
  fit_final_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/after_refitting_a.pdf");


TCanvas *projection_d_pt = new TCanvas("projection_d_pt","projection d",1200,1200);
TCanvas *projection_p_pt = new TCanvas("projection_p_pt","projection p",1200,1200);

projection_d_pt->Divide(3,3);
projection_p_pt->Divide(3,3);

TF1 *f1[7],*f2[7];
Double_t sigma_before[7], sigma_before_err[7], sigma_after[7],sigma_after_err[7],gain[7],gain_err[7];
Int_t l;


for(Int_t fl=1;fl<8;fl++){
   l = fl-1;
     TString gaus1 = TString::Format("f1[%d]",l);
  f1[l] = new TF1(gaus1,"gaus",-pt_range[l],pt_range[l]);
  projection_d_pt->cd(fl);
  h[l]->Fit(gaus1);
  gStyle->SetOptStat(1101);
  h[l]->Draw();

  TString gaus2 = TString::Format("f2[%d]",l);
  f2[l] = new TF1(gaus2,"gaus",-pt_range[l],pt_range[l]);
  projection_p_pt->cd(fl);
  pt_projection[l]->Fit(gaus2);
  sigma_before[l] =  f2[l]->GetParameter(2);
  sigma_before_err[l] = f2[l]->GetParError(2);

  gStyle->SetOptStat(1101);
  pt_projection[l]->Draw();

}


TCanvas *projection_d_a = new TCanvas("projection_d_a","projection d",1200,1200);
TCanvas *projection_a_a = new TCanvas("projection_a_a","projection a",1200,1200);

projection_d_a->Divide(3,3);
projection_a_a->Divide(3,3);

TF1 *f1_a[7],*f2_a[7];
Double_t sigma_before_a[7], sigma_after_a[7],gain_a[7], sigma_before_a_err[7], sigma_after_a_err[7],gain_a_err[7];


for(Int_t fl=1;fl<8;fl++){
   l = fl-1;
     TString gaus1 = TString::Format("f1_a[%d]",l);
  f1_a[l] = new TF1(gaus1,"gaus",-pt_range[l],pt_range[l]);
  projection_d_a->cd(fl);
  h_a[l]->Fit(gaus1);
  gStyle->SetOptStat(1101);
  h_a[l]->Draw();

  TString gaus2 = TString::Format("f2_a[%d]",l);
  f2_a[l] = new TF1(gaus2,"gaus",-pt_range[l],pt_range[l]);
  projection_a_a->cd(fl);
  a_projection[l]->Fit(gaus2);
  sigma_before_a[l] =  f2_a[l]->GetParameter(2);
  sigma_before_a_err[l] = f2_a[l]->GetParError(2);


  gStyle->SetOptStat(1101);
  a_projection[l]->Draw();

}


TCanvas *projection_d_fit = new TCanvas("projection_d_fit","projection d after correction",1200,1200);
TCanvas *projection_p_fit = new TCanvas("projection_p_fit","projection p after correction",1200,1200);

projection_d_fit->Divide(3,3);
projection_p_fit->Divide(3,3);
TF1 *f1_fit[7],*f2_fit[7];


for(Int_t fl=1;fl<8;fl++){
   l = fl-1;
  
  TString gaus1_fit = TString::Format("f1_fit[%d]",l);
  f1_fit[l] = new TF1(gaus1_fit,"gaus",-pt_range[l],pt_range[l]);

  projection_d_fit->cd(fl);
  h_fit[l]->Fit(gaus1_fit);
  gStyle->SetOptStat(1101);
  h_fit[l]->Draw();

  TString gaus2_fit = TString::Format("f2_fit[%d]",l);
  f2_fit[l] = new TF1(gaus2_fit,"gaus",-pt_range[l],pt_range[l]);
  projection_p_fit->cd(fl);
  pt_projection_fit[l]->Fit(gaus2_fit);
  gStyle->SetOptStat(1101);
  sigma_after[l] =  f2_fit[l]->GetParameter(2);
  sigma_after_err[l] = f2_fit[l]->GetParError(2);
 
  pt_projection_fit[l]->Draw();

}

projection_d_fit->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_after_fit_d_pt.pdf");
projection_p_fit->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_after_fit_p_pt.pdf");
projection_d_pt->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_before_fit_d_pt.pdf");
projection_p_pt->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_before_fit_p_pt.pdf");



TCanvas *projection_d_fit_a = new TCanvas("projection_d_fit_a","projection d after correction",1200,1200);
TCanvas *projection_a_fit_a = new TCanvas("projection_a_fit_a","projection a after correction",1200,1200);

projection_d_fit_a->Divide(3,3);
projection_a_fit_a->Divide(3,3);
TF1 *f1_fit_a[7],*f2_fit_a[7];


for(Int_t fl=1;fl<8;fl++){
   l = fl-1;
  
  TString gaus1_fit = TString::Format("f1_fit_a[%d]",l);
  f1_fit_a[l] = new TF1(gaus1_fit,"gaus",-pt_range[l],pt_range[l]);

  projection_d_fit_a->cd(fl);
  h_a_fit[l]->Fit(gaus1_fit);
  gStyle->SetOptStat(1101);
  h_a_fit[l]->Draw();

  TString gaus2_fit = TString::Format("f2_fit_a[%d]",l);
  f2_fit_a[l] = new TF1(gaus2_fit,"gaus",-pt_range[l],pt_range[l]);
  projection_a_fit_a->cd(fl);
  a_projection_fit[l]->Fit(gaus2_fit);
  gStyle->SetOptStat(1101);
  sigma_after_a[l] =  f2_fit_a[l]->GetParameter(2);
  sigma_after_a_err[l] = f2_fit_a[l]->GetParError(2);
  a_projection_fit[l]->Draw();

}

projection_d_fit_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_after_fit_d_a.pdf");
projection_a_fit_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_after_fit_a_a.pdf");
projection_d_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_before_fit_d_a.pdf");
projection_a_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/projection_before_fit_a_a.pdf");

for(Int_t i = 0; i<7;i++){
  gain[i] = ((sigma_before[i] - sigma_after[i])/sigma_before[i]) *100 ; 
  gain_err[i] = ((sqrt(pow(sigma_before[i],2) * pow(sigma_after_err[i],2) + pow(sigma_after[i],2) * pow(sigma_before_err[i],2)))/pow(sigma_before[i],2));

 // gain_err[i] = ((sqrt(sigma_before[i]^2 * sigma_after_err[i]^2 + sigma_after[i]^2 * sigma_before_err[i]^2))/sigma_before[i]^2);
  //gain[i] = ((sigma_before[i] - sigma_after[i])/sigma_before[i]) *100 ; 

  cout<<sigma_before[i]<<"   "<<sigma_after[i]<<"  "<<gain[i]<<endl;

    gain_a[i] = ((sigma_before_a[i] - sigma_after_a[i])/sigma_before_a[i]) *100 ; 
      gain_a_err[i] = ((sqrt(pow(sigma_before_a[i],2) * pow(sigma_after_a_err[i],2) + pow(sigma_after_a[i],2) * pow(sigma_before_a_err[i],2)))/pow(sigma_before_a[i],2));

  cout<<sigma_before_a[i]<<"   "<<sigma_after_a[i]<<"  "<<gain_a[i]<<endl;
}

Double_t null[7] = {0,0,0,0,0,0,0};
gain_plot = new TGraphErrors(7,pT,gain,null,gain_err);
sigma_pt = new TGraph(7,pT,sigma_before);
sigma_pt_fit = new TGraph(7,pT,sigma_after);

gain_plot_a = new TGraphErrors(7,par,gain_a,null,gain_a_err);
sigma_a = new TGraph(7,par,sigma_before_a);
sigma_a_fit = new TGraph(7,par,sigma_after_a);

TCanvas *fi = new TCanvas("fi","gain with pt",500,500);
fi->cd();
 gain_plot->SetMarkerColor(kBlue);
 gain_plot->SetMarkerSize(.8);
 gain_plot->SetMarkerStyle(20);
 gain_plot->SetLineColor(kBlue);
 gain_plot->GetXaxis()->SetTitle("pT (GeV)");
 gain_plot->GetYaxis()->SetTitle(" gain %");
 gain_plot->SetTitle(" gain with pT");
 gain_plot->GetYaxis()->SetRangeUser(0,35);
 gain_plot->Draw("AP*");

 fi->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/gain_pT.pdf");

  TF1 *sigma_st_line = new TF1("sigma_st_line","[0]*x+[1]",0,150);
  sigma_st_line->SetParameter(0,1);
  sigma_st_line->SetParameter(1,0.0);
  sigma_st_line->SetParName(1,"intercept");  
  sigma_st_line->SetParName(0,"slope");  


TCanvas *si1 = new TCanvas("si1","sigma dpT/pT before correction",800,800);
si1->cd();
sigma_pt->SetMarkerColor(kBlue);
sigma_pt->SetMarkerSize(.8);
sigma_pt->SetMarkerStyle(20);
sigma_pt->SetLineColor(kBlue);
sigma_pt->GetXaxis()->SetTitle("pT (GeV) ");
sigma_pt->GetYaxis()->SetTitle(" #sigma(dpT/pT)");
sigma_pt->GetYaxis()->SetRangeUser(0.0,.05);

sigma_pt->SetTitle(" resolution with pT (before correction)");
sigma_pt->Fit(sigma_st_line);
sigma_pt->Draw("AP");

si1->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/sigma_dpT.pdf");

  TF1 *sigma_after_st_line = new TF1("sigma_after_st_line","[0]*x+[1]",0,150);
  sigma_after_st_line->SetParameter(0,1);
  sigma_after_st_line->SetParameter(1,0.0);
  sigma_after_st_line->SetParName(1,"intercept");  
  sigma_after_st_line->SetParName(0,"slope"); 

TCanvas *si2= new TCanvas("si2","sigma dpT/pT after correction",800,800);
si2->cd();
sigma_pt_fit->SetMarkerColor(kBlue);
sigma_pt_fit->SetMarkerSize(.8);
sigma_pt_fit->SetMarkerStyle(20);
sigma_pt_fit->SetLineColor(kBlue);
sigma_pt_fit->GetXaxis()->SetTitle("pT (GeV)");
sigma_pt_fit->GetYaxis()->SetTitle(" #sigma(dpT/pT)");
sigma_pt_fit->GetYaxis()->SetRangeUser(0.0,.05);

sigma_pt_fit->SetTitle(" resolution with a (after correction)");
sigma_pt_fit->Fit(sigma_after_st_line);
sigma_pt_fit->Draw("AP*");
si2->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/sigma_dpT_fit.pdf");


TCanvas *fi_a = new TCanvas("fi_a","gain with a",500,500);
fi_a->cd();
 gain_plot_a->SetMarkerColor(kBlue);
 gain_plot_a->SetMarkerSize(.8);
 gain_plot_a->SetMarkerStyle(20);
 gain_plot_a->SetLineColor(kBlue);
 gain_plot_a->GetXaxis()->SetTitle("a (cm^{-1})");
 gain_plot_a->GetYaxis()->SetTitle(" gain %");
 gain_plot_a->GetYaxis()->SetRangeUser(0,35);
 gain_plot_a->SetTitle(" gain with a");
 gain_plot_a->Draw("AP*");

 fi_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/gain_a.pdf");

  TF1 *sigma_st_line_a = new TF1("sigma_st_line_a","[0]*x+[1]",0,150);
  sigma_st_line_a->SetParameter(0,1);
  sigma_st_line_a->SetParameter(1,0.0);
  sigma_st_line_a->SetParName(1,"intercept");  
  sigma_st_line_a->SetParName(0,"slope");  


TCanvas *si1_a = new TCanvas("si1_a","sigma da/a before correction",800,800);
si1_a->cd();
sigma_a->SetMarkerColor(kBlue);
sigma_a->SetMarkerSize(.8);
sigma_a->SetMarkerStyle(20);
sigma_a->SetLineColor(kBlue);
sigma_a->GetXaxis()->SetTitle("a (cm^{-1})");
sigma_a->GetYaxis()->SetTitle(" #sigma(da/a)");
sigma_a->GetYaxis()->SetRangeUser(0.0,.05);

sigma_a->SetTitle(" resolution with a (before correction)");
sigma_a->Fit(sigma_st_line_a);
sigma_a->Draw("AP");
si1_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/sigma_da.pdf");

  TF1 *sigma_fit_st_line_a = new TF1("sigma_fit_st_line_a","[0]*x+[1]",0,10);
  sigma_fit_st_line_a->SetParameter(0,1);
  sigma_fit_st_line_a->SetParameter(1,0.0);
  sigma_fit_st_line_a->SetParName(1,"intercept");  
  sigma_fit_st_line_a->SetParName(0,"slope");  

TCanvas *si2_a = new TCanvas("si2_a","sigma da/a after correction",800,800);
si2_a->cd();
sigma_a_fit->SetMarkerColor(kBlue);
sigma_a_fit->SetMarkerSize(.8);
sigma_a_fit->SetMarkerStyle(20);
sigma_a_fit->SetLineColor(kBlue);
sigma_a_fit->GetXaxis()->SetTitle("a (cm^{-1})");
sigma_a_fit->GetYaxis()->SetTitle(" #sigma(da/a)");
sigma_a_fit->GetYaxis()->SetRangeUser(0.0,.05);

sigma_a_fit->SetTitle(" resolution with a (after correction)");
sigma_a_fit->Fit(sigma_fit_st_line_a);
sigma_a_fit->Draw("AP*");
si2_a->SaveAs("/home/neha/Documents/work/adhoc_toy_model_study_latest/plots_test/sigma_da_fit.pdf");*/

 }

 
 

  
