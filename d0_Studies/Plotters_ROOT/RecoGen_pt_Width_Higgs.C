#include "TFile.h"
#include <iostream>
#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TROOT.h"
#include "TMath.h"
#include "TStyle.h"
#include "TColor.h"
#include "TLegend.h"
#include "TPad.h"
#include "TLine.h"
#include "TH1F.h"
#include "TPaveStats.h"
#include "TPad.h"
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include "Higgs_Mass_setup.h"


void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TString nome_canvas, TString save, TString x_axis, TLegend *legend, int fit, bool LogY = 0);
void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TString x_axis, TLegend *legend, int fit, bool LogY = 0);
void Draw(TH1F *h1, TH1F *h2, TH1F *h3, RooRealVar* rv_Mass_1, RooRealVar* rv_Mass_2, RooRealVar* rv_Mass_3, RooDataSet* Mass_1, RooDataSet* Mass_2, RooDataSet* Mass_3, TString nome_canvas, TString save, TString x_name, TLegend *legend, int fit, bool LogY);
void Draw_TH3F(TH3F* h1, TString nome_canvas, TString save, TLegend *legend);
void Draw_TH2F(TH2F* h1, TString nome_canvas, TString save, TString x_axis, TString y_axis);
void Draw_TH1(RooRealVar* rv_Mass_1, RooDataSet* Mass_1, TString title, TString save, float minRange, float maxRange);
void FillHisto(bool deltaR, TH3F* h1[5][3], int A, float eta1, float eta2, float phi1, float phi2, float pT1, float pT2, float mass, TString resonance);
void FillHisto(bool deltaR, TH3F* h1[5][3], int A, float eta1, float eta2, float phi1, float phi2, float pT1, float pT2, float val1, float val2, float mass, TString resonance);
void FillHisto(bool deltaR, TH1F* h1[3], int A, float eta1, float eta2, float phi1, float phi2, float pT1, float pT2, float val1, float val2, float mass, TString resonance);
float MCGen(float MC, float GEN);
float DeltaR1(float eta1, float phi1);
float DeltaR2(float eta1, float phi1);
float Max(float a, float b);
std::vector<float> FitMass(TH1F* h1, TString title, TString save_name);
std::vector<float> FitMass(RooRealVar* rv_Mass, RooDataSet* Data_Mass, TString title, TString save_name);
void Draw(bool dscb, bool sigma, int lepton_size, TH1F* FinalSigma_DSCB[3], TLegend* legend_width, TLegend* legend_comparison, TString save_name, int m, TLine* line);

bool muon;
TString save_nome;

using namespace std;

void RecoGen_pt_Width_Higgs(bool isDATA, int FS, TString year, bool isDR, TString channel, int Base) {
	
	TString file_dir;
	std::vector<int> *lep_id = 0; 
	std::vector<float> *GENlep_mass = 0;
	std::vector<float> *GENlep_pt = 0; 
	std::vector<float> *lep_pt = 0;
	std::vector<float> *lepFSR_pt = 0;
	std::vector<float> *vtxLep_pt = 0;
	std::vector<float> *vtxLep_BS_pt = 0;
	std::vector<float> *vtxLepFSR_pt = 0;
	std::vector<float> *vtxLepFSR_BS_pt = 0;
	std::vector<float> *GENlep_eta = 0;
	std::vector<float> *lep_eta = 0;
	std::vector<float> *lepFSR_eta = 0;
	std::vector<float> *vtxLep_eta = 0;
	std::vector<float> *vtxLep_BS_eta = 0;
	std::vector<float> *vtxLepFSR_eta = 0;
	std::vector<float> *vtxLepFSR_BS_eta = 0;
	std::vector<float> *GENlep_phi = 0;
	std::vector<float> *lep_phi = 0;
	std::vector<float> *lepFSR_phi = 0;
	std::vector<float> *vtxLep_phi = 0;
	std::vector<float> *vtxLep_BS_phi = 0;
	std::vector<float> *vtxLepFSR_phi = 0;
	std::vector<float> *vtxLepFSR_BS_phi = 0;
	std::vector<float> *lep_mass = 0;
	std::vector<float> *lepFSR_mass = 0;
	std::vector<float> *vtxLep_mass = 0;
	std::vector<float> *vtxLep_BS_mass = 0;
	std::vector<float> *vtxLepFSR_mass = 0;
	std::vector<float> *vtxLepFSR_BS_mass = 0;
	std::vector<int> *lep_numberOfValidPixelHits = 0;
	
	std::vector<float> *commonPV_x = 0;
	std::vector<float> *commonPV_y = 0;
	std::vector<float> *commonPV_z = 0;
	std::vector<float> *commonBS_x = 0;
	std::vector<float> *commonBS_y = 0;
	std::vector<float> *commonBS_z = 0;
	
	float PV_x;
	float PV_y;
	float PV_z;
	float BS_x;
	float BS_y;
	float BS_z;
	
	float pTL1;
	float etaL1;
	float phiL1;
	float mL1;
	float pTL2;
	float etaL2;
	float phiL2;
	float mL2;
	float pTL3;
	float etaL3;
	float phiL3;
	float mL3;
	float pTL4;
	float etaL4;
	float phiL4;
	float mL4;
	
	float GENMH;
	
	int finalState;
	bool passedFullSelection;
	float mass4l, mass4lREFIT, mass4mu, mass2e2mu, mass4l_vtx, mass4l_vtxFSR, mass4lREFIT_vtx, mass4l_vtx_BS, mass4l_vtxFSR_BS, mass4lREFIT_vtx_BS;
	float mass4lErr, mass4lErrREFIT, mass4lErr_vtx, mass4lErr_vtx_BS, mass4lErrREFIT_vtx, mass4lErrREFIT_vtx_BS;
	Float_t weight;	
	
	float newWeight;
// 	TH1F* LeptonDeltaRDistribution[4][2];
	TH1F* LeptonDeltaRDistribution[4][3];
	
	TString years[3] = {"2016", "2017", "2018"};
	
	gROOT->Reset();
	gROOT->SetBatch();
	
	pT_bins.clear();
	pT_bins.push_back(5);
	pT_bins.push_back(20);
	pT_bins.push_back(30);
	pT_bins.push_back(40);
	pT_bins.push_back(50);
	pT_bins.push_back(60);
	pT_bins.push_back(100);
	pT_bins.push_back(200);
	
	eta_bins.clear();
	eta_bins.push_back(0);
	eta_bins.push_back(0.9);
// 				eta_bins.push_back(1.2);
	eta_bins.push_back(1.4);
// 				eta_bins.push_back(1.8);
// 				eta_bins.push_back(2.1);
	eta_bins.push_back(2.4);

	eta_bins_name.clear();
	eta_bins_name.push_back("0");
	eta_bins_name.push_back("0p9");
// 				eta_bins_name.push_back("1p2");
	eta_bins_name.push_back("1p4");
// 				eta_bins_name.push_back("1p8");
// 				eta_bins_name.push_back("2p1");
	eta_bins_name.push_back("2p4");
	
	TH1F* PtDistribution[5][4];
	TH1F* SigmaPtOverPtDistribution[5][4];
	TH1F* DeltaPtDistribution[5][4];
	TH1F* EtaDistribution[5][4];
	TH1F* PhiDistribution[5][4];
	TH1F* MassDistribution[4];
	TH1F* REFITTEDMassDistribution[4];

	TH1F* W_MassDistribution[4];
	TH1F* W_REFITTEDMassDistribution[4];

	TH2F* SigmaPt_vsSigmaPt[5];
	
	RooRealVar* rv_Mass[4][2];
	RooArgSet*  rastmp_Mass[4][2];
	RooDataSet* Data_Mass[4][2];

	TH1F* MassErrDistribution[3];
	TH1F* REFITTEDMassErrDistribution[3];
	TH1F* DeltaRDistribution[5][3];
	TH3F* Pt[5][3];
	TH3F* Pixel[5][3];
	TH1F* Inclusive[3];
	TH1F* FinalSigma_GAUSS[3];
	TH1F* FinalSigma_DSCB[3];

	TH1F* PixelSigma_GAUSS[3];
	TH1F* PixelSigma_DSCB[3];

	TH1F* FinalSigma3_GAUSS[3][6];
	TH1F* PixelSigma3_GAUSS[3];
	
	TH1F* FinalMean_GAUSS[3];
	TH1F* FinalMean_DSCB[3];
	
	TH2F* Vertex = new TH2F("Vertex", "Vertex", 400, -0.2, 0.2, 400, -0.2, 0.2);				
	TH1F* Vertex_Z = new TH1F("Vertex_Z", "Vertex_Z", 300, -30, 30);
	TH2F* BeamSpot = new TH2F("BeamSpot", "BeamSpot", 400, -0.2, 0.2, 400, -0.2, 0.2);				
	TH1F* BeamSpot_Z = new TH1F("BeamSpot_Z", "BeamSpot_Z", 100, -10, 10);
	TH2F* VX[5];
	TH2F* VX_BS[5];
	TH1F* VX_Z[5];
	TH1F* VX_BS_Z[5];
	
	TH2F* Mass_Base_vs_VXBS = new TH2F("Mass: Base vs VX+BS", "Mass: Base vs VX+BS", 100, 105, 140, 100, 105, 140);
	TH2F* DeltaMass_vs_Base = new TH2F("#DeltaMass vs Base", "#DeltaMass vs Base", 100, 105, 140, 200, -10, 10);
	TH1F* DeltaMass = new TH1F("#DeltaMass", "#DeltaMass", 400, -2, 2);
	RooRealVar* rv_DeltaMass = new RooRealVar("#DeltaMass", "#DeltaMass", 400, -2, 2);
	RooArgSet*  rastmp_DeltaMass = new RooArgSet(*rv_DeltaMass);
	RooDataSet* Data_DeltaMass = new RooDataSet("#DeltaMass", "#DeltaMass", *rastmp_DeltaMass);


	TH2F* refitted_Mass_Base_vs_VXBS = new TH2F("REFITTED Mass: Base vs VX+BS", "REFITTED Mass: Base vs VX+BS", 100, 105, 140, 100, 105, 140);
	TH2F* refitted_DeltaMass_vs_Base = new TH2F("#DeltaREFITTEDMass vs Base", "#DeltaREFITTEDMass vs Base", 100, 105, 140, 200, -10, 10);
	TH1F* refitted_DeltaMass = new TH1F("#DeltaREFITTEDMass", "#DeltaREFITTEDMass", 200, -10, 10);
	RooRealVar* rv_DeltaRefittedMass = new RooRealVar("#DeltaRefittedMass", "#DeltaRefittedMass", 400, -2, 2);
	RooArgSet*  rastmp_DeltaRefittedMass = new RooArgSet(*rv_DeltaRefittedMass);
	RooDataSet* Data_DeltaRefittedMass = new RooDataSet("#DeltaRefittedMass", "#DeltaRefittedMass", *rastmp_DeltaRefittedMass);
	
	float integral_70[3][2] = {0};
	float integral_70_105[3][2]= {0};
	float integral_105_140[3][2] = {0};
	float integral_105_130[3][2] = {0};
	float integral_120_130[3][2] = {0};
	float integral_118_130[3][2] = {0};
	float integral_140_inf[3][2] = {0};

	std::vector<TString> lepton_type;
	lepton_type.push_back("GEN");	
// 	lepton_type.push_back("Selected_");
	lepton_type.push_back("Base");
	lepton_type.push_back("VX");
	lepton_type.push_back("VX+BS");	
	
	std::vector<Double_t> genReco_bins;
	for(int i = 0; i < 2000; i++)
		genReco_bins.push_back(-0.25 + 0.00025*i);
	genReco_bins.push_back(0.25);
// 		genReco_bins.push_back(-0.5 + 0.01*i);
// 	genReco_bins.push_back(0.5);

	std::vector<Double_t> pixel_bins;
	pixel_bins.push_back(0);
	// 	pixel_bins.push_back(1);
	pixel_bins.push_back(2);
	pixel_bins.push_back(3);
	pixel_bins.push_back(4);
	if(year != "2016") pixel_bins.push_back(5);
	if(year != "2016") pixel_bins.push_back(10);
// 	pixel_bins.push_back(20);

	RooRealVar* rv_ptRecoGen[eta_bins.size()][pT_bins.size()][lepton_type.size()-1]; 
	RooArgSet*  rastmp_ptRecoGen[eta_bins.size()][pT_bins.size()][lepton_type.size()-1]; 
	RooDataSet* Data_ptRecoGen[eta_bins.size()][pT_bins.size()][lepton_type.size()-1]; 

	RooRealVar* rv_DeltaMassOverMass[lepton_type.size()-1];
	RooArgSet*  rastmp_DeltaMassOverMass[lepton_type.size()-1];
	RooDataSet* Data_DeltaMassOverMass[lepton_type.size()-1];

	RooRealVar* rv_shift[eta_bins.size()][pT_bins.size()];
	RooArgSet*  rastmp_shift[eta_bins.size()][pT_bins.size()];
	RooDataSet* Data_shift[eta_bins.size()][pT_bins.size()];
	TH1F* FinalShift = new TH1F("Baseline and VX+BS Shift", "Baseline and VX+BS Shift", 23, 0, 23);

    
    float num_sample = 1;
    if(channel == "Full") num_sample = 7;
    
	std::vector<TString> ParodMode_File;
	ProdMode_File.clear();
	ProdMode_File.push_back("GluGluHToZZTo4L");
	ProdMode_File.push_back("VBF_HToZZTo4L");
	ProdMode_File.push_back("WplusH_HToZZTo4L");
	ProdMode_File.push_back("WminusH_HToZZTo4L");
	ProdMode_File.push_back("ZH_HToZZ_4LFilter");
	ProdMode_File.push_back("ttH_HToZZ_4LFilter");
	ProdMode_File.push_back("bbH_HToZZTo4L");

	std::vector<float> XS;
	XS.clear();
	XS.push_back(0.01333521);
	XS.push_back(0.001038159);
	XS.push_back(0.000146235);
	XS.push_back(0.0002305562);
	XS.push_back(0.000662058);
	XS.push_back(0.0003901903);
	XS.push_back(0.000134688);

	std::vector<float> numberEvents;
	if(year == "2016"){
		numberEvents.push_back(992224);
		numberEvents.push_back(500000);
		numberEvents.push_back(277590);
		numberEvents.push_back(197930);
		numberEvents.push_back(487840);
		numberEvents.push_back(492346);
		numberEvents.push_back(999224);
	}
	if(year == "2017"){
		numberEvents.push_back(965198);
		numberEvents.push_back(1000000);
		numberEvents.push_back(33400);
		numberEvents.push_back(195600);
		numberEvents.push_back(494864);
		numberEvents.push_back(136114);
		numberEvents.push_back(499392);
	}
	if(year == "2018"){
		numberEvents.push_back(958000);
		numberEvents.push_back(500000);
		numberEvents.push_back(299100);
		numberEvents.push_back(200000);
		numberEvents.push_back(466101);
		numberEvents.push_back(483860);
		numberEvents.push_back(488200);
	}

	int run_over_scale = 2;
	if(channel == "Full") run_over_scale = 1;

	for(int m = 0; m < run_over_scale; m++){	
		gStyle->SetOptStat(0);
		TString histo_name;
		for(int t = 0; t < 4; t++){
			for(int lep = 0; lep < 5; lep++){
				if(lep == 0){
					histo_name = lepton_type.at(t) + "Pt";					
					PtDistribution[lep][t] = new TH1F(histo_name, histo_name, 200, 0, 200);
					histo_name = lepton_type.at(t) + "SigmaPtOverPt";					
					SigmaPtOverPtDistribution[lep][t] = new TH1F(histo_name, histo_name, 100, 0, 0.1);
					histo_name = lepton_type.at(t) + "DeltaPt";					
					DeltaPtDistribution[lep][t] = new TH1F(histo_name, histo_name, 250, -0.25, 0.25);
					histo_name = lepton_type.at(t) + "Eta";					
					EtaDistribution[lep][t] = new TH1F(histo_name, histo_name, 48, -2.4, 2.4);
					histo_name = lepton_type.at(t) + "Phi";					
					PhiDistribution[lep][t] = new TH1F(histo_name, histo_name, 62, -3.14, 3.14);
					histo_name = lepton_type.at(t) + "Mass";					
					MassDistribution[t] = new TH1F(histo_name, histo_name, 100, 105, 140);
					histo_name = lepton_type.at(t) + "RefittedMass";					
					REFITTEDMassDistribution[t] = new TH1F(histo_name, histo_name, 100, 105, 140);
					histo_name = lepton_type.at(t) + "Weighted Mass";					
					W_MassDistribution[t] = new TH1F(histo_name, histo_name, 5000, 0, 5000);
					histo_name = lepton_type.at(t) + "WeightedRefittedMass";					
					W_REFITTEDMassDistribution[t] = new TH1F(histo_name, histo_name, 5000, 0, 5000);
					histo_name = lepton_type.at(t) + "MassErr";					
					MassErrDistribution[t] = new TH1F(histo_name, histo_name, 200, 0, 5);
					histo_name = lepton_type.at(t) + "RefittedMassErr";					
					REFITTEDMassErrDistribution[t] = new TH1F(histo_name, histo_name, 200, 0, 5);
					
					if(t == 0){
						histo_name = "SigmaPt_vsSigmaPt";					
						SigmaPt_vsSigmaPt[lep] = new TH2F(histo_name, histo_name, 100, 0, 0.1, 100, 0, 0.1);
					}

					if(t != 3){
						for(int k = 0; k < 2; k++){
							if(k == 0) histo_name = Form("MassDistribution_%s", lepton_type[t].Data());
							else histo_name = Form("RefittedMassDistribution_%s", lepton_type[t].Data());
							rv_Mass[t][k] = new RooRealVar(histo_name, histo_name, 105, 140);
							rastmp_Mass[t][k] = new RooArgSet(*rv_Mass[t][k]);
							Data_Mass[t][k] = new RooDataSet("Data" + histo_name, "Data" + histo_name, *rastmp_Mass[t][k]);
						}
					}

					if(t == 0 && m == 0){
						VX[lep] = new TH2F("Inclusive Common VX", "Inclusive Common VX", 400, -0.2, 0.2, 400, -0.2, 0.2);				
						VX_BS[lep] = new TH2F("Inclusive Common VX_BS", "Inclusive Common VX_BS", 400, -0.2, 0.2, 400, -0.2, 0.2);
						VX_Z[lep] = new TH1F("Inclusive Common VX: Z", "Inclusive Common VX: Z", 300, -30, 30);
						VX_BS_Z[lep] = new TH1F("Inclusive Common VX_BS: Z", "Inclusive Common VX_BS: Z", 300, -30, 30);
					}
				
					if(t != 0){
					
						histo_name = lepton_type.at(t) + "Inclusive DeltaM over M";	
						rv_DeltaMassOverMass[t-1] = new RooRealVar(histo_name, histo_name, genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1));
						rastmp_DeltaMassOverMass[t-1] = new RooArgSet(*rv_DeltaMassOverMass[t-1]);
						Data_DeltaMassOverMass[t-1] = new RooDataSet("Data" + histo_name, "Data" + histo_name,  *rastmp_DeltaMassOverMass[t-1]);
					
						histo_name = lepton_type.at(t) + "DeltaR";	
						DeltaRDistribution[0][t-1] = new TH1F("DeltaR_" + lepton_type.at(t), "DeltaR_" + lepton_type.at(t), 100, 0, 0.1);
						Inclusive[t-1] = new TH1F("Inclusive_" + lepton_type.at(t), "Inclusive_" + lepton_type.at(t), genReco_bins.size()-1, &genReco_bins[0]);

						FinalSigma_GAUSS[t-1] = new TH1F("Gauss_" + lepton_type.at(t), "Gauss_" + lepton_type.at(t), 23, 0, 23);
						FinalSigma_DSCB[t-1] = new TH1F("DSCB_" + lepton_type.at(t), "DSCB_" + lepton_type.at(t), 23, 0, 23);
						PixelSigma_GAUSS[t-1] = new TH1F("Gauss_" + lepton_type.at(t), "Gauss_" + lepton_type.at(t), 23, 0, 23);
						PixelSigma3_GAUSS[t-1] = new TH1F("Gauss3_" + lepton_type.at(t), "Gauss3_" + lepton_type.at(t), 23, 0, 23);
						PixelSigma_DSCB[t-1] = new TH1F("DSCB_" + lepton_type.at(t), "DSCB_" + lepton_type.at(t), 23, 0, 23);
						if(m == 0){
							FinalMean_GAUSS[t-1] = new TH1F("MeanGauss_" + lepton_type.at(t), "MeanGauss_" + lepton_type.at(t), 23, 0, 23);
							FinalMean_DSCB[t-1] = new TH1F("MeanDSCB_" + lepton_type.at(t), "MeanDSCB_" + lepton_type.at(t), 23, 0, 23);
						}
	
						for(int h = 0; h < 6; h++){
							TString nome_canvas;
							nome_canvas = "Gauss_" + lepton_type.at(t);
							nome_canvas += Form("%d", h+1);										
							FinalSigma3_GAUSS[t-1][h] = new TH1F(nome_canvas, nome_canvas, 23, 0, 23);
						}

						histo_name = lepton_type.at(t) + "Pt";					
						Pt[0][t-1] = new TH3F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0], pT_bins.size()-1, &pT_bins[0], genReco_bins.size()-1, &genReco_bins[0]);
						histo_name = lepton_type.at(t) + "Pixel";					
						Pixel[0][t-1] = new TH3F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0], pixel_bins.size()-1, &pixel_bins[0], genReco_bins.size()-1, &genReco_bins[0]);
					}

				}
				else{
					histo_name = lepton_type.at(t) + Form("Pt_Lep_%d", lep);
					PtDistribution[lep][t] = new TH1F(histo_name, histo_name, 200, 0, 200);
					histo_name = lepton_type.at(t) + Form("SigmaPtOverPt_Lep_%d", lep);
					SigmaPtOverPtDistribution[lep][t] = new TH1F(histo_name, histo_name, 100, 0, 0.1);
					histo_name = lepton_type.at(t) + Form("DeltaPt_Lep_%d", lep);
					DeltaPtDistribution[lep][t] = new TH1F(histo_name, histo_name, 250, -0.25, 0.25);
					histo_name = lepton_type.at(t) + Form("Eta_Lep_%d", lep);
					EtaDistribution[lep][t] = new TH1F(histo_name, histo_name, 48, -2.4, 2.4);
					histo_name = lepton_type.at(t) + Form("Phi_Lep_%d", lep);
					PhiDistribution[lep][t] = new TH1F(histo_name, histo_name, 62, -3.14, 3.14);

					if(t == 0){
						histo_name = Form("SigmaPt_vsSigmaPt_Lep_%d", lep); 
						SigmaPt_vsSigmaPt[lep] = new TH2F(histo_name, histo_name, 100, 0, 0.1, 100, 0, 0.1);
					}

					if(t != 0){
						histo_name = lepton_type.at(t) + Form("DeltaR_Lep_%d", lep);
						DeltaRDistribution[lep][t-1] = new TH1F("DeltaR_" + lepton_type.at(t), "DeltaR_" + lepton_type.at(t), 100, 0, 0.1);
						histo_name = lepton_type.at(t) + Form("Pt_Lep_%d", lep);
						Pt[lep][t-1] = new TH3F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0],  pT_bins.size()-1, &pT_bins[0], genReco_bins.size()-1, &genReco_bins[0]);
						histo_name = lepton_type.at(t) + Form("Pixel_Lep_%d", lep);
						Pixel[lep][t-1] = new TH3F(histo_name, histo_name, eta_bins.size()-1, &eta_bins[0],  pixel_bins.size()-1, &pixel_bins[0], genReco_bins.size()-1, &genReco_bins[0]);
					}
					if(t == 0 && m == 0){
						TString lepton = Form(" %d", lep);
						VX[lep] = new TH2F("Inclusive Common VX"+lepton, "Inclusive Common VX"+lepton, 400, -0.2, 0.2, 400, -0.2, 0.2);				
						VX_BS[lep] = new TH2F("Inclusive Common VX_BS"+lepton, "Inclusive Common VX_BS"+lepton, 400, -0.2, 0.2, 400, -0.2, 0.2);				
						VX_Z[lep] = new TH1F("Inclusive Common VX"+lepton+": Z", "Inclusive Common VX"+lepton+": Z", 300, -30, 30);
						VX_BS_Z[lep] = new TH1F("Inclusive Common VX_BS"+lepton+": Z", "Inclusive Common VX_BS"+lepton+": Z", 300, -30, 30);
					}

				}
			}		

			if(t != 0){
// 				std::cout<<Pt[0][0]->GetNbinsX()<<"\t"<<Pt[0][0]->GetNbinsY()<<std::endl;
				for(int x = 0; x <= Pt[0][0]->GetNbinsX(); x++){
					for(int y = 0; y <= Pt[0][0]->GetNbinsY(); y++){
						if(x == 0 && y == 0)
							histo_name = Form("UnbinnedInclusive_%s", lepton_type[t].Data());
						else if((x == 0 && y != 0) || (x != 0 && y == 0)) continue;
						else
							histo_name = Form("UnbinnedPt_%s_%.0f_%.0f_%s_%s", lepton_type[t].Data(), pT_bins.at(y-1), pT_bins.at(y), eta_bins_name[x-1].Data(), eta_bins_name[x].Data());
						rv_ptRecoGen[x][y][t-1] = new RooRealVar(histo_name, histo_name, genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1));
						rastmp_ptRecoGen[x][y][t-1] = new RooArgSet(*rv_ptRecoGen[x][y][t-1]);
						Data_ptRecoGen[x][y][t-1] = new RooDataSet("Data" + histo_name, "Data" + histo_name, *rastmp_ptRecoGen[x][y][t-1]);

						if(t == 1){
							if(x == 0 && y == 0)
								histo_name = Form("UnbinnedInclusive_%s", lepton_type[t].Data());
							else if((x == 0 && y != 0) || (x != 0 && y == 0)) continue;
							else
								histo_name = Form("UnbinnedPt_%.0f_%.0f_%s_%s", pT_bins.at(y-1), pT_bins.at(y), eta_bins_name[x-1].Data(), eta_bins_name[x].Data());
							rv_shift[x][y] = new RooRealVar(histo_name, histo_name, genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1));
							rastmp_shift[x][y] = new RooArgSet(*rv_shift[x][y]);
							Data_shift[x][y] = new RooDataSet("Shift" + histo_name, "Shift" + histo_name, *rastmp_shift[x][y]);
						}
					}
				}
			}

		}
		
	    
	for(int sample = 0; sample < num_sample; sample++){
	    
		if(channel == "Full")
		    filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/" + ProdMode_File.at(sample) + "_M125_" + year + "_skimmed.root";

	    if(channel == "ggF"){
	    	if(Base == 1) filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Production_10_2_18/Higgs_base/125/GluGluHToZZTo4L_M125_" + year + "_skimmed.root";
			else filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Production_10_2_18/Higgs_VX_BS/125/GluGluHToZZTo4L_M125_" + year + "_skimmed.root";
		}
	    if(channel == "ttH")
		    filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Vtx/ttH_HToZZ_4LFilter_M125_" + year + "_skimmed.root";
	    if(channel == "qqZZ")
		    filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Vtx/ZZTo4L_powheg_" + year + "_Vtx.root";
	    if(channel == "ggZZ")
		    filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Vtx/GluGluToContinToZZTo4mu_" + year + "_Vtx.root";

// 	    if(isDATA) filename = "/eos/user/f/ferrico/LXPLUS/ROOT_FILE_HZZ/Full_RunII/Vtx/Data_m2mu_" + year + ".root";

		std::cout<<filename<<std::endl;

	    _file0 = new TFile(filename);

		if(_file0)
			tree = (TTree*)_file0->Get("passedEvents");		
		else std::cout<<"ERROR could not find the file"<<std::endl;
		
		tree->SetBranchStatus("mass4l",1);               
		tree->SetBranchStatus("mass4lErr",1);            
		tree->SetBranchStatus("mass4lREFIT",1);          
		tree->SetBranchStatus("mass4lErrREFIT",1);       
		tree->SetBranchStatus("mass4mu",1);              
		tree->SetBranchStatus("mass2e2mu",1);            
                                                 
 		tree->SetBranchStatus("mass4l_vtx",1);       
 		tree->SetBranchStatus("mass4l_vtxFSR",1);       
 		tree->SetBranchStatus("mass4lErr_vtx",1);    
 		tree->SetBranchStatus("mass4lREFIT_vtx",1);  
 		tree->SetBranchStatus("mass4lErrREFIT_vtx",1);
//  		tree->SetBranchStatus("massH_vtx_chi2",1);   
                                                                                           
		tree->SetBranchStatus("mass4l_vtx_BS",1);        
 		tree->SetBranchStatus("mass4l_vtxFSR_BS",1);       
		tree->SetBranchStatus("mass4lErr_vtx_BS",1);     
		tree->SetBranchStatus("mass4lREFIT_vtx_BS",1);   
		tree->SetBranchStatus("mass4lErrREFIT_vtx_BS",1);
// 		tree->SetBranchStatus("massH_vtx_chi2_BS",1);    
 
		tree->SetBranchStatus("GENMH",1);        	
		tree->SetBranchStatus("GENlep_pt",1);                
		tree->SetBranchStatus("GENlep_eta",1);               
		tree->SetBranchStatus("GENlep_phi",1);               
		tree->SetBranchStatus("GENlep_mass",1);                
                                                     
		tree->SetBranchStatus("lepFSR_pt",1);                   
		tree->SetBranchStatus("lepFSR_eta",1);                  
		tree->SetBranchStatus("lepFSR_phi",1);                

		tree->SetBranchStatus("lep_pt",1);                   
		tree->SetBranchStatus("lep_pterr", 1);
		tree->SetBranchStatus("lep_eta",1);                  
		tree->SetBranchStatus("lep_phi",1);                
		tree->SetBranchStatus("lep_id",1);                
		tree->SetBranchStatus("lep_mass",1);                

		tree->SetBranchStatus("vtxLep_pt",1);                
		tree->SetBranchStatus("vtxLep_ptError", 1);
		tree->SetBranchStatus("vtxLep_eta",1);               
		tree->SetBranchStatus("vtxLep_phi",1);               
		tree->SetBranchStatus("vtxLep_mass",1);              
                                                     
		tree->SetBranchStatus("vtxLep_BS_pt",1);             
		tree->SetBranchStatus("vtxLep_BS_ptError", 1);
		tree->SetBranchStatus("vtxLep_BS_eta",1);            
		tree->SetBranchStatus("vtxLep_BS_phi",1);            
		tree->SetBranchStatus("vtxLep_BS_mass",1);           

		tree->SetBranchStatus("vtxLepFSR_pt",1);                
		tree->SetBranchStatus("vtxLepFSR_eta",1);               
		tree->SetBranchStatus("vtxLepFSR_phi",1);               
		tree->SetBranchStatus("vtxLepFSR_mass",1);              
                                                     
		tree->SetBranchStatus("vtxLepFSR_BS_pt",1);             
		tree->SetBranchStatus("vtxLepFSR_BS_eta",1);            
		tree->SetBranchStatus("vtxLepFSR_BS_phi",1);            
		tree->SetBranchStatus("vtxLepFSR_BS_mass",1);           

		tree->SetBranchStatus("passedFullSelection",1);                                                           
		tree->SetBranchStatus("finalState",1); 
		
		tree->SetBranchStatus("pTL1",1); 
		tree->SetBranchStatus("etaL1",1); 
		tree->SetBranchStatus("phiL1",1);               
		tree->SetBranchStatus("mL1",1);               

		tree->SetBranchStatus("pTL2",1); 
		tree->SetBranchStatus("etaL2",1); 
		tree->SetBranchStatus("phiL2",1);               
		tree->SetBranchStatus("mL3",1);               

		tree->SetBranchStatus("pTL3",1); 
		tree->SetBranchStatus("etaL3",1); 
		tree->SetBranchStatus("phiL3",1);               
		tree->SetBranchStatus("mL3",1);               

		tree->SetBranchStatus("pTL4",1); 
		tree->SetBranchStatus("etaL4",1); 
		tree->SetBranchStatus("phiL4",1);               
		tree->SetBranchStatus("mL4",1);               
                                                     
		tree->SetBranchStatus("lep_numberOfValidPixelHits",1);  

		tree->SetBranchStatus("PV_x",1);  
		tree->SetBranchStatus("PV_y",1);  
		tree->SetBranchStatus("PV_z",1);  
		tree->SetBranchStatus("BS_x",1);  
		tree->SetBranchStatus("BS_y",1);  
		tree->SetBranchStatus("BS_z",1);  

		tree->SetBranchStatus("commonPV_x",1);  
		tree->SetBranchStatus("commonPV_y",1);  
		tree->SetBranchStatus("commonPV_z",1);  
		tree->SetBranchStatus("commonBS_x",1);  
		tree->SetBranchStatus("commonBS_y",1);  
		tree->SetBranchStatus("commonBS_z",1);  

		tree->SetBranchStatus("eventWeight",1);  


		tree->SetBranchAddress("mass4l", &mass4l);
		tree->SetBranchAddress("mass4lREFIT", &mass4lREFIT);
		tree->SetBranchAddress("mass4mu", &mass4mu);
		tree->SetBranchAddress("mass2e2mu", &mass2e2mu);

		tree->SetBranchAddress("mass4l_vtx", &mass4l_vtx);
		tree->SetBranchAddress("mass4l_vtxFSR", &mass4l_vtxFSR);
		tree->SetBranchAddress("mass4lREFIT_vtx", &mass4lREFIT_vtx);

		tree->SetBranchAddress("mass4l_vtx_BS", &mass4l_vtx_BS);
		tree->SetBranchAddress("mass4l_vtxFSR_BS", &mass4l_vtxFSR_BS);
		tree->SetBranchAddress("mass4lREFIT_vtx_BS", &mass4lREFIT_vtx_BS);

		tree->SetBranchAddress("GENMH", &GENMH);        	

		tree->SetBranchAddress("GENlep_pt", &GENlep_pt);
		tree->SetBranchAddress("GENlep_eta", &GENlep_eta);
		tree->SetBranchAddress("GENlep_phi", &GENlep_phi);
		tree->SetBranchAddress("GENlep_mass", &GENlep_mass);

		tree->SetBranchAddress("lep_pt", &lep_pt);
		tree->SetBranchAddress("lep_pterr", &lep_pterr);
		tree->SetBranchAddress("lep_eta", &lep_eta);
		tree->SetBranchAddress("lep_phi", &lep_phi);
		tree->SetBranchAddress("lep_id", &lep_id);                
		tree->SetBranchAddress("lep_mass", &lep_mass);

		tree->SetBranchAddress("lepFSR_pt", &lepFSR_pt);
		tree->SetBranchAddress("lepFSR_eta", &lepFSR_eta);
		tree->SetBranchAddress("lepFSR_phi", &lepFSR_phi);
		tree->SetBranchAddress("lepFSR_mass", &lepFSR_mass);

		tree->SetBranchAddress("vtxLep_pt", &vtxLep_pt);
		tree->SetBranchAddress("vtxLep_ptError", &vtxLep_ptError);
		tree->SetBranchAddress("vtxLep_eta", &vtxLep_eta);
		tree->SetBranchAddress("vtxLep_phi", &vtxLep_phi);
		tree->SetBranchAddress("vtxLep_mass", &vtxLep_mass);

		tree->SetBranchAddress("vtxLep_BS_pt", &vtxLep_BS_pt);
		tree->SetBranchAddress("vtxLep_BS_ptError", &vtxLep_BS_ptError);
		tree->SetBranchAddress("vtxLep_BS_eta", &vtxLep_BS_eta);
		tree->SetBranchAddress("vtxLep_BS_phi", &vtxLep_BS_phi);
		tree->SetBranchAddress("vtxLep_BS_mass", &vtxLep_BS_mass);

		tree->SetBranchAddress("vtxLepFSR_pt", &vtxLepFSR_pt);
		tree->SetBranchAddress("vtxLepFSR_eta", &vtxLepFSR_eta);
		tree->SetBranchAddress("vtxLepFSR_phi", &vtxLepFSR_phi);
		tree->SetBranchAddress("vtxLepFSR_mass", &vtxLepFSR_mass);

		tree->SetBranchAddress("vtxLepFSR_BS_pt", &vtxLepFSR_BS_pt);
		tree->SetBranchAddress("vtxLepFSR_BS_eta", &vtxLepFSR_BS_eta);
		tree->SetBranchAddress("vtxLepFSR_BS_phi", &vtxLepFSR_BS_phi);
		tree->SetBranchAddress("vtxLepFSR_BS_mass", &vtxLepFSR_BS_mass);

		tree->SetBranchAddress("passedFullSelection", &passedFullSelection);
		tree->SetBranchAddress("finalState", &finalState);

		tree->SetBranchAddress("pTL1", &pTL1);
		tree->SetBranchAddress("etaL1", &etaL1);
		tree->SetBranchAddress("phiL1", &phiL1);
		tree->SetBranchAddress("mL1", &mL1);

		tree->SetBranchAddress("pTL2", &pTL2);
		tree->SetBranchAddress("etaL2", &etaL2);
		tree->SetBranchAddress("phiL2", &phiL2);
		tree->SetBranchAddress("mL2", &mL2);

		tree->SetBranchAddress("pTL3", &pTL3);
		tree->SetBranchAddress("etaL3", &etaL3);
		tree->SetBranchAddress("phiL3", &phiL3);
		tree->SetBranchAddress("mL3", &mL3);

		tree->SetBranchAddress("pTL4", &pTL4);
		tree->SetBranchAddress("etaL4", &etaL4);
		tree->SetBranchAddress("phiL4", &phiL4);
		tree->SetBranchAddress("mL4", &mL4);
		
		tree->SetBranchAddress("lep_numberOfValidPixelHits", &lep_numberOfValidPixelHits);

		tree->SetBranchAddress("PV_x", &PV_x);
		tree->SetBranchAddress("PV_y", &PV_y);
		tree->SetBranchAddress("PV_z", &PV_z);
		tree->SetBranchAddress("BS_x", &BS_x);
		tree->SetBranchAddress("BS_y", &BS_y);
		tree->SetBranchAddress("BS_z", &BS_z);

		tree->SetBranchAddress("commonPV_x", &commonPV_x);
		tree->SetBranchAddress("commonPV_y", &commonPV_y);
		tree->SetBranchAddress("commonPV_z", &commonPV_z);
		tree->SetBranchAddress("commonBS_x", &commonBS_x);
		tree->SetBranchAddress("commonBS_y", &commonBS_y);
		tree->SetBranchAddress("commonBS_z", &commonBS_z);
 
		tree->SetBranchAddress("mass4lErr", &mass4lErr);  
		tree->SetBranchAddress("mass4lErrREFIT", &mass4lErrREFIT);  
		tree->SetBranchAddress("mass4lErr_vtx", &mass4lErr_vtx);  
		tree->SetBranchAddress("mass4lErr_vtx_BS", &mass4lErr_vtx_BS);  
		tree->SetBranchAddress("mass4lErrREFIT_vtx", &mass4lErrREFIT_vtx);  
		tree->SetBranchAddress("mass4lErrREFIT_vtx_BS", &mass4lErrREFIT_vtx_BS);  

		tree->SetBranchAddress("eventWeight", &weight);  

		
		Long64_t nentries_MU = tree->GetEntries();

		if(year == "2017" && channel == "qqZZ") 	nentries_MU = 25000000;	
		
		std::cout<<nentries_MU<<std::endl;	
		
		float minDR;
		if(!isDR) minDR = 1000;
		else minDR = 0.02;
				
		float new_pt;
		int keta, kpt, bin;
		
		float lumi;
				
 		for(int entry = 0; entry < nentries_MU; entry++){
// 		for(int entry = 0; entry < (int)nentries_MU/10; entry++){
// 		for(int entry = 0; entry < (int)nentries_MU/100; entry++){
// 		for(int entry = 0; entry < 1000; entry++){

			tree->GetEntry(entry);
									
			float deltaEta, deltaPhi, DR, DR2; 
			int index;

			if(entry % 100000 == 0)       
				std::cout<<entry<<" --- Dentro il TREE --- "<<m<<std::endl;      	

				if(!passedFullSelection) continue;
				if(FS != finalState) continue;
												
// 				if(channel == "ttH" && lepFSR_pt->size() > 4) continue;
// 				if(channel == "ttH" && lepFSR_pt->size() <= 4) continue;
				
				if(year == "2016") lumi = 35920;
				if(year == "2017") lumi = 41530;
				if(year == "2018") lumi = 59470;
				
				weight = weight * XS.at(sample) * lumi / numberEvents[sample];
				
				
				std::vector<float> PT;
				PT.push_back(pTL1);
				PT.push_back(pTL2);
				PT.push_back(pTL3);
				PT.push_back(pTL4);
// 				PT.push_back(lepFSR_pt->at(0));
// 				PT.push_back(lepFSR_pt->at(1));
// 				PT.push_back(lepFSR_pt->at(2));
// 				PT.push_back(lepFSR_pt->at(3));
				std::vector<float> ETA;
				ETA.push_back(etaL1);
				ETA.push_back(etaL2);
				ETA.push_back(etaL3);
				ETA.push_back(etaL4);
				std::vector<float> PHI;
				PHI.push_back(phiL1);
				PHI.push_back(phiL2);
				PHI.push_back(phiL3);
				PHI.push_back(phiL4);
				std::vector<float> MASS;
				MASS.push_back(mL1);
				MASS.push_back(mL2);
				MASS.push_back(mL3);
				MASS.push_back(mL4);
												
				REFITTEDMassDistribution[0]->Fill(GENMH);
				REFITTEDMassDistribution[1]->Fill(mass4lREFIT);
				REFITTEDMassDistribution[2]->Fill(mass4lREFIT_vtx);
				REFITTEDMassDistribution[3]->Fill(mass4lREFIT_vtx_BS);

				W_REFITTEDMassDistribution[1]->Fill(mass4lREFIT, weight);
				W_REFITTEDMassDistribution[2]->Fill(mass4lREFIT_vtx, weight);
				W_REFITTEDMassDistribution[3]->Fill(mass4lREFIT_vtx_BS, weight);

				REFITTEDMassErrDistribution[0]->Fill(mass4lErrREFIT);
				REFITTEDMassErrDistribution[1]->Fill(mass4lErrREFIT_vtx);
				REFITTEDMassErrDistribution[2]->Fill(mass4lErrREFIT_vtx_BS);
				
				refitted_Mass_Base_vs_VXBS->Fill(mass4lREFIT, mass4lREFIT_vtx_BS);
				refitted_DeltaMass_vs_Base->Fill(mass4lREFIT, mass4lREFIT - mass4lREFIT_vtx_BS);
				refitted_DeltaMass->Fill(mass4lREFIT - mass4lREFIT_vtx_BS);
				
				rv_DeltaRefittedMass->setVal(mass4lREFIT - mass4lREFIT_vtx_BS);
				Data_DeltaRefittedMass->add(*rastmp_DeltaRefittedMass);

				rv_Mass[0][1]->setVal(mass4lREFIT);
				Data_Mass[0][1]->add(*rastmp_Mass[0][1]);
				rv_Mass[1][1]->setVal(mass4lREFIT_vtx);
				Data_Mass[1][1]->add(*rastmp_Mass[1][1]);
				rv_Mass[2][1]->setVal(mass4lREFIT_vtx_BS);
				Data_Mass[2][1]->add(*rastmp_Mass[2][1]);

				TLorentzVector tmp, refit, refit_vtx, refit_vtx_BS;
// 					for(int lep = 0; lep < 4; lep++){
// 						tmp.SetPtEtaPhiM(lepFSR_pt->at(lep), lepFSR_eta->at(lep), lepFSR_phi->at(lep), lepFSR_mass->at(lep));
// 							refit += tmp;														
// 						tmp.SetPtEtaPhiM(vtxLep_pt->at(lep), vtxLep_eta->at(lep), vtxLep_phi->at(lep), vtxLep_mass->at(lep));
// 							refit_vtx += tmp;														
// 						tmp.SetPtEtaPhiM(vtxLep_BS_pt->at(lep), vtxLep_BS_eta->at(lep), vtxLep_BS_phi->at(lep), vtxLep_BS_mass->at(lep));
// 							refit_vtx_BS += tmp;														
// 					}
// 					std::cout<<refit.M()-mass4l<<"\t"<<refit_vtx.M()-mass4l_vtx<<"\t"<<refit_vtx_BS.M()-mass4l_vtx_BS<<std::endl;
				
				
				if(m == 1){
					for(int lep = 0; lep < 4; lep++){
// 						tmp.SetPtEtaPhiM(lep_pt->at(lep), lep_eta->at(lep), lep_phi->at(lep), lep_mass->at(lep));
// 						tmp.SetPtEtaPhiM(PT.at(lep), ETA.at(lep), PHI.at(lep), MASS.at(lep));
						for(int pt = 1; pt < pT_bins.size(); pt++){
// 							if(lepFSR_pt->at(lep) < pT_bins.at(pt)){
							if(lep_pt->at(lep) < pT_bins.at(pt)){
								kpt = pt;
								break;
							}
						}
						for(int eta = 1; eta < eta_bins.size(); eta++){
// 							if(fabs(lepFSR_eta->at(lep)) < fabs(eta_bins.at(eta))){
							if(fabs(lep_eta->at(lep)) < fabs(eta_bins.at(eta))){
								keta = eta;
								break;
							}
						}								
						bin = 2 + (kpt - 1) + (pT_bins.size()-1) * (keta - 1);
						newWeight = FinalMean_GAUSS[0]->GetBinContent((int)bin + 1);
// 						std::cout<<"BBB:\t\t"<<"\tbaseline = "<<kpt<<"\t"<<lep_pt->at(lep)<<"\t"<<fabs(lep_eta->at(lep))<<"\t"<<keta<<"\t"<<bin<<"\t"<<newWeight<<std::endl;
						new_pt = lepFSR_pt->at(lep)*(1-newWeight);
						tmp.SetPtEtaPhiM(new_pt, lepFSR_eta->at(lep), lepFSR_phi->at(lep), lepFSR_mass->at(lep));
						refit += tmp;					



						for(int pt = 1; pt < pT_bins.size(); pt++){
							if(vtxLep_pt->at(lep) < pT_bins.at(pt)){
								kpt = pt;
								break;
							}
						}
						for(int eta = 1; eta < eta_bins.size(); eta++){
							if(fabs(vtxLep_eta->at(lep)) < fabs(eta_bins.at(eta))){
								keta = eta;
								break;
							}
						}								
						bin = 2 + (kpt - 1) + (pT_bins.size()-1) * (keta - 1);
						newWeight = FinalMean_GAUSS[1]->GetBinContent((int)bin + 1);
						new_pt = vtxLep_pt->at(lep)*(1-newWeight);
						tmp.SetPtEtaPhiM(new_pt, vtxLepFSR_eta->at(lep), vtxLepFSR_phi->at(lep), vtxLepFSR_mass->at(lep));
						refit_vtx += tmp;



						for(int pt = 1; pt < pT_bins.size(); pt++){
// 							if(vtxLepFSR_BS_pt->at(lep) < pT_bins.at(pt)){
							if(vtxLep_BS_pt->at(lep) < pT_bins.at(pt)){
								kpt = pt;
								break;
							}
						}
						for(int eta = 1; eta < eta_bins.size(); eta++){
// 							if(fabs(vtxLepFSR_BS_eta->at(lep)) < fabs(eta_bins.at(eta))){
							if(fabs(vtxLep_BS_eta->at(lep)) < fabs(eta_bins.at(eta))){
								keta = eta;
								break;
							}
						}								
						bin = 2 + (kpt - 1) + (pT_bins.size()-1) * (keta - 1);
						newWeight = FinalMean_GAUSS[2]->GetBinContent((int)bin + 1);
// 						std::cout<<"BBB:\t\tVX+BS = "<<kpt<<"\t"<<vtxLep_BS_pt->at(lep)<<"\t"<<fabs(vtxLep_BS_eta->at(lep))<<"\t"<<keta<<"\t"<<bin<<"\t"<<newWeight<<std::endl;
						new_pt = vtxLepFSR_BS_pt->at(lep)*(1-newWeight);
						tmp.SetPtEtaPhiM(new_pt, vtxLepFSR_BS_eta->at(lep), vtxLepFSR_BS_phi->at(lep), vtxLepFSR_BS_mass->at(lep));
						refit_vtx_BS += tmp;
					}
				
					mass4l = refit.M();
					mass4l_vtxFSR = refit_vtx.M();
					mass4l_vtxFSR_BS = refit_vtx_BS.M();
				}

				rv_DeltaMassOverMass[0]->setVal((mass4l - 125)/125);
				Data_DeltaMassOverMass[0]->add(*rastmp_DeltaMassOverMass[0]);
				rv_DeltaMassOverMass[1]->setVal((mass4l_vtxFSR - 125)/125);
				Data_DeltaMassOverMass[1]->add(*rastmp_DeltaMassOverMass[1]);
				rv_DeltaMassOverMass[2]->setVal((mass4l_vtxFSR_BS - 125)/125);
				Data_DeltaMassOverMass[2]->add(*rastmp_DeltaMassOverMass[2]);
				
// 				std::cout<<(mass4l - 125)/125<<"\t\t"<<(mass4l_vtxFSR - 125)/125<<"\t\t"<<(mass4l_vtxFSR_BS - 125)/125<<std::endl;


								
				MassDistribution[0]->Fill(GENMH);
				MassDistribution[1]->Fill(mass4l);
				MassDistribution[2]->Fill(mass4l_vtxFSR);
				MassDistribution[3]->Fill(mass4l_vtxFSR_BS);

				W_MassDistribution[1]->Fill(mass4l, weight);
				W_MassDistribution[2]->Fill(mass4l_vtxFSR, weight);
				W_MassDistribution[3]->Fill(mass4l_vtxFSR_BS, weight);

				MassErrDistribution[0]->Fill(mass4lErr);
				MassErrDistribution[1]->Fill(mass4lErr_vtx);
				MassErrDistribution[2]->Fill(mass4lErr_vtx_BS);
								
				Mass_Base_vs_VXBS->Fill(mass4l, mass4l_vtxFSR_BS);
				DeltaMass_vs_Base->Fill(mass4l, mass4l - mass4l_vtxFSR_BS);
				DeltaMass->Fill(mass4l - mass4l_vtxFSR_BS);

				rv_DeltaMass->setVal(mass4l - mass4l_vtxFSR_BS);
				Data_DeltaMass->add(*rastmp_DeltaMass);
			
				rv_Mass[0][0]->setVal(mass4l);
				Data_Mass[0][0]->add(*rastmp_Mass[0][0]);
				rv_Mass[1][0]->setVal(mass4l_vtxFSR);
				Data_Mass[1][0]->add(*rastmp_Mass[1][0]);
				rv_Mass[2][0]->setVal(mass4l_vtxFSR_BS);
				Data_Mass[2][0]->add(*rastmp_Mass[2][0]);
						
				Vertex->Fill(PV_x, PV_y);
				Vertex_Z->Fill(PV_z);
				BeamSpot->Fill(BS_x, BS_y);
				BeamSpot_Z->Fill(BS_z);
				
				int n_lep;
								
				for(int lep = 0; lep < 4; lep++){
				
				
// 					if(FS == 3 || FS == 4){
// 						std::cout<<PT.at(lep)<<"\t"<<vtxLep_BS_pt->at(lep)<<std::endl;
// 						std::cout<<ETA.at(lep)<<"\t"<<vtxLep_BS_eta->at(lep)<<std::endl;
// 						std::cout<<PHI.at(lep)<<"\t"<<vtxLep_BS_phi->at(lep)<<std::endl;
// 					}
				
									
					if(mass4l > 105 && mass4l < 140 && PT.at(lep) > 5 && PT.at(lep) < pT_bins.at(pT_bins.size()-1) && fabs(ETA.at(lep)) < 2.4){
						if(mass4l_vtxFSR_BS > 105 && mass4l_vtxFSR_BS < 140 && vtxLep_BS_pt->at(lep) > 5 && vtxLep_BS_pt->at(lep) < pT_bins.at(pT_bins.size()-1) && fabs(vtxLep_BS_eta->at(lep)) < 2.4){ 
							SigmaPt_vsSigmaPt[0]->Fill(lep_pterr->at(lep)/PT.at(lep), vtxLep_BS_ptError->at(lep)/vtxLep_BS_pt->at(lep));
							SigmaPt_vsSigmaPt[lep+1]->Fill(lep_pterr->at(lep)/PT.at(lep), vtxLep_BS_ptError->at(lep)/vtxLep_BS_pt->at(lep));
						}
					}
					
					
					
					
					if(mass4l > 105 && mass4l < 140 && PT.at(lep) > 5 && PT.at(lep) < pT_bins.at(pT_bins.size()-1) && fabs(ETA.at(lep)) < 2.4){ 
// 						if(fabs(lep_id->at(lep)) == 13){
						if(13 == 13){
							DR = 999;
							for(int lep2 = 0; lep2 < GENlep_pt->size(); lep2++){
								deltaEta = ETA.at(lep) - GENlep_eta->at(lep2);
								deltaPhi = PHI.at(lep) - GENlep_phi->at(lep2);
								DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
								if(DR2 < DR){
									DR = DR2;
									index = lep2;
								}							
							}

							if(DR < minDR){
								for(int pt = 1; pt < pT_bins.size(); pt++){
									if(PT.at(lep) < pT_bins.at(pt)){
										kpt = pt;
										break;
									}
								}
								for(int eta = 1; eta < eta_bins.size(); eta++){
									if(fabs(ETA.at(lep)) < fabs(eta_bins.at(eta))){
										keta = eta;
										break;
									}
								}
								
								if(m == 0)
									newWeight = 0;
								if(m == 1){
									bin = 2 + (kpt - 1) + (pT_bins.size()-1) * (keta - 1);
									newWeight = FinalMean_GAUSS[0]->GetBinContent((int)bin + 1);
								}
								
								float new_pt = PT.at(lep)*(1-newWeight);
						
// 								PtDistribution[0][1]->Fill(PT.at(lep));
// 								PtDistribution[lep+1][1]->Fill(PT.at(lep));
// 								SigmaPtOverPtDistribution[0][1]->Fill(lep_pterr->at(lep)/PT.at(lep));
// 								SigmaPtOverPtDistribution[lep+1][1]->Fill(lep_pterr->at(lep)/PT.at(lep));
								PtDistribution[0][1]->Fill(new_pt);
								PtDistribution[lep+1][1]->Fill(new_pt);
								SigmaPtOverPtDistribution[0][1]->Fill(lep_pterr->at(lep)/new_pt);
								SigmaPtOverPtDistribution[lep+1][1]->Fill(lep_pterr->at(lep)/new_pt);
								EtaDistribution[0][1]->Fill(ETA.at(lep));
								EtaDistribution[lep+1][1]->Fill(ETA.at(lep));
								PhiDistribution[0][1]->Fill(PHI.at(lep));
								PhiDistribution[lep+1][1]->Fill(PHI.at(lep));
							
								PtDistribution[0][0]->Fill(GENlep_pt->at(index));
								PtDistribution[lep+1][0]->Fill(GENlep_pt->at(index));
								EtaDistribution[0][0]->Fill(GENlep_eta->at(index));
								EtaDistribution[lep+1][0]->Fill(GENlep_eta->at(index));
								PhiDistribution[0][0]->Fill(GENlep_phi->at(index));
								PhiDistribution[lep+1][0]->Fill(GENlep_phi->at(index));
								
								DeltaRDistribution[0][0]->Fill(DR);
								DeltaRDistribution[lep+1][0]->Fill(DR);

								DeltaPtDistribution[0][0]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
								DeltaPtDistribution[lep+1][0]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));

								if(fabs(lep_id->at(lep)) == 13){
									Pt[0][0]->Fill(fabs(ETA.at(lep)), new_pt, (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Pt[lep+1][0]->Fill(fabs(ETA.at(lep)), new_pt, (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Inclusive[0]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
// 									if(FS == 1) Pixel[0][0]->Fill(fabs(ETA.at(lep)), lep_numberOfValidPixelHits->at(lep), (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
// 									if(FS == 1) Pixel[lep+1][0]->Fill(fabs(ETA.at(lep)), lep_numberOfValidPixelHits->at(lep), (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));

									rv_ptRecoGen[keta][kpt][0]->setVal((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Data_ptRecoGen[keta][kpt][0]->add(*rastmp_ptRecoGen[keta][kpt][0]);
									rv_ptRecoGen[0][0][0]->setVal((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Data_ptRecoGen[0][0][0]->add(*rastmp_ptRecoGen[0][0][0]);

									if(mass4l_vtxFSR_BS > 105 && mass4l_vtxFSR_BS < 140 && vtxLep_BS_pt->at(lep) > 5 && vtxLep_BS_pt->at(lep) < pT_bins.at(pT_bins.size()-1) && fabs(vtxLep_BS_eta->at(lep)) < 2.4){ 
										rv_shift[keta][kpt]->setVal((vtxLep_BS_pt->at(lep) - PT.at(lep))/PT.at(lep));
										Data_shift[keta][kpt]->add(*rastmp_shift[keta][kpt]);
										rv_shift[0][0]->setVal((vtxLep_BS_pt->at(lep) - PT.at(lep))/PT.at(lep));
										Data_shift[0][0]->add(*rastmp_shift[0][0]);
									}
								}
								
							}

						}
						else{
							PtDistribution[0][1]->Fill(-999);
							PtDistribution[lep+1][1]->Fill(-999);
							EtaDistribution[0][1]->Fill(-999);
							EtaDistribution[lep+1][1]->Fill(-999);
							PhiDistribution[0][1]->Fill(-999);
							PhiDistribution[lep+1][1]->Fill(-999);
							DeltaRDistribution[0][0]->Fill(-999);
							DeltaRDistribution[lep+1][0]->Fill(-999);
							Pt[0][0]->Fill(-999, -999, -999);
							Pt[lep+1][0]->Fill(-999, -999, -999);
							Inclusive[0]->Fill(-999);
							if(FS == 1) Pixel[0][0]->Fill(-999, -999, -999);
							if(FS == 1) Pixel[lep+1][0]->Fill(-999, -999, -999);
						}
					}

					if(mass4l_vtxFSR > 105 && mass4l_vtxFSR < 140 && vtxLep_pt->at(lep) > 5 && vtxLep_pt->at(lep) < pT_bins.at(pT_bins.size()-1) && fabs(vtxLep_eta->at(lep)) < 2.4){ 
// 						if(fabs(lep_id->at(lep)) == 13){
						if(13 == 13){
							DR = 999;
							for(int lep2 = 0; lep2 < GENlep_pt->size(); lep2++){
								deltaEta = vtxLep_eta->at(lep) - GENlep_eta->at(lep2);
								deltaPhi = vtxLep_phi->at(lep) - GENlep_phi->at(lep2);
								DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
								DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
								if(DR2 < DR){
									DR = DR2;
									index = lep2;
								}							
							}

							if(DR < minDR){
								int keta, kpt;
								for(int pt = 1; pt < pT_bins.size(); pt++){
									if(vtxLep_pt->at(lep) < pT_bins.at(pt)){
										kpt = pt;
										break;
									}
								}
								for(int eta = 1; eta < eta_bins.size(); eta++){
									if(fabs(vtxLep_eta->at(lep)) < fabs(eta_bins.at(eta))){
										keta = eta;
										break;
									}
								}

								if(m == 0)
									newWeight = 0;
								if(m == 1){
									bin = 2 + (kpt - 1) + (pT_bins.size()-1) * (keta - 1);
									newWeight = FinalMean_GAUSS[1]->GetBinContent((int)bin + 1);
								}
								
								float new_pt = vtxLep_pt->at(lep)*(1-newWeight);
								
								PtDistribution[0][2]->Fill(vtxLep_pt->at(lep)*(1-newWeight));
								PtDistribution[lep+1][2]->Fill(vtxLep_pt->at(lep)*(1-newWeight));
								SigmaPtOverPtDistribution[0][2]->Fill(vtxLep_ptError->at(lep)/vtxLep_pt->at(lep)*(1-newWeight));
								SigmaPtOverPtDistribution[lep+1][2]->Fill(vtxLep_ptError->at(lep)/vtxLep_pt->at(lep)*(1-newWeight));
								EtaDistribution[0][2]->Fill(vtxLep_eta->at(lep));
								EtaDistribution[lep+1][2]->Fill(vtxLep_eta->at(lep));
								PhiDistribution[0][2]->Fill(vtxLep_phi->at(lep));
								PhiDistribution[lep+1][2]->Fill(vtxLep_phi->at(lep));

								DeltaRDistribution[0][1]->Fill(DR);
								DeltaRDistribution[lep+1][1]->Fill(DR);

								DeltaPtDistribution[0][1]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
								DeltaPtDistribution[lep+1][1]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));

								if(fabs(lep_id->at(lep)) == 13){
									Pt[0][1]->Fill(fabs(vtxLep_eta->at(lep)), new_pt, (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Pt[lep+1][1]->Fill(fabs(vtxLep_eta->at(lep)), new_pt, (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Inclusive[1]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
// 									if(FS == 1) Pixel[0][1]->Fill(fabs(vtxLep_eta->at(lep)), lep_numberOfValidPixelHits->at(lep), (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
// 									if(FS == 1) Pixel[lep+1][1]->Fill(fabs(vtxLep_eta->at(lep)), lep_numberOfValidPixelHits->at(lep), (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));

									rv_ptRecoGen[keta][kpt][1]->setVal((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Data_ptRecoGen[keta][kpt][1]->add(*rastmp_ptRecoGen[keta][kpt][1]);
									rv_ptRecoGen[0][0][1]->setVal((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Data_ptRecoGen[0][0][1]->add(*rastmp_ptRecoGen[0][0][1]);
								}

								VX[0]->Fill(commonPV_x->at(lep), commonPV_y->at(lep));
								VX[lep+1]->Fill(commonPV_x->at(lep), commonPV_y->at(lep));
								VX_Z[0]->Fill(commonPV_z->at(lep));
								VX_Z[lep+1]->Fill(commonPV_z->at(lep));
							}
						}
						else{
							PtDistribution[0][2]->Fill(-999);
							PtDistribution[lep+1][2]->Fill(-999);
							EtaDistribution[0][2]->Fill(-999);
							EtaDistribution[lep+1][2]->Fill(-999);
							PhiDistribution[0][2]->Fill(-999);
							PhiDistribution[lep+1][2]->Fill(-999);
							DeltaRDistribution[0][1]->Fill(-999);
							DeltaRDistribution[lep+1][1]->Fill(-999);
							Pt[0][1]->Fill(-999, -999, -999);
							Pt[lep+1][1]->Fill(-999, -999, -999);
							Inclusive[1]->Fill(-999);
							if(FS == 1) Pixel[0][1]->Fill(-999, -999, -999);
							if(FS == 1) Pixel[lep+1][1]->Fill(-999, -999, -999);
						}
					}
										
					if(mass4l_vtxFSR_BS > 105 && mass4l_vtxFSR_BS < 140 && vtxLep_BS_pt->at(lep) > 5 && vtxLep_BS_pt->at(lep) < pT_bins.at(pT_bins.size()-1) && fabs(vtxLep_BS_eta->at(lep)) < 2.4){ 
// 						if(fabs(lep_id->at(lep)) == 13){
						if(13 == 13){
							DR = 999;
							for(int lep2 = 0; lep2 < GENlep_pt->size(); lep2++){
								deltaEta = vtxLep_BS_eta->at(lep) - GENlep_eta->at(lep2);
								deltaPhi = vtxLep_BS_phi->at(lep) - GENlep_phi->at(lep2);
								DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
								DR2 = sqrt(pow(deltaEta, 2) + pow(deltaPhi, 2));
								if(DR2 < DR){
									DR = DR2;
									index = lep2;
								}							
							}
	
							if(DR < minDR){
								int keta, kpt;
								for(int pt = 1; pt < pT_bins.size(); pt++){
									if(vtxLep_BS_pt->at(lep) < pT_bins.at(pt)){
										kpt = pt;
										break;
									}
								}
								for(int eta = 1; eta < eta_bins.size(); eta++){
									if(fabs(vtxLep_BS_eta->at(lep)) < fabs(eta_bins.at(eta))){
										keta = eta;
										break;
									}
								}

								if(m == 0)
									newWeight = 0;
								if(m == 1){
									bin = 2 + (kpt - 1) + (pT_bins.size()-1) * (keta - 1);
									newWeight = FinalMean_GAUSS[2]->GetBinContent((int)bin + 1);
								}
								
								float new_pt = vtxLep_BS_pt->at(lep)*(1-newWeight);
						
								PtDistribution[0][3]->Fill(vtxLep_BS_pt->at(lep)*(1-newWeight));
								PtDistribution[lep+1][3]->Fill(vtxLep_BS_pt->at(lep)*(1-newWeight));
								SigmaPtOverPtDistribution[0][3]->Fill(vtxLep_BS_ptError->at(lep)/vtxLep_BS_pt->at(lep)*(1-newWeight));
								SigmaPtOverPtDistribution[lep+1][3]->Fill(vtxLep_BS_ptError->at(lep)/vtxLep_BS_pt->at(lep)*(1-newWeight));
								EtaDistribution[0][3]->Fill(vtxLep_BS_eta->at(lep));
								EtaDistribution[lep+1][3]->Fill(vtxLep_BS_eta->at(lep));
								PhiDistribution[0][3]->Fill(vtxLep_BS_phi->at(lep));
								PhiDistribution[lep+1][3]->Fill(vtxLep_BS_phi->at(lep));

								DeltaRDistribution[0][2]->Fill(DR);
								DeltaRDistribution[lep+1][2]->Fill(DR);

								DeltaPtDistribution[0][2]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
								DeltaPtDistribution[lep+1][2]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
	
								if(fabs(lep_id->at(lep)) == 13){
									Pt[0][2]->Fill(fabs(vtxLep_BS_eta->at(lep)), new_pt, (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Pt[lep+1][2]->Fill(fabs(vtxLep_BS_eta->at(lep)), new_pt, (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Inclusive[2]->Fill((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
// 								if(FS == 1) Pixel[0][2]->Fill(fabs(vtxLep_BS_eta->at(lep)), lep_numberOfValidPixelHits->at(lep), (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
// 								if(FS == 1) Pixel[lep+1][2]->Fill(fabs(vtxLep_BS_eta->at(lep)), lep_numberOfValidPixelHits->at(lep), (new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));

									rv_ptRecoGen[keta][kpt][2]->setVal((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Data_ptRecoGen[keta][kpt][2]->add(*rastmp_ptRecoGen[keta][kpt][2]);
									rv_ptRecoGen[0][0][2]->setVal((new_pt - GENlep_pt->at(index))/GENlep_pt->at(index));
									Data_ptRecoGen[0][0][2]->add(*rastmp_ptRecoGen[0][0][2]);
								}
								
								VX_BS[0]->Fill(commonBS_x->at(lep), commonBS_y->at(lep));
								VX_BS[lep+1]->Fill(commonBS_x->at(lep), commonBS_y->at(lep));
								VX_BS_Z[0]->Fill(commonBS_z->at(lep));
								VX_BS_Z[lep+1]->Fill(commonBS_z->at(lep));

							}
						}
						else{
							PtDistribution[0][3]->Fill(-999);
							PtDistribution[lep+1][3]->Fill(-999);
							EtaDistribution[0][3]->Fill(-999);
							EtaDistribution[lep+1][3]->Fill(-999);
							PhiDistribution[0][3]->Fill(-999);
							PhiDistribution[lep+1][3]->Fill(-999);
							DeltaRDistribution[0][2]->Fill(-999);
							DeltaRDistribution[lep+1][2]->Fill(-999);
							Pt[0][2]->Fill(-999, -999, -999);
							Pt[lep+1][2]->Fill(-999, -999, -999);
							Inclusive[2]->Fill(-999);
							if(FS == 1) Pixel[0][2]->Fill(-999, -999, -999);
							if(FS == 1) Pixel[lep+1][2]->Fill(-999, -999, -999);
						}
					}

				}					
// 				std::cout<<" ---- "<<std::endl;
		} // for on entry	
		
	}
		
// 		return;

		if(channel != "Full"){
		
		TLegend *legend_reco_gen = new TLegend(0.75,0.75,0.9,0.9);
		TLegend *legend_year = new TLegend(0.75,0.7,0.9,0.9);
		TLegend *legend_width = new TLegend(0.1,0.7,0.25,0.9);
		TLegend *legend_comparison = new TLegend(0.1,0.6,0.2,0.9);

		if(isDR) directory = "./RecoGen_pt_Width_check_Higgs_DeltaR";
// 		else directory = "./RecoGen_pt_Width_check_Higgs_VtxBsFsr";
		else{ 
			if(Base == 1) directory = "./Production_10_2_18_RecoGen_pt_Width_check_Higgs_base";
			else directory = "./Production_10_2_18_RecoGen_pt_Width_check_Higgs_VX_BS";
		}
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		if(FS == 1) directory += "/4mu";
		if(FS == 2) directory += "/4e";
		if(FS == 3) directory += "/2e2mu";
		if(FS == 4) directory += "/2mu2e";
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/" + year;
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		if(channel == "ggF") directory += "/ggF";
		if(channel == "ttH") directory += "/ttH";
		if(channel == "qqZZ") directory += "/qqZZ_pohweg";
		if(channel == "ggZZ") directory += "/ggZZ";
		execute = "mkdir " + directory;
		gSystem->Exec(execute);
		directory += "/";
		
		save_nome = "Higgs_VtxContraint_Studies_" + year;
		if(isDATA)
			save_nome = "Higgs_VtxContraint_Studies_DATA_" + year;	
		TString nome_canvas;
		if(m == 1) save_nome += "_shifted";
			
		Draw(h_blank, h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf[", "", legend_reco_gen, 100);
		for(int i = 0; i < 1; i++){
			
			PtDistribution[0][0]->SetLineColor(kBlue);
			PtDistribution[0][1]->SetLineColor(1);
			PtDistribution[0][2]->SetLineColor(2);
			PtDistribution[0][3]->SetLineColor(3);
			
			legend_year->AddEntry(PtDistribution[0][0], lepton_type.at(0));
			legend_year->AddEntry(PtDistribution[0][1], lepton_type.at(1));
// 			legend_year->AddEntry(PtDistribution[0][2], lepton_type.at(2));
			legend_year->AddEntry(PtDistribution[0][3], lepton_type.at(3));

			legend_reco_gen->AddEntry(PtDistribution[0][1], lepton_type.at(1));
// 			legend_reco_gen->AddEntry(PtDistribution[0][2], lepton_type.at(2));
			legend_reco_gen->AddEntry(PtDistribution[0][3], lepton_type.at(3));

			float fitMin, fitMax;
			if(FS == 1){	fitMin = -0.5;	fitMax = 0.5;}
			if(FS == 3){	fitMin = -0.2;	fitMax = 0.2;}
			if(FS == 4){	fitMin = -0.4;	fitMax = 0.4;}

			nome_canvas = "MassBase vs Mass VX+BS";
			Draw_TH2F(Mass_Base_vs_VXBS, nome_canvas, directory + save_nome + ".pdf", "m_{4l} base [GeV]", "m_{4l} VX+BS [GeV]");
			nome_canvas = "#DeltaMass vs Mass VX+BS";
			Draw_TH2F(DeltaMass_vs_Base, nome_canvas, directory + save_nome + ".pdf", "m_{4l} base [GeV]", "m_{4l} base - m_{4l} VX+BS [GeV]");
			nome_canvas = "#DeltaMass";
			Draw_TH1(rv_DeltaMass, Data_DeltaMass, nome_canvas, directory + save_nome + ".pdf", fitMin, fitMax);
// 			TCanvas *c1 = new TCanvas("#DeltaMass", "#DeltaMass", 700, 700);
// 			DeltaMass->Draw();
// 			DeltaMass->Fit("gaus","","",fitMin,fitMax);
// 			gStyle->SetOptFit(1);
// 			DeltaMass->GetXaxis()->SetTitle("m_{4l} base [GeV] - m_{4l} VX+BS [GeV]");
// 			c1->Print(directory + save_nome + ".pdf");
			nome_canvas = "REFITTED: MassBase vs Mass VX+BS";
			Draw_TH2F(refitted_Mass_Base_vs_VXBS, nome_canvas, directory + save_nome + ".pdf", "m_{4l} base [GeV]", "m_{4l} VX+BS [GeV]");
			nome_canvas = "REFITTED: #DeltaMass vs Mass VX+BS";
			Draw_TH2F(refitted_DeltaMass_vs_Base, nome_canvas, directory + save_nome + ".pdf", "m_{4l} base [GeV]", "m_{4l} base - m_{4l} VX+BS [GeV]");
// 			c1 = new TCanvas("#DeltaRefittedMass", "#DeltaRefittedMass", 700, 700);
// 			refitted_DeltaMass->Draw();
// 			refitted_DeltaMass->GetXaxis()->SetTitle("m_{4l} base [GeV] - m_{4l} VX+BS [GeV]");
// 			if(FS == 1) fitMin = -fitMax = -1;
// 			if(FS == 3) fitMin = -fitMax = -0.25;
// 			if(FS == 4) fitMin = -fitMax = -0.5;			
// 			refitted_DeltaMass->Fit("gaus","","",fitMin,fitMax);
// 			gStyle->SetOptFit(1);
// 			c1->Print(directory + save_nome + ".pdf");
			nome_canvas = "#DeltaRefittedMass";
			Draw_TH1(rv_DeltaRefittedMass, Data_DeltaRefittedMass, nome_canvas, directory + save_nome + ".pdf", fitMin, fitMax);


			for(int lep = 0; lep < 5; lep++){
					if(lep == 0){
						
						Draw(PtDistribution[lep][0], PtDistribution[lep][1], PtDistribution[lep][2],PtDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "p_{T} [GeV]", legend_year, 1000, 1);
						nome_canvas = "SigmaPtOverPtDistribution";
						Draw(SigmaPtOverPtDistribution[lep][1], SigmaPtOverPtDistribution[lep][2],SigmaPtOverPtDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "#sigma_{p_{T}}/p_{T}", legend_reco_gen, 1000, 0);
						nome_canvas = "DeltaPtDistribution";
						Draw(DeltaPtDistribution[0][0], DeltaPtDistribution[0][1], DeltaPtDistribution[0][2], nome_canvas,  directory + save_nome + ".pdf", "(p_{T}^{reco} - p_{T}^{GEN})/p_{T}^{GEN}", legend_reco_gen, 1000, 0);
						nome_canvas = "EtaDistribution";
						Draw(EtaDistribution[lep][0], EtaDistribution[lep][1], EtaDistribution[lep][2],EtaDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "#eta", legend_year, 1000, 1);
						nome_canvas = "PhiDistribution";
						Draw(PhiDistribution[lep][0], PhiDistribution[lep][1], PhiDistribution[lep][2],PhiDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "#phi", legend_year, 1000, 1);
						nome_canvas = "DeltaRDistribution";
						Draw(DeltaRDistribution[lep][0], DeltaRDistribution[lep][1], DeltaRDistribution[lep][2], nome_canvas,  directory + save_nome + ".pdf", "#DeltaR", legend_reco_gen, 1000, 1);
						nome_canvas = "MassDistribution";
						Draw(MassDistribution[1], MassDistribution[2], MassDistribution[3], rv_Mass[0][0], rv_Mass[1][0], rv_Mass[2][0], Data_Mass[0][0], Data_Mass[1][0], Data_Mass[2][0], nome_canvas,  directory + save_nome + ".pdf", "m_{4l} [GeV]", legend_reco_gen, 1, 0);
						nome_canvas = "RefittedMassDistribution";
						Draw(REFITTEDMassDistribution[1], REFITTEDMassDistribution[2], REFITTEDMassDistribution[3], rv_Mass[0][1], rv_Mass[1][1], rv_Mass[2][1], Data_Mass[0][1], Data_Mass[1][1], Data_Mass[2][1], nome_canvas,  directory + save_nome + ".pdf", "m_{4l} [GeV]", legend_reco_gen, 1, 0);
						nome_canvas = "MassErrDistribution";
						Draw(MassErrDistribution[0], MassErrDistribution[1], MassErrDistribution[2], nome_canvas,  directory + save_nome + ".pdf", "m_{4l} [GeV]", legend_reco_gen, 1000, 0);
						nome_canvas = "RefittedMassErrDistribution";
						Draw(REFITTEDMassErrDistribution[0], REFITTEDMassErrDistribution[1], REFITTEDMassErrDistribution[2], nome_canvas,  directory + save_nome + ".pdf", "m_{4l} [GeV]", legend_reco_gen, 1000, 0);
						nome_canvas = "SigmaPt_vsSigmaPt";
						Draw_TH2F(SigmaPt_vsSigmaPt[0], nome_canvas, directory + save_nome + ".pdf", "#sigma_{p_{T}}/p_{T} baseline", "#sigma_{p_{T}}/p_{T} VX+BS");

// 						for(int t = 1; t < 4; t++){
// 							nome_canvas = "DeltaPt_" + lepton_type.at(t);						
// 							Draw_TH3F(Pt[0][t-1], nome_canvas, directory + save_nome + ".pdf", legend_year);
// 							nome_canvas = "DeltaPixel_" + lepton_type.at(t);						
// 							Draw_TH3F(Pixel[0][t-1], nome_canvas, directory + save_nome + ".pdf", legend_year);
// 						}

					
					}
					else{
						nome_canvas = Form("PtDistribution_Lep_%d", lep);
						Draw(PtDistribution[lep][0], PtDistribution[lep][1], PtDistribution[lep][2],PtDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "p_{T} [GeV]", legend_year, 1000, 1);
						nome_canvas = Form("SigmaPtOverPtDistribution_Lep_%d", lep);
						Draw(SigmaPtOverPtDistribution[lep][1], SigmaPtOverPtDistribution[lep][2],SigmaPtOverPtDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "#sigma_{p_{T}}/p_{T}", legend_reco_gen, 1000, 1);
						nome_canvas = Form("DeltaPtDistribution_Lep_%d", lep);
						Draw(DeltaPtDistribution[lep][0], DeltaPtDistribution[lep][1], DeltaPtDistribution[lep][2], nome_canvas,  directory + save_nome + ".pdf", "(p_{T}^{reco} - p_{T}^{GEN})/p_{T}^{GEN}", legend_reco_gen, 1000, 0);
						nome_canvas = Form("EtaDistribution_Lep_%d", lep);
						Draw(EtaDistribution[lep][0], EtaDistribution[lep][1], EtaDistribution[lep][2],EtaDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "#eta", legend_year, 1000, 1);
						nome_canvas = Form("PhiDistribution_Lep_%d", lep);
						Draw(PhiDistribution[lep][0], PhiDistribution[lep][1], PhiDistribution[lep][2],PhiDistribution[lep][3], nome_canvas,  directory + save_nome + ".pdf", "#phi", legend_year, 1000, 1);
						nome_canvas = Form("DeltaRDistribution_Lep_%d", lep);
						Draw(DeltaRDistribution[lep][0], DeltaRDistribution[lep][1], DeltaRDistribution[lep][2], nome_canvas,  directory + save_nome + ".pdf", "#DeltaR", legend_reco_gen, 1000, 1);
						nome_canvas = Form("SigmaPt_vsSigmaPt_Lep_%d", lep);
						Draw_TH2F(SigmaPt_vsSigmaPt[lep], nome_canvas, directory + save_nome + ".pdf", "#sigma_{p_{T}}/p_{T} baseline", "#sigma_{p_{T}}/p_{T} VX+BS");

// 						for(int t = 1; t < 4; t++){
// 							nome_canvas = Form("DeltaPt_Lep_%d_%s", lep, lepton_type[t].Data());						
// 							Draw_TH3F(Pt[lep][t-1], nome_canvas, directory + save_nome + ".pdf", legend_year);
// 							nome_canvas = Form("DeltaPixel_Lep_%d_%s", lep, lepton_type[t].Data());						
// 							Draw_TH3F(Pixel[lep][t-1], nome_canvas, directory + save_nome + ".pdf", legend_year);
// 						}
					}
			}
		}
		Draw(h_blank, h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", "", legend_reco_gen, 100);
	
		
		if(m == 0){
			gStyle->SetOptStat(1);
			save_nome += "_Vertex";
			Draw(h_blank, h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", "", legend_reco_gen, 100);
			TCanvas *cV = new TCanvas("PV", "PV", 700, 500);
			Vertex->Draw("COLZ");
			Vertex->SetStats(1111);
			Vertex->GetXaxis()->SetTitle("x [cm]");
			Vertex->GetYaxis()->SetTitle("y [cm]");
			cV->Print(directory + save_nome + ".pdf[");
			cV->Print(directory + save_nome + ".pdf");
			cV = new TCanvas("BS", "BS", 700, 500);
			BeamSpot->Draw("COLZ");
			BeamSpot->SetStats(1111);
			BeamSpot->GetXaxis()->SetTitle("x [cm]");
			BeamSpot->GetYaxis()->SetTitle("y [cm]");
			cV->Print(directory + save_nome + ".pdf");
			cV = new TCanvas("VX", "VX", 700, 500);
			VX[0]->Draw("COLZ");
			VX[0]->SetStats(1111);
			VX[0]->GetXaxis()->SetTitle("x [cm]");
			VX[0]->GetYaxis()->SetTitle("y [cm]");
			cV->Print(directory + save_nome + ".pdf");
			cV = new TCanvas("VX+BS","VX+BS", 700, 500);
			VX_BS[0]->Draw("COLZ");
			VX_BS[0]->SetStats(1111);
			VX_BS[0]->GetXaxis()->SetTitle("x [cm]");
			VX_BS[0]->GetYaxis()->SetTitle("y [cm]");
			cV->Print(directory + save_nome + ".pdf");
			for(int i = 1; i < 5; i++){
				TString lepton = Form(" %d", i);
				cV = new TCanvas("VX"+lepton, "VX"+lepton, 700, 500);
				VX[i]->Draw("COLZ");
				VX[i]->SetStats(1111);
				VX[i]->GetXaxis()->SetTitle("x [cm]");
				VX[i]->GetYaxis()->SetTitle("y [cm]");
				cV->Print(directory + save_nome + ".pdf");
				cV = new TCanvas("VX+BS"+lepton,"VX+BS"+lepton, 700, 500);
				VX_BS[i]->Draw("COLZ");
				VX_BS[i]->SetStats(1111);
				VX_BS[i]->GetXaxis()->SetTitle("x [cm]");
				VX_BS[i]->GetYaxis()->SetTitle("y [cm]");
				cV->Print(directory + save_nome + ".pdf");
			}
			
			cV = new TCanvas("PV: z", "PV: z", 700, 500);
			Vertex_Z->Draw();
			Vertex_Z->SetStats(1111);
			Vertex_Z->GetXaxis()->SetTitle("z [cm]");
			cV->Print(directory + save_nome + ".pdf");
			cV = new TCanvas("BS: z", "BS: z", 700, 500);
			BeamSpot_Z->Draw();
			BeamSpot_Z->SetStats(1111);
			BeamSpot_Z->GetXaxis()->SetTitle("z [cm]");
			cV->Print(directory + save_nome + ".pdf");
			cV = new TCanvas("VX: z", "VX: z", 700, 500);
			VX_Z[0]->Draw();
			VX_Z[0]->SetStats(1111);
			VX_Z[0]->GetXaxis()->SetTitle("z [cm]");
			cV->Print(directory + save_nome + ".pdf");
			cV = new TCanvas("VX+BS: z","VX+BS: z", 700, 500);
			VX_BS_Z[0]->Draw();
			VX_BS_Z[0]->SetStats(1111);
			VX_BS_Z[0]->GetXaxis()->SetTitle("z [cm]");
			cV->Print(directory + save_nome + ".pdf");
			for(int i = 1; i < 5; i++){
				TString lepton = Form(" %d", i);
				cV = new TCanvas("VX"+lepton+": z", "VX"+lepton+": z", 700, 500);
				VX_Z[i]->Draw();
				VX_Z[i]->SetStats(1111);
				VX_Z[i]->GetXaxis()->SetTitle("z [cm]");
				cV->Print(directory + save_nome + ".pdf");
				cV = new TCanvas("VX+BS"+lepton+": z","VX+BS"+lepton+": z", 700, 500);
				VX_BS_Z[i]->Draw();
				VX_BS_Z[i]->SetStats(1111);
				VX_BS_Z[i]->GetXaxis()->SetTitle("z [cm]");
				cV->Print(directory + save_nome + ".pdf");
				if(i == 5)
					cV->Print(directory + save_nome + ".pdf]");				
			}
			Draw(h_blank, h_blank, h_blank, h_blank, "h_blank", directory + save_nome + ".pdf]", "", legend_reco_gen, 100);
		}		
/*				
		save_nome = "DeltaR_check_" + year;
*/
		
		int t;
		int conuter_entry[5] = {0};
		int inclusive_entry[5] = {0};
		for(int h = 1; h < lepton_type.size(); h++){
			t = h-1;
// 			std::cout<<Pt[0][t]->GetNbinsX()<<"\t"<<Pt[0][t]->GetNbinsY()<<std::endl;
			int bin = 0;
			for(int x = 0; x <= Pt[0][t]->GetNbinsX(); x++){
				for(int y = 0; y <= Pt[0][t]->GetNbinsY(); y++){
					float rangeMin, rangeMax;
					rangeMin = -0.05;
					rangeMax = 0.05;
					if((x == 0 && y == 0) || x == 1){
						rangeMin = -0.03;
						rangeMax = 0.03;
					}

					TString nome;
					TCanvas *c1;

					nome = "Inclusive";
					RooRealVar var("var", "var", genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1));

					RooRealVar Mean("Mean", "Mean", 0, -0.01, 0.01);
					RooRealVar Sigma_Gauss("Sigma_Gauss", "Sigma_Gauss", 0.01, 0.005, 0.1);
					RooRealVar MeanA("MeanA", "MeanA", 0, -0.01, 0.01);
					RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.005, 0.1);

					RooRealVar Sigma_DSCB("Sigma_DSCB", "Sigma_DSCB", 0.01, 0.0005, 0.1);
					RooRealVar AlphaL("AlphaL", "AlphaL", 5, 0, 50);
					RooRealVar ExpL("ExpL", "ExpL", 10, -50, 500);
					RooRealVar AlphaR("AlphaR", "AlphaR", 5, 0, 50);
					RooRealVar ExpR("ExpR", "ExpR", 10, -50, 500);

// 					RooGaussian Gauss("Gauss", "Gauss", var, Mean, Sigma_Gauss);
// 					RooGaussian GaussA("GaussA", "GaussA", var, Mean, Sigma_Gauss);
// 					RooMyPDF_DSCB DSCB("DSCB", "DSCB", var, Mean, Sigma_DSCB, AlphaL, ExpL, AlphaR, ExpR);										
					
					if(x == 0 && y == 0){
						std::cout<<nome<<std::endl;
// 						RooDataHist histo(nome, nome, var, Inclusive[t]);
						RooDataSet histo = RooDataSet(Data_ptRecoGen[0][0][t]->GetName(), Data_ptRecoGen[0][0][t]->GetTitle(), Data_ptRecoGen[0][0][t], *Data_ptRecoGen[0][0][t]->get());//, "1", "weight");

						RooGaussian Gauss("Gauss", "Gauss", *rv_ptRecoGen[x][y][t], Mean, Sigma_Gauss);
						RooGaussian GaussA("GaussA", "GaussA", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussB("GaussB", "GaussB", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussC("GaussC", "GaussC", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussD("GaussD", "GaussD", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussE("GaussE", "GaussE", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussF("GaussF", "GaussF", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
// 						RooMyPDF_DSCB DSCB("DSCB", "DSCB", *rv_ptRecoGen[x][y][t], Mean, Sigma_DSCB, AlphaL, ExpL, AlphaR, ExpR);										
						RooCBShape CB("CB", "CB", *rv_ptRecoGen[x][y][t], Mean, Sigma_DSCB, AlphaL, ExpL);	

						c1 = new TCanvas(nome, nome, 700, 500);
// 						RooPlot* xframe = var.frame(Title(nome));
						RooPlot* xframe = rv_ptRecoGen[0][0][t]->frame(Title(nome));
						histo.plotOn(xframe);

// 						Gauss.fitTo(histo, Range(rangeMin, rangeMax), Save());
// 						Gauss.plotOn(xframe,RooFit::LineColor(kBlack));
// 						Gauss.paramOn(xframe, RooFit::Layout(0.1, 0.45, 0.9));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						FinalSigma_GAUSS[t]->SetBinContent(1, Sigma_Gauss.getVal());
// 						FinalSigma_GAUSS[t]->SetBinError(1, Sigma_Gauss.getError());
// 						FinalMean_GAUSS[t]->SetBinContent(1, Mean.getVal());
// 						FinalMean_GAUSS[t]->SetBinError(1, Mean.getError());

						GaussA.fitTo(histo, Range(genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1)));
// 						GaussA.plotOn(xframe,RooFit::LineColor(kRed+2));
// 						Gauss.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.75));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kRed+2);						
						FinalSigma3_GAUSS[t][0]->SetBinContent(1, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][0]->SetBinError(1, Sigma_GaussA.getError());
						GaussB.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussB.plotOn(xframe,RooFit::LineColor(kGreen+2));
// 						GaussB.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.63));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kGreen+2);	
						FinalSigma3_GAUSS[t][1]->SetBinContent(1, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][1]->SetBinError(1, Sigma_GaussA.getError());
						GaussC.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussC.plotOn(xframe,RooFit::LineColor(kBlue));
// 						GaussC.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.51));						
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kBlue);						
						FinalSigma3_GAUSS[t][2]->SetBinContent(1, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][2]->SetBinError(1, Sigma_GaussA.getError());
						GaussD.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussD.plotOn(xframe,RooFit::LineColor(kYellow+2));
// 						GaussD.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.39));						
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kYellow+2);						
						FinalSigma3_GAUSS[t][3]->SetBinContent(1, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][3]->SetBinError(1, Sigma_GaussA.getError());
						GaussE.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussE.plotOn(xframe,RooFit::LineColor(kBlack));
// 						GaussE.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.27));						
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kBlack);						
						FinalSigma3_GAUSS[t][4]->SetBinContent(1, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][4]->SetBinError(1, Sigma_GaussA.getError());
						GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
						GaussF.plotOn(xframe,RooFit::LineColor(kRed+2));
						GaussF.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
						xframe->getAttText()->SetTextSize(0.025);
						xframe->getAttText()->SetTextColor(kRed+2);						
						FinalSigma3_GAUSS[t][5]->SetBinContent(1, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][5]->SetBinError(1, Sigma_GaussA.getError());
						FinalMean_GAUSS[t]->SetBinContent(1, MeanA.getVal());
						FinalMean_GAUSS[t]->SetBinError(1, MeanA.getError());

// 						DSCB.fitTo(histo, Range(-0.1, 0.1));
// 						DSCB.plotOn(xframe,RooFit::LineColor(kGreen+2));
// 						DSCB.paramOn(xframe, RooFit::Layout(0.6, 0.9, 0.9));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kGreen+2);
// 						FinalSigma_DSCB[t]->SetBinContent(1, Sigma_DSCB.getVal());
// 						FinalSigma_DSCB[t]->SetBinError(1, Sigma_DSCB.getError());
// 						FinalMean_DSCB[t]->SetBinContent(1, Mean.getVal());
// 						FinalMean_DSCB[t]->SetBinError(1, Mean.getError());
	
// 						CB.fitTo(histo, Range(-0.1, 0.1));
// 						CB.plotOn(xframe,RooFit::LineColor(kBlue));
// 						CB.paramOn(xframe, RooFit::Layout(0.6, 0.9, 0.5));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kBlue);						
					
						Double_t integral = Inclusive[t]->Integral();
						inclusive_entry[t] = integral;
						TString integ = Form("Integral = %.0f", integral);
						TLatex* tex1 = new TLatex(0.6,0.3, integ);
						tex1->SetNDC();
					
						xframe->Draw();

// 						tex1->Draw();										

						bin++;
					}
					else if((x == 0 && y != 0) || (x != 0 && y == 0)) continue;
					else{
						nome = Form("Pt_%s_%.0f_%.0f_%s_%s", lepton_type[h].Data(), pT_bins.at(y-1), pT_bins.at(y), eta_bins_name[x-1].Data(), eta_bins_name[x].Data());
						TH1F* new_1 = (TH1F*) Pt[0][t]->ProjectionZ(nome, x, x, y, y);
						new_1->SetTitle(nome);
						new_1->GetXaxis()->SetTitle("(p_{T}^{reco} - p_{T}^{GEN})/p_{T}^{GEN}");

// 						RooDataHist histo(nome, nome, var, new_1);
						RooDataSet histo = RooDataSet(Data_ptRecoGen[x][y][t]->GetName(), Data_ptRecoGen[x][y][t]->GetTitle(), Data_ptRecoGen[x][y][t], *Data_ptRecoGen[x][y][t]->get());//, "1", "weight");

						RooGaussian Gauss("Gauss", "Gauss", *rv_ptRecoGen[x][y][t], Mean, Sigma_Gauss);
						RooGaussian GaussA("GaussA", "GaussA", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussB("GaussB", "GaussB", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussC("GaussC", "GaussC", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussD("GaussD", "GaussD", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussE("GaussE", "GaussE", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
						RooGaussian GaussF("GaussF", "GaussF", *rv_ptRecoGen[x][y][t], MeanA, Sigma_GaussA);
// 						RooMyPDF_DSCB DSCB("DSCB", "DSCB", *rv_ptRecoGen[x][y][t], Mean, Sigma_DSCB, AlphaL, ExpL, AlphaR, ExpR);										
						RooCBShape CB("CB", "CB", *rv_ptRecoGen[x][y][t], Mean, Sigma_DSCB, AlphaL, ExpL);	

						c1 = new TCanvas(nome, nome, 700, 500);
// 						RooPlot* xframe = var.frame(Title(nome));
						RooPlot* xframe = rv_ptRecoGen[x][y][t]->frame(Title(nome));
						histo.plotOn(xframe);

// 						Gauss.fitTo(histo, Range(rangeMin, rangeMax));
// 						Gauss.plotOn(xframe,RooFit::LineColor(kBlack));
// 						Gauss.paramOn(xframe, RooFit::Layout(0.1, 0.45, 0.9));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						FinalSigma_GAUSS[t]->SetBinContent(bin+2, Sigma_Gauss.getVal());
// 						FinalSigma_GAUSS[t]->SetBinError(bin+2, Sigma_Gauss.getError());
// 						FinalMean_GAUSS[t]->SetBinContent(bin+2, Mean.getVal());
// 						FinalMean_GAUSS[t]->SetBinError(bin+2, Mean.getError());

						GaussA.fitTo(histo, Range(genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1)));
// 						GaussA.plotOn(xframe,RooFit::LineColor(kRed+2));
// 						Gauss.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.75));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kRed+2);						
						FinalSigma3_GAUSS[t][0]->SetBinContent(bin+2, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][0]->SetBinError(bin+2, Sigma_GaussA.getError());
						GaussB.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussB.plotOn(xframe,RooFit::LineColor(kGreen+2));
// 						GaussB.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.63));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kGreen+2);	
						FinalSigma3_GAUSS[t][1]->SetBinContent(bin+2, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][1]->SetBinError(bin+2, Sigma_GaussA.getError());
						GaussC.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussC.plotOn(xframe,RooFit::LineColor(kBlue));
// 						GaussC.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.51));						
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kBlue);						
						FinalSigma3_GAUSS[t][2]->SetBinContent(bin+2, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][2]->SetBinError(bin+2, Sigma_GaussA.getError());
						GaussD.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussD.plotOn(xframe,RooFit::LineColor(kYellow+2));
// 						GaussD.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.39));						
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kYellow+2);						
						FinalSigma3_GAUSS[t][3]->SetBinContent(bin+2, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][3]->SetBinError(bin+2, Sigma_GaussA.getError());
						GaussE.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 						GaussE.plotOn(xframe,RooFit::LineColor(kBlack));
// 						GaussE.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.27));						
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kBlack);						
						FinalSigma3_GAUSS[t][4]->SetBinContent(bin+2, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][4]->SetBinError(bin+2, Sigma_GaussA.getError());
						GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
						GaussF.plotOn(xframe,RooFit::LineColor(kRed+2));
						GaussF.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
						xframe->getAttText()->SetTextSize(0.025);
						xframe->getAttText()->SetTextColor(kRed+2);						
						FinalSigma3_GAUSS[t][5]->SetBinContent(bin+2, Sigma_GaussA.getVal());
						FinalSigma3_GAUSS[t][5]->SetBinError(bin+2, Sigma_GaussA.getError());
						FinalMean_GAUSS[t]->SetBinContent(bin+2, MeanA.getVal());
						FinalMean_GAUSS[t]->SetBinError(bin+2, MeanA.getError());

// 						DSCB.fitTo(histo, Range(-0.1, 0.1));
// 						DSCB.plotOn(xframe,RooFit::LineColor(kGreen+2));
// 						DSCB.paramOn(xframe, RooFit::Layout(0.6, 0.9, 0.9));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kGreen+2);						
// 						FinalSigma_DSCB[t]->SetBinContent(bin+2, Sigma_DSCB.getVal());
// 						FinalSigma_DSCB[t]->SetBinError(bin+2, Sigma_DSCB.getError());
// 						FinalMean_DSCB[t]->SetBinContent(bin+2, Mean.getVal());
// 						FinalMean_DSCB[t]->SetBinError(bin+2, Mean.getError());

// 						CB.fitTo(histo, Range(-0.1, 0.1));
// 						CB.plotOn(xframe,RooFit::LineColor(kBlue));
// 						CB.paramOn(xframe, RooFit::Layout(0.6, 0.9, 0.5));
// 						xframe->getAttText()->SetTextSize(0.025);
// 						xframe->getAttText()->SetTextColor(kBlue);
					
						Double_t integral = new_1->Integral();
						conuter_entry[t] += integral;
						TString integ = Form("Integral = %.0f", integral);
						TLatex* tex1 = new TLatex(0.6,0.3, integ);
						tex1->SetNDC();

						xframe->Draw();

// 						tex1->Draw();										
											
						bin++;					
					}
					
					if(m == 0){
						if((x == 0 && y == 0) || (x != 0 && y != 0)){	
							if(x == 0 && y == 0)
								c1->Print(directory + lepton_type.at(h) + ".pdf[");
							c1->Print(directory + lepton_type.at(h) + ".pdf");
							if(x == Pt[0][t]->GetNbinsX() && y == Pt[0][t]->GetNbinsY())
								c1->Print(directory + lepton_type.at(h) + ".pdf]");
						}
					}
					if(m == 1){
						if((x == 0 && y == 0) || (x != 0 && y != 0)){	
							if(x == 0 && y == 0)
								c1->Print(directory + lepton_type.at(h) + "_shifted.pdf[");
							c1->Print(directory + lepton_type.at(h) + "_shifted.pdf");
							if(x == Pt[0][t]->GetNbinsX() && y == Pt[0][t]->GetNbinsY())
								c1->Print(directory + lepton_type.at(h) + "_shifted.pdf]");
						}
					}
					
				}
			}
		}
					
		for(int h = 1; h < lepton_type.size(); h++){
			t = h-1;

			TString nome;
			TCanvas *c1;

			nome = "Inclusive";
			RooRealVar var("var", "var", genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1));

			RooRealVar Mean("Mean", "Mean", 0, -0.01, 0.01);
			RooRealVar Sigma_Gauss("Sigma_Gauss", "Sigma_Gauss", 0.01, 0.005, 0.1);
			RooRealVar MeanA("MeanA", "MeanA", 0, -0.01, 0.01);
			RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.005, 0.1);

			std::cout<<nome<<std::endl;
			RooDataSet histo = RooDataSet(Data_DeltaMassOverMass[t]->GetName(), Data_DeltaMassOverMass[t]->GetTitle(), Data_DeltaMassOverMass[t], *Data_DeltaMassOverMass[t]->get());//, "1", "weight");
			RooGaussian Gauss("Gauss", "Gauss", *rv_DeltaMassOverMass[t], Mean, Sigma_Gauss);
			RooGaussian GaussA("GaussA", "GaussA", *rv_DeltaMassOverMass[t], MeanA, Sigma_GaussA);
			RooGaussian GaussB("GaussB", "GaussB", *rv_DeltaMassOverMass[t], MeanA, Sigma_GaussA);
			RooGaussian GaussC("GaussC", "GaussC", *rv_DeltaMassOverMass[t], MeanA, Sigma_GaussA);
			RooGaussian GaussD("GaussD", "GaussD", *rv_DeltaMassOverMass[t], MeanA, Sigma_GaussA);
			RooGaussian GaussE("GaussE", "GaussE", *rv_DeltaMassOverMass[t], MeanA, Sigma_GaussA);
			RooGaussian GaussF("GaussF", "GaussF", *rv_DeltaMassOverMass[t], MeanA, Sigma_GaussA);

			c1 = new TCanvas(nome, nome, 700, 500);
			RooPlot* xframe = rv_DeltaMassOverMass[t]->frame(Title(nome));
			histo.plotOn(xframe);
			GaussA.fitTo(histo, Range(genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1)));
// 			FinalSigma3_GAUSS[t][0]->SetBinContent(1, Sigma_GaussA.getVal());
// 			FinalSigma3_GAUSS[t][0]->SetBinError(1, Sigma_GaussA.getError());
			GaussB.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 			FinalSigma3_GAUSS[t][1]->SetBinContent(1, Sigma_GaussA.getVal());
// 			FinalSigma3_GAUSS[t][1]->SetBinError(1, Sigma_GaussA.getError());
			GaussC.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 			FinalSigma3_GAUSS[t][2]->SetBinContent(1, Sigma_GaussA.getVal());
// 			FinalSigma3_GAUSS[t][2]->SetBinError(1, Sigma_GaussA.getError());
			GaussD.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 			FinalSigma3_GAUSS[t][3]->SetBinContent(1, Sigma_GaussA.getVal());
// 			FinalSigma3_GAUSS[t][3]->SetBinError(1, Sigma_GaussA.getError());
			GaussE.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
// 			FinalSigma3_GAUSS[t][4]->SetBinContent(1, Sigma_GaussA.getVal());
// 			FinalSigma3_GAUSS[t][4]->SetBinError(1, Sigma_GaussA.getError());
			GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
			GaussF.plotOn(xframe,RooFit::LineColor(kRed+2));
			GaussF.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
			xframe->getAttText()->SetTextSize(0.025);
			xframe->getAttText()->SetTextColor(kRed+2);						
// 			FinalSigma3_GAUSS[t][5]->SetBinContent(1, Sigma_GaussA.getVal());
// 			FinalSigma3_GAUSS[t][5]->SetBinError(1, Sigma_GaussA.getError());
// 			FinalMean_GAUSS[t]->SetBinContent(1, MeanA.getVal());
// 			FinalMean_GAUSS[t]->SetBinError(1, MeanA.getError());

			Double_t integral = Inclusive[t]->Integral();
			inclusive_entry[t] = integral;
			TString integ = Form("Integral = %.0f", integral);
			TLatex* tex1 = new TLatex(0.6,0.3, integ);
			tex1->SetNDC();
					
			xframe->Draw();

// 			tex1->Draw();										
					
			if(m == 0){
// 				if(h == 1)
// 					c1->Print(directory + lepton_type.at(h) + "DeltaMOverM.pdf[");
				c1->Print(directory + lepton_type.at(h) + "DeltaMOverM.pdf");
// 				if(h == lepton_type.size()-1)
// 					c1->Print(directory + lepton_type.at(h) + "DeltaMOverM.pdf]");
			}
			if(m == 1){
// 				if(h == 1)
// 						c1->Print(directory + lepton_type.at(h) + "DeltaMOverM_shifted.pdf[");
				c1->Print(directory + lepton_type.at(h) + "DeltaMOverM_shifted.pdf");
// 				if(h == lepton_type.size()-1)
// 					c1->Print(directory + lepton_type.at(h) + "DeltaMOverM_shifted.pdf]");
			}
		}

					
										
		int bin = 0;
		for(int x = 0; x <= Pt[0][t]->GetNbinsX(); x++){
				for(int y = 0; y <= Pt[0][t]->GetNbinsY(); y++){
				
					float rangeMin, rangeMax;
					rangeMin = -0.05;
					rangeMax = 0.05;
					if((x == 0 && y == 0) || x == 1){
						rangeMin = -0.03;
						rangeMax = 0.03;
					}

					TString nome;
					TCanvas *c1;

					nome = "Inclusive";
					RooRealVar var("var", "var", genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1));

					RooRealVar Mean("Mean", "Mean", 0, -0.01, 0.01);
					RooRealVar Sigma_Gauss("Sigma_Gauss", "Sigma_Gauss", 0.01, 0.005, 0.1);
					RooRealVar MeanA("MeanA", "MeanA", 0, -0.01, 0.01);
					RooRealVar Sigma_GaussA("Sigma_GaussA", "Sigma_GaussA", 0.01, 0.005, 0.1);

					RooRealVar Sigma_DSCB("Sigma_DSCB", "Sigma_DSCB", 0.01, 0.0005, 0.1);
					RooRealVar AlphaL("AlphaL", "AlphaL", 5, 0, 50);
					RooRealVar ExpL("ExpL", "ExpL", 10, -50, 500);
					RooRealVar AlphaR("AlphaR", "AlphaR", 5, 0, 50);
					RooRealVar ExpR("ExpR", "ExpR", 10, -50, 500);

						if(x == 0 && y == 0){
							std::cout<<nome<<std::endl;
// 							RooDataHist histo(nome, nome, var, Inclusive[t]);
							RooDataSet histo = RooDataSet(Data_shift[0][0]->GetName(), Data_shift[0][0]->GetTitle(), Data_shift[0][0], *Data_shift[0][0]->get());//, "1", "weight");

							RooGaussian Gauss("Gauss", "Gauss", *rv_shift[x][y], Mean, Sigma_Gauss);
							RooGaussian GaussA("GaussA", "GaussA", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussB("GaussB", "GaussB", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussC("GaussC", "GaussC", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussD("GaussD", "GaussD", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussE("GaussE", "GaussE", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussF("GaussF", "GaussF", *rv_shift[x][y], MeanA, Sigma_GaussA);

							c1 = new TCanvas(nome, nome, 700, 500);
							RooPlot* xframe = rv_shift[0][0]->frame(Title(nome));
							histo.plotOn(xframe);

							GaussA.fitTo(histo, Range(genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1)));
							GaussB.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussC.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussD.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussE.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussF.plotOn(xframe,RooFit::LineColor(kRed+2));
							GaussF.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
							xframe->getAttText()->SetTextSize(0.025);
							xframe->getAttText()->SetTextColor(kRed+2);						
							FinalShift->SetBinContent(1, MeanA.getVal());
							FinalShift->SetBinError(1, MeanA.getError());
					
							Double_t integral = Inclusive[t]->Integral();
							inclusive_entry[t] = integral;
							TString integ = Form("Integral = %.0f", integral);
							TLatex* tex1 = new TLatex(0.6,0.3, integ);
							tex1->SetNDC();
						
							xframe->Draw();

//	 						tex1->Draw();										
	
							bin++;
						}
						else if((x == 0 && y != 0) || (x != 0 && y == 0)) continue;
						else{
							nome = Form("Pt_%.0f_%.0f_%s_%s", pT_bins.at(y-1), pT_bins.at(y), eta_bins_name[x-1].Data(), eta_bins_name[x].Data());
							TH1F* new_1 = (TH1F*) Pt[0][t]->ProjectionZ(nome, x, x, y, y);
							new_1->SetTitle(nome);
							new_1->GetXaxis()->SetTitle("(p_{T}^{baseline} - p_{T}^{VX+BS})/p_{T}^{baseline}");

							RooDataSet histo = RooDataSet(Data_shift[x][y]->GetName(), Data_shift[x][y]->GetTitle(), Data_shift[x][y], *Data_shift[x][y]->get());//, "1", "weight");

							RooGaussian Gauss("Gauss", "Gauss", *rv_shift[x][y], Mean, Sigma_Gauss);
							RooGaussian GaussA("GaussA", "GaussA", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussB("GaussB", "GaussB", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussC("GaussC", "GaussC", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussD("GaussD", "GaussD", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussE("GaussE", "GaussE", *rv_shift[x][y], MeanA, Sigma_GaussA);
							RooGaussian GaussF("GaussF", "GaussF", *rv_shift[x][y], MeanA, Sigma_GaussA);

							c1 = new TCanvas(nome, nome, 700, 500);
							RooPlot* xframe = rv_shift[x][y]->frame(Title(nome));
							histo.plotOn(xframe);
	
// 							FinalMean_GAUSS[t]->SetBinContent(bin+2, Mean.getVal());
// 							FinalMean_GAUSS[t]->SetBinError(bin+2, Mean.getError());
// 
							GaussA.fitTo(histo, Range(genReco_bins.at(0), genReco_bins.at(genReco_bins.size()-1)));
							GaussB.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussC.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussD.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussE.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussF.fitTo(histo, Range(-2*Sigma_GaussA.getVal(), 2*Sigma_GaussA.getVal()));
							GaussF.plotOn(xframe,RooFit::LineColor(kRed+2));
							GaussF.paramOn(xframe, RooFit::Layout(0.1, 0.4, 0.7));						
							xframe->getAttText()->SetTextSize(0.025);
							xframe->getAttText()->SetTextColor(kRed+2);						
							FinalShift->SetBinContent(bin+2, MeanA.getVal());
							FinalShift->SetBinError(bin+2, MeanA.getError());
					
							Double_t integral = new_1->Integral();
							conuter_entry[t] += integral;
							TString integ = Form("Integral = %.0f", integral);
							TLatex* tex1 = new TLatex(0.6,0.3, integ);
							tex1->SetNDC();
	
							xframe->Draw();
	
// 							tex1->Draw();										
												
							bin++;					
						}
					
						if(m == 0){
							if((x == 0 && y == 0) || (x != 0 && y != 0)){	
								if(x == 0 && y == 0)
									c1->Print(directory + "Shift.pdf[");
								c1->Print(directory + "Shift.pdf");
								if(x == Pt[0][t]->GetNbinsX() && y == Pt[0][t]->GetNbinsY())
									c1->Print(directory + "Shift.pdf]");
							}
						}
						if(m == 1){
							if((x == 0 && y == 0) || (x != 0 && y != 0)){	
								if(x == 0 && y == 0)
									c1->Print(directory + "Shift_shifted.pdf[");
								c1->Print(directory + "Shift_shifted.pdf");
								if(x == Pt[0][t]->GetNbinsX() && y == Pt[0][t]->GetNbinsY())
									c1->Print(directory + "Shift_shifted.pdf]");
							}
						}
				}
			}

		for(int t = 0; t < lepton_type.size(); t++)
			std::cout<<inclusive_entry[t]<<"\t"<<conuter_entry[t]<<std::endl;
				
		TCanvas* canvas = new TCanvas("Gauss Width", "Gauss Width", 750, 500);
	   	canvas->cd();
	   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
	   	pad11->SetGrid();
	   	pad11->SetBottomMargin(0.1);
	   	pad11->Draw();
	   	pad11->cd();
	   	pad11->SetTicks();
		
		canvas->Update();

		for(int h = 1; h < lepton_type.size(); h++){
			t = h-1;
// 			FinalSigma_GAUSS[t]->Sumw2();
			FinalSigma_GAUSS[t]->SetLineColor(h);
			FinalSigma_GAUSS[t]->SetMarkerColor(h);
			FinalSigma_GAUSS[t]->SetMarkerStyle(20);
			FinalSigma_GAUSS[t]->GetYaxis()->SetRangeUser(0.005, 0.04);
			FinalSigma_GAUSS[t]->GetYaxis()->SetTitle("Gauss #sigma");
			FinalSigma_GAUSS[t]->GetXaxis()->SetTitle("p_{T} and #eta bins");
			FinalSigma_GAUSS[t]->SetTitle("Gauss #sigma");

			if(t == 0)
				FinalSigma_GAUSS[t]->Draw();
			else{
				FinalSigma_GAUSS[t]->Draw("same");
			}

			legend_width->AddEntry(FinalSigma_GAUSS[t], lepton_type.at(h));

		}
		legend_width->Draw();

		canvas->Update();		
		canvas->cd();
		
		TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
		pad22->SetGrid();
		pad22->Draw();
		pad22->cd();
		pad22->SetTicks();
		
		for(int h = 2; h < lepton_type.size(); h++){
			t = h-1;
			TH1F* ratio = (TH1F*) FinalSigma_GAUSS[t]->Clone();
			ratio->Divide(FinalSigma_GAUSS[0]);
			ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
			ratio->SetTitle("");
			ratio->SetStats(0);
			ratio->GetYaxis()->SetTitleSize(0.1);
			ratio->GetYaxis()->SetLabelSize(0.14);
			ratio->GetYaxis()->SetTitle("X / Baseline");
			ratio->GetYaxis()->SetTitleOffset(0.40);
			ratio->GetYaxis()->SetNdivisions(506); 
			ratio->SetLineColor(h);
			ratio->SetMarkerColor(h);
			if(t == 1)
				ratio->Draw();
			else
				ratio->Draw("same");
			legend_comparison->AddEntry(ratio, lepton_type.at(h));
		}
		legend_comparison->Draw();
	
		canvas->Update();
		pad22->Update();

		TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
		line->SetLineColor(kBlack);
		line->SetLineWidth(1);
		line->Draw();
	
		if(m == 0)
			canvas->Print(directory + "GaussWidth.pdf");
		if(m == 1)
			canvas->Print(directory + "GaussWidth_shifted.pdf");
			
		///// VX+BS shift wrt base /////
		canvas = new TCanvas("VX+BS shift wrt baseline", "VX+BS shift wrt baseline", 750, 500);
	   	canvas->cd();

		FinalShift->SetMarkerStyle(20);
		FinalShift->GetYaxis()->SetRangeUser(-0.005, 0.005);
		FinalShift->GetYaxis()->SetTitle("Gauss #sigma");
		FinalShift->GetXaxis()->SetTitle("p_{T} and #eta bins");
		FinalShift->SetTitle("Gauss #sigma");
		FinalShift->Draw();
	
		if(m == 0)
			canvas->Print(directory + "VX_BS_shift.pdf");
		if(m == 1)
			canvas->Print(directory + "VX_BS_shift_shifted.pdf");
		///// VX+BS shift wrt base /////		
		
		//////// recoursive guassian /////
		for(int k = 0; k < 6; k++){
			canvas = new TCanvas("Gauss Width", "Gauss Width", 750, 500);
		   	canvas->cd();
	   		pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
		   	pad11->SetGrid();
		   	pad11->SetBottomMargin(0.1);
	   		pad11->Draw();
		   	pad11->cd();
		   	pad11->SetTicks();
		
			canvas->Update();
			for(int h = 1; h < lepton_type.size(); h++){
				t = h-1;
				FinalSigma3_GAUSS[t][k]->SetLineColor(h);			
				FinalSigma3_GAUSS[t][k]->SetMarkerColor(h);
				FinalSigma3_GAUSS[t][k]->SetMarkerStyle(20);
				FinalSigma3_GAUSS[t][k]->GetYaxis()->SetRangeUser(0.005, 0.04);
				FinalSigma3_GAUSS[t][k]->GetYaxis()->SetTitle("Gauss #sigma");
				FinalSigma3_GAUSS[t][k]->GetXaxis()->SetTitle("p_{T} and #eta bins");
				FinalSigma3_GAUSS[t][k]->SetTitle("Gauss #sigma");

				if(t == 0)
					FinalSigma3_GAUSS[t][k]->Draw();
				else{
					FinalSigma3_GAUSS[t][k]->Draw("same");
				}

			}
			legend_width->Draw();

			canvas->Update();		
			canvas->cd();
		
			pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
			pad22->SetGrid();
			pad22->Draw();
			pad22->cd();
			pad22->SetTicks();
		
			for(int h = 2; h < lepton_type.size(); h++){
				t = h-1;
				TH1F* ratio = (TH1F*) FinalSigma3_GAUSS[t][k]->Clone();
				ratio->Divide(FinalSigma3_GAUSS[0][k]);
				ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
				ratio->SetTitle("");
				ratio->SetStats(0);
				ratio->GetYaxis()->SetTitleSize(0.1);
				ratio->GetYaxis()->SetLabelSize(0.14);
				ratio->GetYaxis()->SetTitle("X / Baseline");
				ratio->GetYaxis()->SetTitleOffset(0.40);
				ratio->GetYaxis()->SetNdivisions(506); 
				ratio->SetLineColor(h);
				ratio->SetMarkerColor(h);
				if(t == 1)
					ratio->Draw();
				else
					ratio->Draw("same");
			}
			legend_comparison->Draw();
	
			canvas->Update();
			pad22->Update();

			line->Draw();

			if(m == 0){
				if(k == 0)
					canvas->Print(directory + "GaussWidthRecoursive.pdf[");
				canvas->Print(directory + "GaussWidthRecoursive.pdf");
				if(k == 5)
					canvas->Print(directory + "GaussWidthRecoursive.pdf]");				
			}
			if(m == 1){
				if(k == 0)
					canvas->Print(directory + "GaussWidthRecoursive_shifted.pdf[");
				canvas->Print(directory + "GaussWidthRecoursive_shifted.pdf");
				if(k == 5)
					canvas->Print(directory + "GaussWidthRecoursive_shifted.pdf]");				
			}
		}
		//////// ricorsive guassian /////

//  		Draw(bool dscb, bool sigma, int lepton_size, TH1F* FinalSigma_DSCB[3], TLegend* legend_width, TLegend* legend_comparison, TString save_name, int m);
//  		Draw(1, 1, lepton_type.size(), FinalSigma_DSCB, legend_width, legend_comparison, directory, m, line);
//  		Draw(0, 0, lepton_type.size(), FinalMean_GAUSS, legend_width, legend_comparison, directory, m, line);
//  		Draw(1, 0, lepton_type.size(), FinalMean_DSCB, legend_width, legend_comparison, directory, m, line);


		canvas = new TCanvas("DSCB Width", "DSCB Width", 750, 500);
	   	canvas->cd();
	   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
	   	pad11->SetGrid();
	   	pad11->SetBottomMargin(0.1);
	   	pad11->Draw();
	   	pad11->cd();
	   	pad11->SetTicks();
		
		canvas->Update();


		//////// mean guassian /////
		canvas = new TCanvas("Gauss mean", "Gauss mean", 750, 500);
	   	canvas->cd();
	   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
	   	pad11->SetGrid();
	   	pad11->SetBottomMargin(0.1);
	   	pad11->Draw();
	   	pad11->cd();
	   	pad11->SetTicks();
		
		canvas->Update();

		for(int h = 1; h < lepton_type.size(); h++){
			t = h-1;
			FinalMean_GAUSS[t]->SetLineColor(h);
			FinalMean_GAUSS[t]->SetMarkerColor(h);
			FinalMean_GAUSS[t]->SetMarkerColor(h);
			FinalMean_GAUSS[t]->SetMarkerStyle(20);
			FinalMean_GAUSS[t]->GetYaxis()->SetRangeUser(-0.01, 0.01);
			FinalMean_GAUSS[t]->GetYaxis()->SetTitle("Gauss mean");
			FinalMean_GAUSS[t]->GetXaxis()->SetTitle("p_{T} and #eta bins");
			FinalMean_GAUSS[t]->SetTitle("Gauss mean");

			if(t == 0)
				FinalMean_GAUSS[t]->Draw();
			else{
				FinalMean_GAUSS[t]->Draw("same");
			}
		}
// 		legend_year->Draw();
// 		legend_reco_gen->Draw();
		legend_width->Draw();

		canvas->Update();		
		canvas->cd();
		
		pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
		pad22->SetGrid();
		pad22->Draw();
		pad22->cd();
		pad22->SetTicks();
		
		for(int h = 2; h < lepton_type.size(); h++){
			t = h-1;
			TH1F* ratio = (TH1F*) FinalMean_GAUSS[t]->Clone();
			ratio->Add(FinalMean_GAUSS[0], -1);
// 			ratio->Divide(Fin.alMean_GAUSS[0]);
			ratio->GetYaxis()->SetRangeUser(-0.01, 0.01);
			ratio->SetTitle("");
			ratio->SetStats(0);
			ratio->GetYaxis()->SetTitleSize(0.1);
			ratio->GetYaxis()->SetLabelSize(0.14);
			ratio->GetYaxis()->SetTitle("X - Base");
			ratio->GetYaxis()->SetTitleOffset(0.40);
			ratio->GetYaxis()->SetNdivisions(506); 
			ratio->SetLineColor(h);
			ratio->SetMarkerColor(h);
			if(t == 1)
				ratio->Draw();
			else
				ratio->Draw("same");
		}
		legend_comparison->Draw();
	
		canvas->Update();
		pad22->Update();

		line = new TLine(pad22->GetUxmin(), 0, pad22->GetUxmax(), 0);
		line->Draw();
		
		if(m == 0)
			canvas->Print(directory + "GaussMean.pdf");
		if(m == 1)
			canvas->Print(directory + "GaussMean_shifted.pdf");	
		


		canvas = new TCanvas("DSCB mean", "DSCB mean", 750, 500);
	   	canvas->cd();
	   	pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
	   	pad11->SetGrid();
	   	pad11->SetBottomMargin(0.1);
	   	pad11->Draw();
	   	pad11->cd();
	   	pad11->SetTicks();
		
		canvas->Update();


		if(m == 0){
			for(int i = 1; i < 4; i++){
				for(int bin = 1; bin < W_MassDistribution[i]->GetNbinsX(); bin++){
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 70)
						integral_70[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 70 && W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 105)
						integral_70_105[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 105 && W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 140)
						integral_105_140[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 105 && W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 130)
						integral_105_130[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 120 && W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 130)
						integral_120_130[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 118 && W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 130)
						integral_118_130[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);
					if(W_MassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 140)
						integral_140_inf[i-1][0] += W_MassDistribution[i]->GetBinContent(bin);					
				}
			}

			for(int i = 1; i < 4; i++){
				for(int bin = 1; bin < W_REFITTEDMassDistribution[i]->GetNbinsX(); bin++){
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 70)
						integral_70[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 70 && W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 105)
						integral_70_105[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 105 && W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 140)
						integral_105_140[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 105 && W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 130)
						integral_105_130[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 120 && W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 130)
						integral_120_130[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 118 && W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) < 130)
						integral_118_130[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);
					if(W_REFITTEDMassDistribution[i]->GetXaxis()->GetBinLowEdge(bin) >= 140)
						integral_140_inf[i-1][1] += W_REFITTEDMassDistribution[i]->GetBinContent(bin);					
				}
			}
		}


	} // for on m

			std::cout<<"Mass\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t"<<W_MassDistribution[i]->Integral()<<"\t\t"<<W_REFITTEDMassDistribution[i]->Integral()<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [70, INF] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_70[i-1][0]<<"\t\t"<<integral_70[i-1][1]<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [70, 105] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_70_105[i-1][0]<<"\t\t"<<integral_70_105[i-1][1]<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [105, 140] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_105_140[i-1][0]<<"\t\t"<<integral_105_140[i-1][1]<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [105, 130] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_105_130[i-1][0]<<"\t\t"<<integral_105_130[i-1][1]<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [120, 130] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_120_130[i-1][0]<<"\t\t"<<integral_120_130[i-1][1]<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [118, 130] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_118_130[i-1][0]<<"\t\t"<<integral_118_130[i-1][1]<<std::endl;
			std::cout<<" ----------------------------------- "<<std::endl;		
			std::cout<<"Mass --> [140, INF] GeV\t"<<"MASS\t"<<"REFITTED\t"<<std::endl;
			for(int i = 1; i < 4; i++)
				std::cout<<lepton_type.at(i)<<"\t\t\t\t"<<integral_140_inf[i-1][0]<<"\t\t"<<integral_140_inf[i-1][1]<<std::endl;
	
}


void Draw_TH3F(TH3F* h1, TString nome_canvas, TString save, TLegend *legend){

// 	h1->Scale(1/h1->Integral());
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
	h1->Draw("BOX");
// 	h1->GetZaxis()->SetRangeUser(0, 2);
// 	legend_year->Draw();

	canvas->Print(save);
}

void Draw_TH2F(TH2F* h1, TString nome_canvas, TString save, TString x_axis, TString y_axis){

// 	h1->Scale(1/h1->Integral());
	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0, 1, 1);
   	pad11->SetGrid();
   	pad11->SetLeftMargin(0.15);
   	pad11->Draw();
   	pad11->cd();

	h1->Draw();
	h1->GetXaxis()->SetTitle(x_axis);
	h1->GetYaxis()->SetTitle(y_axis);
// 	h1->GetZaxis()->SetRangeUser(0, 2);
// 	legend_year->Draw();

	canvas->Print(save);
}

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TH1F *h4, TString nome_canvas, TString save, TString x_name, TLegend *legend, int fit, bool LogY){
	
	std::vector<float> result_baseline;
	std::vector<float> result_CommonVertex;
	std::vector<float> result_CommonVertex_FSR;
	std::vector<float> result_CommonVertex_BS;
	std::vector<float> result_CommonVertex_BS_FSR;

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
   		min = 0.000001;
   	}
	else
		min = 0;
	
// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());
// 	h4->Scale(1/h4->Integral());

	float tmp_max;
	tmp_max = Max(h1->GetMaximum(), h2->GetMaximum());
	max = Max(h3->GetMaximum(), tmp_max);
	tmp_max = Max(h4->GetMaximum(), max);
		
	max = 2 * max;
	
	if(max == 0) max = 1;	
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h3->SetTitle(nome_canvas);
	h4->SetTitle(nome_canvas);
	h1->SetLineColor(kBlue);
	h2->SetLineColor(1);
	h3->SetLineColor(2);
	h4->SetLineColor(3);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h3->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h4->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	h1->GetXaxis()->SetTitle(x_name);
	h2->GetXaxis()->SetTitle(x_name);
	h3->GetXaxis()->SetTitle(x_name);
	h4->GetXaxis()->SetTitle(x_name);
	
	canvas->Update();
	
	h1->Draw("E");
	h2->Draw("Same E");
// 	h3->Draw("Same E");
	h4->Draw("Same E");

	legend->Draw();	

	if(fit == 0){

		TLatex* tex1 = new TLatex(0.6,0.67, "RMS:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", h1->GetRMS());
		tex1 = new TLatex(0.6,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetRMS());
		tex1 = new TLatex(0.6,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetRMS());
		tex1 = new TLatex(0.6,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kRed+2);
		tex1->Draw();
		
		ciao = Form("Vertex: BS = %.3f", h4->GetRMS());
		tex1 = new TLatex(0.6,0.47, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kBlue);
		tex1->Draw();

		tex1 = new TLatex(0.15,0.67, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		ciao = Form("FSR = %.3f", h1->GetMean());
		tex1 = new TLatex(0.15,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetMean());
		tex1 = new TLatex(0.15,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetMean());
		tex1 = new TLatex(0.15,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kRed+2);
		tex1->Draw();
		
		ciao = Form("Vertex: BS = %.3f", h4->GetMean());
		tex1 = new TLatex(0.15,0.47, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kBlue);
		tex1->Draw();
	}
	if(fit == 1){

		TLatex* tex1 = new TLatex(0.6,0.47, "#sigma:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", result_baseline.at(1));
		tex1 = new TLatex(0.6,0.42, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", result_CommonVertex.at(1));
		tex1 = new TLatex(0.6,0.37, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", result_CommonVertex_FSR.at(1));
		tex1 = new TLatex(0.6,0.32, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kRed+2);
		tex1->Draw();
		
		ciao = Form("Vertex: BS = %.3f", result_CommonVertex_BS.at(1));
		tex1 = new TLatex(0.6,0.27, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kBlue);
		tex1->Draw();

		ciao = Form("Vertex FSR: BS = %.3f", result_CommonVertex_BS_FSR.at(1));
		tex1 = new TLatex(0.6,0.22, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kYellow+2);
		tex1->Draw();
		
		tex1 = new TLatex(0.15,0.47, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		ciao = Form("FSR = %.3f", result_baseline.at(0));
		tex1 = new TLatex(0.15,0.42, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", result_CommonVertex.at(0));
		tex1 = new TLatex(0.15,0.37, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", result_CommonVertex_FSR.at(0));
		tex1 = new TLatex(0.15,0.32, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kRed+2);
		tex1->Draw();
		
		ciao = Form("Vertex: BS = %.3f", result_CommonVertex_BS.at(0));
		tex1 = new TLatex(0.15,0.27, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kBlue);
		tex1->Draw();

		ciao = Form("Vertex FSR: BS = %.3f", result_CommonVertex_BS_FSR.at(0));
		tex1 = new TLatex(0.15,0.22, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kYellow+2);
		tex1->Draw();
	}	
	if(fit == 2){

		TLatex* tex1 = new TLatex(0.5,0.67, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", h1->GetMean());
		tex1 = new TLatex(0.5,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetMean());
		tex1 = new TLatex(0.5,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetMean());
		tex1 = new TLatex(0.5,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kRed+2);
		tex1->Draw();
		
		ciao = Form("Vertex: BS = %.3f", h4->GetMean());
		tex1 = new TLatex(0.5,0.47, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kBlue);
		tex1->Draw();
	}
	canvas->Update();

	TH1F *ratio_2 = (TH1F*) h2->Clone();
	TH1F *ratio_3 = (TH1F*) h3->Clone();
	TH1F *ratio_4 = (TH1F*) h4->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2->Divide(h1);
	ratio_2->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->SetLineColor(1);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / GEN");//Reco / Gen");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	
	ratio_3->Divide(h1);
	ratio_3->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_3->SetTitle("");
	ratio_3->SetStats(0);
	ratio_3->SetLineColor(2);
	ratio_3->GetYaxis()->SetTitleSize(0.2);
	ratio_3->GetYaxis()->SetLabelSize(0.14);
	ratio_3->GetYaxis()->SetTitle("X / GEN");//Reco / Gen");
	ratio_3->GetYaxis()->SetTitleOffset(0.50);
	ratio_3->GetYaxis()->SetNdivisions(506); 
	
	ratio_4->Divide(h1);
	ratio_4->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_4->SetTitle("");
	ratio_4->SetStats(0);
	ratio_4->SetLineColor(3);
	ratio_4->GetYaxis()->SetTitleSize(0.2);
	ratio_4->GetYaxis()->SetLabelSize(0.14);
	ratio_4->GetYaxis()->SetTitle("X / GEN");//Reco / Gen");
	ratio_4->GetYaxis()->SetTitleOffset(0.50);
	ratio_4->GetYaxis()->SetNdivisions(506); 

	ratio_2->Draw();
// 	ratio_3->Draw("same");
	ratio_4->Draw("same");
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
	TLegend *legend_comparison = new TLegend(0.75,0.55,0.9,0.9);
	legend_comparison->AddEntry(ratio_2, "Base");
// 	legend_comparison->AddEntry(ratio_3, "VX");
	legend_comparison->AddEntry(ratio_4, "VX+BS");
	legend_comparison->Draw();

	canvas->Print(save);
	
}

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, TString nome_canvas, TString save, TString x_name, TLegend *legend, int fit, bool LogY){
	
	std::vector<float> result_baseline;
	std::vector<float> result_CommonVertex;
	std::vector<float> result_CommonVertex_BS;

	if(fit == 1){
		result_baseline = FitMass(h1, "Baseline", save);
		result_CommonVertex = FitMass(h2, "VX", save);
		result_CommonVertex_BS = FitMass(h3, "VX+BS", save);
	}

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
   		min = 0.000001;
   	}
	else
		min = 0;
	
// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());

	float tmp_max;
	tmp_max = Max(h1->GetMaximum(), h2->GetMaximum());
	max = Max(h3->GetMaximum(), tmp_max);
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;	
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h3->SetTitle(nome_canvas);
	h1->SetLineColor(1);
	h2->SetLineColor(2);
	h3->SetLineColor(3);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h3->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	h1->GetXaxis()->SetTitle(x_name);
	h2->GetXaxis()->SetTitle(x_name);
	h3->GetXaxis()->SetTitle(x_name);
	
	canvas->Update();
	
	h1->Draw("E");
// 	h2->Draw("Same E");
	h3->Draw("Same E");

	legend->Draw();	

	if(fit == 0){

		TLatex* tex1 = new TLatex(0.6,0.67, "RMS:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", h1->GetRMS());
		tex1 = new TLatex(0.6,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetRMS());
		tex1 = new TLatex(0.6,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetRMS());
		tex1 = new TLatex(0.6,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();

		tex1 = new TLatex(0.15,0.67, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		ciao = Form("FSR = %.3f", h1->GetMean());
		tex1 = new TLatex(0.15,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetMean());
		tex1 = new TLatex(0.15,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetMean());
		tex1 = new TLatex(0.15,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();
	}
	if(fit == 1){

		TLatex* tex1 = new TLatex(0.65,0.47, "#sigma:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("Base = %.3f", result_baseline.at(1));
		tex1 = new TLatex(0.65,0.42, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(1);
		tex1->Draw();

		ciao = Form("VX = %.3f", result_CommonVertex.at(1));
		tex1 = new TLatex(0.65,0.37, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
// 		tex1->Draw();

		ciao = Form("VX+BS = %.3f", result_CommonVertex_BS.at(1));
		tex1 = new TLatex(0.65,0.32, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();
		
		tex1 = new TLatex(0.15,0.47, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		ciao = Form("Base = %.3f", result_baseline.at(0));
		tex1 = new TLatex(0.15,0.42, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(1);
		tex1->Draw();

		ciao = Form("VX = %.3f", result_CommonVertex.at(0));
		tex1 = new TLatex(0.15,0.37, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
// 		tex1->Draw();

		ciao = Form("VX+BS = %.3f", result_CommonVertex_BS.at(0));
		tex1 = new TLatex(0.15,0.32, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();

	}	
	if(fit == 2){

		TLatex* tex1 = new TLatex(0.5,0.67, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", h1->GetMean());
		tex1 = new TLatex(0.5,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetMean());
		tex1 = new TLatex(0.5,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetMean());
		tex1 = new TLatex(0.5,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();

	}
	canvas->Update();

	TH1F *ratio_2 = (TH1F*) h2->Clone();
	TH1F *ratio_3 = (TH1F*) h3->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2->Divide(h1);
	ratio_2->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_3->SetLineColor(2);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / Baseline");//Reco / Gen");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	
	ratio_3->Divide(h1);
	ratio_3->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_3->SetTitle("");
	ratio_3->SetStats(0);
	ratio_3->SetLineColor(3);
	ratio_3->GetYaxis()->SetTitleSize(0.2);
	ratio_3->GetYaxis()->SetLabelSize(0.14);
	ratio_3->GetYaxis()->SetTitle("X / Baseline");//Reco / Gen");
	ratio_3->GetYaxis()->SetTitleOffset(0.50);
	ratio_3->GetYaxis()->SetNdivisions(506); 

// 	ratio_2->Draw();
// 	ratio_3->Draw("same");
	ratio_3->Draw();
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
	TLegend *legend_comparison = new TLegend(0.75,0.55,0.9,0.9);
// 	legend_comparison->AddEntry(ratio_2, "VX");
	legend_comparison->AddEntry(ratio_3, "VX+BS");
	legend_comparison->Draw();

	canvas->Print(save);
	
}








void Draw_TH1(RooRealVar* rv_Mass_1, RooDataSet* Mass_1, TString title, TString save, float minRange, float maxRange){
			
	std::vector<float> param;
	
	RooDataSet h_Higgs = RooDataSet(Mass_1->GetName(), Mass_1->GetTitle(), Mass_1, *Mass_1->get());//, "1", "weight");

	RooRealVar Mean("Mean", "Mean", 0, -0.5, 0.5);
	RooRealVar Sigma("Sigma", "Sigma", 0.05, -0.5, 0.5);

	RooGaussian Gauss("Gauss", "Gauss", *rv_Mass_1, Mean, Sigma);

	TCanvas *c_MC = new TCanvas(title, title, 900, 700);
	c_MC->SetFrameFillColor(0);
	c_MC->cd(1)->SetBottomMargin(0.2);
// 	c_MC->SetLogy();
	RooPlot* xframe = rv_Mass_1->frame(Title(title));
	h_Higgs.plotOn(xframe);
	Gauss.fitTo(h_Higgs, Range(minRange, maxRange));
	Gauss.plotOn(xframe,RooFit::LineColor(kBlue));
	RooPlot *framePull_DATA = rv_Mass_1->frame("");
	framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
	Gauss.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

	xframe->Draw();
	
	c_MC->Print(save);

}

void Draw(TH1F *h1, TH1F *h2, TH1F *h3, RooRealVar* rv_Mass_1, RooRealVar* rv_Mass_2, RooRealVar* rv_Mass_3, RooDataSet* Mass_1, RooDataSet* Mass_2, RooDataSet* Mass_3, TString nome_canvas, TString save, TString x_name, TLegend *legend, int fit, bool LogY){
	
	std::vector<float> result_baseline;
	std::vector<float> result_CommonVertex;
	std::vector<float> result_CommonVertex_BS;

	if(fit == 1){
		result_baseline = FitMass(rv_Mass_1, Mass_1, "Baseline", save);
		result_CommonVertex = FitMass(rv_Mass_2, Mass_2, "VX", save);
		result_CommonVertex_BS = FitMass(rv_Mass_3, Mass_3, "VX+BS", save);
	}

	TCanvas *canvas = new TCanvas(nome_canvas, nome_canvas, 750, 750);
   	canvas->cd();
   	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
   	
	Double_t min, max;

   	if(LogY){
   		pad11->SetLogy();
   		min = 0.000001;
   	}
	else
		min = 0;
	
// 	h1->Scale(1/h1->Integral());
// 	h2->Scale(1/h2->Integral());
// 	h3->Scale(1/h3->Integral());

	float tmp_max;
	tmp_max = Max(h1->GetMaximum(), h2->GetMaximum());
	max = Max(h3->GetMaximum(), tmp_max);
		
	max = 1.1 * max;
	
	if(max == 0) max = 1;	
	
	h1->SetTitle(nome_canvas);
	h2->SetTitle(nome_canvas);
	h3->SetTitle(nome_canvas);
	h1->SetLineColor(1);
	h2->SetLineColor(2);
	h3->SetLineColor(3);

	h1->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h2->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	h3->GetYaxis()->SetRangeUser(min, max);//SetMaximum(max);
	
	h1->GetXaxis()->SetTitle(x_name);
	h2->GetXaxis()->SetTitle(x_name);
	h3->GetXaxis()->SetTitle(x_name);
	
	canvas->Update();
	
	h1->Draw("E");
// 	h2->Draw("Same E");
	h3->Draw("Same E");

	legend->Draw();	

	if(fit == 0){

		TLatex* tex1 = new TLatex(0.6,0.67, "RMS:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", h1->GetRMS());
		tex1 = new TLatex(0.6,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetRMS());
		tex1 = new TLatex(0.6,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetRMS());
		tex1 = new TLatex(0.6,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();

		tex1 = new TLatex(0.15,0.67, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		ciao = Form("FSR = %.3f", h1->GetMean());
		tex1 = new TLatex(0.15,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetMean());
		tex1 = new TLatex(0.15,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetMean());
		tex1 = new TLatex(0.15,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();
	}
	if(fit == 1){

		TLatex* tex1 = new TLatex(0.65,0.47, "#sigma:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("Base = %.3f", result_baseline.at(1));
		tex1 = new TLatex(0.65,0.42, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(1);
		tex1->Draw();

		ciao = Form("VX = %.3f", result_CommonVertex.at(1));
		tex1 = new TLatex(0.65,0.37, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
// 		tex1->Draw();

		ciao = Form("VX+BS = %.3f", result_CommonVertex_BS.at(1));
		tex1 = new TLatex(0.65,0.32, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();
		
		tex1 = new TLatex(0.15,0.47, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		ciao = Form("Base = %.3f", result_baseline.at(0));
		tex1 = new TLatex(0.15,0.42, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(1);
		tex1->Draw();

		ciao = Form("VX = %.3f", result_CommonVertex.at(0));
		tex1 = new TLatex(0.15,0.37, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
// 		tex1->Draw();

		ciao = Form("VX+BS = %.3f", result_CommonVertex_BS.at(0));
		tex1 = new TLatex(0.15,0.32, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(kGreen+2);
		tex1->Draw();

	}	
	if(fit == 2){

		TLatex* tex1 = new TLatex(0.5,0.67, "Mean:");
		tex1->SetNDC();
		tex1->Draw();
	
		TString ciao = Form("FSR = %.3f", h1->GetMean());
		tex1 = new TLatex(0.5,0.62, ciao);
		tex1->SetNDC();
		tex1->Draw();

		ciao = Form("Vertex = %.3f", h2->GetMean());
		tex1 = new TLatex(0.5,0.57, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(2);
		tex1->Draw();

		ciao = Form("Vertex FSR = %.3f", h3->GetMean());
		tex1 = new TLatex(0.5,0.52, ciao);
		tex1->SetNDC();
		tex1->SetTextColor(3);
		tex1->Draw();

	}
	canvas->Update();

	TH1F *ratio_2 = (TH1F*) h2->Clone();
	TH1F *ratio_3 = (TH1F*) h3->Clone();

	canvas->Update();
	canvas->cd();

	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.25);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
	
	ratio_2->Divide(h1);
	ratio_2->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_2->SetTitle("");
	ratio_2->SetStats(0);
	ratio_2->SetLineColor(2);
	ratio_2->GetYaxis()->SetTitleSize(0.1);
	ratio_2->GetYaxis()->SetLabelSize(0.14);
	ratio_2->GetYaxis()->SetTitle("X / Baseline");//Reco / Gen");
	ratio_2->GetYaxis()->SetTitleOffset(0.50);
	ratio_2->GetYaxis()->SetNdivisions(506); 
	
	ratio_3->Divide(h1);
	ratio_3->GetYaxis()->SetRangeUser(0.0, 2.0);
	ratio_3->SetTitle("");
	ratio_3->SetStats(0);
	ratio_3->SetLineColor(3);
	ratio_3->GetYaxis()->SetTitleSize(0.2);
	ratio_3->GetYaxis()->SetLabelSize(0.14);
	ratio_3->GetYaxis()->SetTitle("X / Baseline");//Reco / Gen");
	ratio_3->GetYaxis()->SetTitleOffset(0.50);
	ratio_3->GetYaxis()->SetNdivisions(506); 

// 	ratio_2->Draw();
// 	ratio_3->Draw("same");
	ratio_3->Draw();
	
	canvas->Update();
	pad22->Update();
	TLine *line = new TLine(pad22->GetUxmin(), 1, pad22->GetUxmax(), 1);
	line->SetLineColor(kBlack);
	line->SetLineWidth(1);
	line->Draw();
	TLegend *legend_comparison = new TLegend(0.75,0.55,0.9,0.9);
// 	legend_comparison->AddEntry(ratio_2, "VX");
	legend_comparison->AddEntry(ratio_3, "VX+BS");
	legend_comparison->Draw();

	canvas->Print(save);
	
}
void FillHisto(bool deltaR, TH3F* h1[3][5], int A, float eta1, float eta2, float phi1, float phi2,  float pT1, float pT2, float mass, TString resonance){

	float massMin, massMax;
	massMin = 110;
	massMax = 130;
	float maxDR = 0.01;
	
	if(!deltaR) maxDR = 1000;
	
	if(pT1 > 5 && pT2 > 5){
		if(mass > massMin && mass < massMax){
			if(DeltaR1(eta1, phi1) < maxDR) h1[0][A]->Fill(eta1, pT1, pT1);
			if(DeltaR2(eta2, phi2) < maxDR) h1[0][A]->Fill(eta2, pT2, pT2);
			if(DeltaR1(eta1, phi1) < maxDR) h1[1][A]->Fill(eta1, pT1, pT1);
			if(DeltaR2(eta2, phi2) < maxDR) h1[2][A]->Fill(eta2, pT2, pT2);
		}
	}	
}
void FillHisto(bool deltaR, TH3F* h1[3][5], int A, float eta1, float eta2, float phi1, float phi2,  float pT1, float pT2, float val1, float val2, float mass, TString resonance){

	float massMin, massMax;
	massMin = 110;
	massMax = 130;
	float maxDR = 0.01;
	
	if(!deltaR) maxDR = 1000;
	
	if(pT1 > 5 && pT2 > 5){
		if(mass > massMin && mass < massMax){
			if(DeltaR1(eta1, phi1) < maxDR) h1[0][A]->Fill(eta1, pT1, val1);
			if(DeltaR2(eta2, phi2) < maxDR) h1[0][A]->Fill(eta2, pT2, val2);
			if(DeltaR1(eta1, phi1) < maxDR) h1[1][A]->Fill(eta1, pT1, val1);
			if(DeltaR2(eta2, phi2) < maxDR) h1[2][A]->Fill(eta2, pT2, val2);
		}
	}	
}
void FillHisto(bool deltaR, TH1F* h1[3], int A, float eta1, float eta2, float phi1, float phi2, float pT1, float pT2, float val1, float val2, float mass, TString resonance){

	float massMin, massMax;
	massMin = 110;
	massMax = 130;
	float maxDR = 0.01;
	
	if(!deltaR) maxDR = 1000;
	
	if(pT1 > 5 && pT2 > 5){
		if(mass > massMin && mass < massMax){
			if(DeltaR1(eta1, phi1) < maxDR) h1[A]->Fill(val1);
			if(DeltaR2(eta2, phi2) < maxDR) h1[A]->Fill(val2);
		}
	}	



}
float MCGen(float MC, float GEN){

	return (MC - GEN)/GEN;
}
float DeltaR1(float eta1, float phi1){

	float deltaEta = eta1-genLep_eta1;
	float DeltaPhi = phi1-genLep_phi1;
	float R = sqrt(pow(deltaEta, 2) + pow(DeltaPhi, 2));
	
	return R;
}
float DeltaR2(float eta1, float phi1){

	float deltaEta = eta1-genLep_eta2;
	float DeltaPhi = phi1-genLep_phi2;
	float R = sqrt(pow(deltaEta, 2) + pow(DeltaPhi, 2));
	
	return R;
}
float Max(float a, float b){

	if(a > b) return a;
	else return b;

}
std::vector<float> FitMass(TH1F* h1, TString title, TString save_name){
			
	std::vector<float> param;
	
	RooRealVar x("x [cm]", "Mass (GeV/c^{2})", 105, 140);
	RooDataHist h_Higgs("h_Hboson", "h_Hboson", x, h1);


	RooRealVar Mean("Mean", "Mean", 125, 120, 130);
	RooRealVar Sigma("Sigma", "Sigma", 1, 0, 30);//sigma[decay]);
	RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0, 30);//alphaL[decay]);
	RooRealVar ExpL("ExpL", "ExpL", 1, 0, 30);//expL[decay]);
	RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0, 30);//alphaR[decay]);
	RooRealVar ExpR("ExpR", "ExpR", 1, 1, 50);//expR[decay]);

	RooMyPDF_DSCB DSCB("DSCB", "DSCB", x, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);

	TCanvas *c_MC = new TCanvas(title, title, 900, 700);
	c_MC->SetFrameFillColor(0);
	c_MC->cd(1)->SetBottomMargin(0.2);
// 	c_MC->SetLogy();
	RooPlot* xframe = x.frame(Title(title));
	h_Higgs.plotOn(xframe);
	DSCB.fitTo(h_Higgs, Range(105, 140));
	DSCB.plotOn(xframe,RooFit::LineColor(kBlue));
	RooPlot *framePull_DATA = x.frame("");
	framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
	DSCB.plotOn(xframe, Components("bkg"), LineColor(kBlue), LineStyle(kDashed));
	DSCB.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

	xframe->Draw();
	
	TPad* padPull_DATA =  new  TPad("padPull","padPull",0.,0.,1.,0.2);
	padPull_DATA->Draw();
	padPull_DATA->cd(0);
	framePull_DATA->GetYaxis()->SetLabelSize(0.1);
	framePull_DATA->GetXaxis()->SetLabelSize(0.1);
	framePull_DATA->SetMinimum(-5.);
	framePull_DATA->SetMaximum(5.);
	framePull_DATA->Draw();
	TLine* lineRef = new TLine(105,0,140,0.);
	lineRef->Draw("same");

	c_MC->Print(save_name);// + ".pdf");
	
	param.push_back(Mean.getVal());
	param.push_back(Sigma.getVal());
	
	return 	param;

}

std::vector<float> FitMass(RooRealVar* rv_Mass, RooDataSet* Data_Mass, TString title, TString save_name){
			
	std::vector<float> param;
	
	RooDataSet h_Higgs = RooDataSet(Data_Mass->GetName(), Data_Mass->GetTitle(), Data_Mass, *Data_Mass->get());//, "1", "weight");


	RooRealVar Mean("Mean", "Mean", 125, 120, 130);
	RooRealVar Sigma("Sigma", "Sigma", 1, 0, 30);//sigma[decay]);
	RooRealVar AlphaL("AlphaL", "AlphaL", 1, 0, 30);//alphaL[decay]);
	RooRealVar ExpL("ExpL", "ExpL", 1, 0, 30);//expL[decay]);
	RooRealVar AlphaR("AlphaR", "AlphaR", 1, 0, 30);//alphaR[decay]);
	RooRealVar ExpR("ExpR", "ExpR", 1, 1, 50);//expR[decay]);

	RooMyPDF_DSCB DSCB("DSCB", "DSCB", *rv_Mass, Mean, Sigma, AlphaL, ExpL, AlphaR, ExpR);

	TCanvas *c_MC = new TCanvas(title, title, 900, 700);
	c_MC->SetFrameFillColor(0);
	c_MC->cd(1)->SetBottomMargin(0.2);
// 	c_MC->SetLogy();
	RooPlot* xframe = rv_Mass->frame(Title(title));
	h_Higgs.plotOn(xframe);
	DSCB.fitTo(h_Higgs, Range(105, 140));
	DSCB.plotOn(xframe,RooFit::LineColor(kBlue));
	RooPlot *framePull_DATA = rv_Mass->frame("");
	framePull_DATA->addObject((TObject*)xframe->pullHist(), "p");
// 	DSCB.plotOn(xframe, Components("bkg"), LineColor(kBlue), LineStyle(kDashed));
	DSCB.paramOn(xframe, RooFit::Layout(0.13, 0.5, 0.80));

	xframe->Draw();

	Double_t integral = h_Higgs.numEntries();
	TString integ = Form("Integral = %.0f", integral);
	TLatex* tex1 = new TLatex(0.6,0.3, integ);
	tex1->SetNDC();
	tex1->Draw();										
	
	TPad* padPull_DATA =  new  TPad("padPull","padPull",0.,0.,1.,0.2);
	padPull_DATA->Draw();
	padPull_DATA->cd(0);
	framePull_DATA->GetYaxis()->SetLabelSize(0.1);
	framePull_DATA->GetXaxis()->SetLabelSize(0.1);
	framePull_DATA->SetMinimum(-5.);
	framePull_DATA->SetMaximum(5.);
	framePull_DATA->Draw();
	TLine* lineRef = new TLine(105,0,140,0.);
	lineRef->Draw("same");

	c_MC->Print(save_name);// + ".pdf");
	
	param.push_back(Mean.getVal());
	param.push_back(Sigma.getVal());
	
	return 	param;

}


void Draw(bool dscb, bool sigma, int lepton_size, TH1F* histo[3], TLegend* legend_width, TLegend* legend_comparison, TString save_name, int m, TLine* line){

	if(!dscb && !sigma){
		for(int i = 0; i < histo[0]->GetNbinsX(); i++){
			std::cout<<histo[0]->GetBinContent(i)<<"\t"<<histo[1]->GetBinContent(i)<<"\t"<<histo[2]->GetBinContent(i)<<"\t"<<std::endl;
		}	
	}
	
	TString gd, sm;
	if(dscb)
		gd = "DSCB";
	else
		gd = "Gauss";
		
	if(sigma)
		sm = "Width";
	else
		sm = "Mean";
		
	TCanvas* canvas = new TCanvas(gd + " " + sm, gd + " " + sm, 750, 500);
   	canvas->cd();
	TPad* pad11 = new TPad("pad1", "pad1", 0, 0.25, 1, 1);
   	pad11->SetGrid();
   	pad11->SetBottomMargin(0.1);
   	pad11->Draw();
   	pad11->cd();
   	pad11->SetTicks();
		
	canvas->Update();
	
	int t;
	for(int h = 1; h < lepton_size; h++){
		t = h-1;
		histo[t]->SetLineColor(h);			
		histo[t]->SetMarkerColor(h);
		histo[t]->SetMarkerStyle(20);
		histo[t]->GetYaxis()->SetRangeUser(0.005, 0.04);
		if(sigma){
			histo[t]->GetYaxis()->SetTitle(gd + " #sigma");
			histo[t]->SetTitle(gd + " #sigma");
		}
		else{
			histo[t]->GetYaxis()->SetTitle(gd + " mean");
			histo[t]->SetTitle(gd + " mean");
		}
		histo[t]->GetXaxis()->SetTitle("p_{T} and #eta bins");

		if(t == 0)
			histo[t]->Draw();
		else{
			histo[t]->Draw("same");
		}
	}
	legend_width->Draw();

	canvas->Update();		
	canvas->cd();
		
	TPad* pad22 = new TPad("pad2", "pad2", 0, 0.01, 1, 0.285);
	pad22->SetGrid();
	pad22->Draw();
	pad22->cd();
	pad22->SetTicks();
		
	for(int h = 2; h < lepton_size; h++){
		t = h-1;
		TH1F* ratio = (TH1F*) histo[t]->Clone();
		ratio->Divide(histo[0]);
		ratio->GetYaxis()->SetRangeUser(0.6, 1.3);
		ratio->SetTitle("");
		ratio->SetStats(0);
		ratio->GetYaxis()->SetTitleSize(0.1);
		ratio->GetYaxis()->SetLabelSize(0.14);
		ratio->GetYaxis()->SetTitle("X / Baseline");
		ratio->GetYaxis()->SetTitleOffset(0.40);
		ratio->GetYaxis()->SetNdivisions(506); 
		ratio->SetLineColor(h);
		ratio->SetMarkerColor(h);
		if(t == 1)
			ratio->Draw();
		else
			ratio->Draw("same");
	}
	legend_comparison->Draw();
	
		canvas->Update();
		pad22->Update();

		line->Draw();

	save_name += gd + sm;
	
	if(m == 0)
		canvas->Print(save_name + ".pdf");
	if(m == 1)
		canvas->Print(save_name + "_shifted.pdf");
}
