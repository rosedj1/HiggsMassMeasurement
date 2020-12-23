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
#include <vector>
#include <fstream>
//
#include "TRandom3.h"
  
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "RooBreitWigner.h"
#include "RooProdPdf.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooCBShape.h"
#include "RooMinuit.h"
#include "RooFormulaVar.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"

#include "RooPlot.h"

using namespace std;

void MioSkim_2L(
  TString fs = "2mu",  TString year = "2018",  
  bool max_entry = true, int set_low = 0, int set_max = 500000, int print_every = 100000,
  TString run_on_gen = "yes",
  TString subnameIn = "",
  TString subnameOut = "",
  TString resonance = "") {

  // TString infilename = "/cmsuf/data/store/user/t2/users/ferrico/SingleBS_studies/After/";
  // TString infilename = "/cmsuf/data/store/user/t2/users/rosedj1/Higgs/HiggsMassMeas/NTuples/Hmumu/";
  TString infilename = "/cmsuf/data/store/user/t2/users/rosedj1/Higgs/HiggsMassMeas/NTuples/Hmumu/GluGluHToMuMu_M125_TuneCP5_PSweights_13TeV_amcatnloFXFX_pythia8_RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2.root";
  TString new_file = "/cmsuf/data/store/user/t2/users/rosedj1/Higgs/HiggsMassMeas/NTuples/Hmumu/MC2018Hmumu_liteskim_fullstats.root";

	// std::vector<TString> DataSetType;
	// DataSetType.clear();
	// DataSetType.push_back("SingleMuon");	
  // DataSetType.push_back("Hmumu");	
	
//   std::vector<float>* Z_mass; 
//   std::vector<float>* Z_noFSR_mass;
//   std::vector<float>* Z_massErr;
// 
//   std::string *triggersPassed;
  ULong64_t Run, LumiSect, Event;
// 
  std::vector<double> *commonPV_x = 0;
  std::vector<double> *commonPV_y = 0;
  std::vector<double> *commonPV_z = 0;
  std::vector<double> *commonBS_x = 0;
  std::vector<double> *commonBS_y = 0;
  std::vector<double> *commonBS_z = 0;
  
  std::vector<float> *lep_d0BS = 0;
  std::vector<float> *lep_d0PV = 0;
  std::vector<int> *lep_id = 0;
  std::vector<int> *lep_Sip = 0;
  std::vector<int> *lep_tightId = 0;
  std::vector<float> *lep_pt = 0;
  // std::vector<float> *lep_pt_FromMuonBestTrack = 0;  
  std::vector<float> *lep_eta = 0;
  std::vector<float> *lep_phi = 0;
  std::vector<float>* lep_mass = 0;
  std::vector<float> *lepFSR_pt = 0; 
  std::vector<float> *lepFSR_eta = 0;
  std::vector<float> *lepFSR_phi = 0;
  std::vector<float>* lepFSR_mass = 0;
  std::vector<int> *lep_genindex = 0;
  std::vector<float> *lep_RelIso = 0;
  std::vector<float> *lep_pterr = 0;
  std::vector<float> *lep_pterrold = 0;
  std::vector<int> *lep_ecalDriven = 0; 

  // std::vector<double> *vtxRecoLep_BS_d0 = 0;
  // std::vector<double> *vtxRecoLep_BS_pt = 0;
  // std::vector<double> *vtxRecoLep_BS_ptError = 0;
  // std::vector<double> *vtxRecoLep_BS_eta = 0;
  // std::vector<double> *vtxRecoLep_BS_phi = 0;
  // std::vector<double> *vtxRecoLep_BS_mass = 0;
  // std::vector<double> *vtxRecoLep_pt = 0;
  // std::vector<double> *vtxRecoLep_ptError = 0;
  // std::vector<double> *vtxRecoLep_eta = 0;
  // std::vector<double> *vtxRecoLep_phi = 0;
  // std::vector<double> *vtxRecoLep_mass = 0;
  
  std::vector<float> *lep_numberOfValidPixelHits = 0;
  std::vector<float> *lep_trackerLayersWithMeasurement = 0;
  
  // float_t mass2l_vtx_BS;
  // float_t mass2l_vtx;
  // float_t massZ_vtx_chi2;
  // float_t massZ_vtx_chi2_BS;

  std::vector<float>  *lep_dataMC = 0;
  
  std::vector<double> *fsrPhotons_pt = 0;
  std::vector<double> *fsrPhotons_eta = 0;
  std::vector<double> *fsrPhotons_phi = 0;
  
  std::vector<int> *lep_matchedR03_MomId = 0;
  std::vector<int> *lep_matchedR03_MomMomId = 0;

  float_t eventWeight;

  std::vector<float>  *GENZ_mass = 0;
  std::vector<float>  *GENlep_pt = 0;
  std::vector<float>  *GENlep_eta = 0;
  std::vector<float>  *GENlep_phi = 0;
  std::vector<float>  *GENlep_mass = 0;  
  float GenWeight;

  // float PV_x, PV_y, PV_z;   
  // float BS_x, BS_y, BS_z; 
  
  double massZ, massZErr, massZErrOld;
  double massZ_FSR, massZErr_FSR;
  // double massZ_vtx, massZ_vtx_FSR, massZErr_vtx, massZErr_vtx_FSR, massZ_vtxChi2;
  // double massZ_vtx_BS, massZ_vtx_BS_FSR, massZErr_vtx_BS, massZErr_vtx_BS_FSR, massZ_vtxChi2_BS;
  double Iso1, Iso2;
  double d0BS1, d0BS2;
  double d0PV1, d0PV2;
  int Id1, Id2;
  int Tight1, Tight2;
  // double PvX, PvY, PvZ;
  // double BsX, BsY, BsZ;
  // double PvX1, PvY1, PvZ1, PvX2, PvY2, PvZ2;
  // double BsX1, BsY1, BsZ1, BsX2, BsY2, BsZ2;
  double Pixel1, Pixel2, Tracker1, Tracker2;
  double pT1, pT2, eta1, eta2, m1, m2, phi1, phi2;
  // double pT1_FromMuonBestTrack, pT2_FromMuonBestTrack;
  double pT_FSR1, pT_FSR2, eta_FSR1, eta_FSR2, m_FSR1, m_FSR2, phi_FSR1, phi_FSR2;
  // double vtx_pT1, vtx_pT2, vtx_eta1, vtx_eta2, vtx_m1, vtx_m2, vtx_phi1, vtx_phi2;
  // double vtx_BS_pT1, vtx_BS_pT2, vtx_BS_eta1, vtx_BS_eta2, vtx_BS_m1, vtx_BS_m2, vtx_BS_phi1, vtx_BS_phi2;
  // double vtx_pT_FSR1, vtx_pT_FSR2, vtx_eta_FSR1, vtx_eta_FSR2, vtx_m_FSR1, vtx_m_FSR2, vtx_phi_FSR1, vtx_phi_FSR2;
  // double vtx_BS_pT_FSR1, vtx_BS_pT_FSR2, vtx_BS_eta_FSR1, vtx_BS_eta_FSR2, vtx_BS_m_FSR1, vtx_BS_m_FSR2, vtx_BS_phi_FSR1, vtx_BS_phi_FSR2;
  // double d0BS_vtx_BS1, d0BS_vtx_BS2;

  double pterr1, pterr2;
  // double pterr1_VX, pterr2_VX;
  // double pterr1_VX_BS, pterr2_VX_BS;
  double pterr1old, pterr2old;  
  double genzm, GENmass2l;
  double weight;
  double genLep_pt1=-999, genLep_pt2=-999;
  double genLep_eta1=-999, genLep_eta2=-999;
  double genLep_phi1=-999, genLep_phi2=-999;
  double genLep_mass1=-999, genLep_mass2=-999;
  int lep1_ecalDriven = -1, lep2_ecalDriven = -1;
  int nFSRPhotons;
  double Met;
  float_t met;
  bool passedTrig;
  bool passedFullSelection; 
  float massZ1, massZ2;
  int i_max;
  int LepMomId1;
  int LepMomId2;
  int LepGrandMomId1;
  int LepGrandMomId2;
  int nVtx, nInt;
  int NVtx, NInt;
  
  int count_Id1, count_Id2;
  count_Id1 = 0;
  count_Id2 = 0;
  
  if(fs!="2e" && fs!="2mu"){
  	cout<<"fs has to be 2e, or 2mu"<<endl;
  	return;
  }

  // if(run_on_gen == "no"){
  // 	infilename = infilename + "Hmumu_" + year + "_" + subnameIn + ".root";  
  // 	new_file =  "./DY_" + year + "_" + subnameOut + "_cancella.root";	
  // }
  // else{
  // 	infilename = infilename + DataSetType.at(0) + "_Run" + year + run_on_gen + ".root";
  // 	new_file = "./Data_" + DataSetType.at(0) + "_Run" + year + run_on_gen + "_skimmed_" + subnameOut + ".root";
  // }
  
  // TString infilename = indir + infile;
  cout<<"---- fs is "<<fs<<endl;
  cout<<"---- reading file:\n"<< infilename << endl;
  TFile* infile = infile = TFile::Open(infilename);
  
  TTree* tree; 
  if(infile){ 
  	std::cout<<"File trovato (opened). All is well."<<std::endl;
    infile->cd("Ana");
    tree = (TTree*)gDirectory->Get("passedEvents");
  }
  else{ std::cout<<"ERROR could not find the file"<<std::endl; return -1;}

  if(!tree) { cout<<"ERROR could not find the tree"<<endl; return -1;}

 std::cout << "Creating new file:\n" << new_file << std::endl; 
  TFile* tmpFile = new TFile(new_file,"RECREATE");

  TTree* newtree = new TTree("passedEvents","passedEvents");

  cout<<"start setting tree "<<endl;
  
  newtree->Branch("GENmass2l", &GENmass2l,"GENmass2l/D");
  // newtree->Branch("genzm", &genzm,"genzm/D");
  newtree->Branch("genLep_pt1", &genLep_pt1, "genLep_pt1/D");
  newtree->Branch("genLep_pt2", &genLep_pt2, "genLep_pt2/D");
  newtree->Branch("genLep_eta1", &genLep_eta1, "genLep_eta1/D");
  newtree->Branch("genLep_eta2", &genLep_eta2, "genLep_eta2/D");
  newtree->Branch("genLep_phi1", &genLep_phi1, "genLep_phi1/D");
  newtree->Branch("genLep_phi2", &genLep_phi2, "genLep_phi2/D");
  newtree->Branch("genLep_mass1", &genLep_mass1, "genLep_mass1/D");
  newtree->Branch("genLep_mass2", &genLep_mass2, "genLep_mass2/D");

  newtree->Branch("lep_genindex", &lep_genindex);
  
  newtree->Branch("NInt",&NInt,"NInt/I");
  newtree->Branch("NVtx",&NVtx,"NVtx/I");  

  // newtree->Branch("PvX",&PvX,"PvX/D");
  // newtree->Branch("PvY",&PvY,"PvY/D");
  // newtree->Branch("PvZ",&PvZ,"PvZ/D");
  // newtree->Branch("BsX",&BsX,"BsX/D");
  // newtree->Branch("BsY",&BsY,"BsY/D");
  // newtree->Branch("BsZ",&BsZ,"BsZ/D");

  newtree->Branch("massZ",&massZ,"massZ/D");
  newtree->Branch("massZErr",&massZErr,"massZErr/D");

  newtree->Branch("massZ_FSR",&massZ_FSR,"massZ_FSR/D");
  newtree->Branch("massZErr_FSR",&massZErr_FSR,"massZErr_FSR/D");

  // newtree->Branch("massZ_vtx",&massZ_vtx,"massZ_vtx/D");
  // newtree->Branch("massZ_vtx_FSR",&massZ_vtx_FSR,"massZ_vtx_FSR/D");
  // newtree->Branch("massZErr_vtx",&massZErr_vtx,"massZErr_vtx/D");
  // newtree->Branch("massZErr_vtx_FSR",&massZErr_vtx_FSR,"massZErr_vtx_FSR/D");
  // newtree->Branch("massZ_vtxChi2",&massZ_vtxChi2,"massZ_vtxChi2/D");

  // newtree->Branch("massZ_vtx_BS",&massZ_vtx_BS,"massZ_vtx_BS/D");
  // newtree->Branch("massZ_vtx_BS_FSR",&massZ_vtx_BS_FSR,"massZ_vtx_BS_FSR/D");
  // newtree->Branch("massZErr_vtx_BS",&massZErr_vtx_BS,"massZErr_vtx_BS/D");
  // newtree->Branch("massZErr_vtx_BS_FSR",&massZErr_vtx_BS_FSR,"massZErr_vtx_BS_FSR/D");
  // newtree->Branch("massZ_vtxChi2_BS",&massZ_vtxChi2_BS,"massZ_vtxChi2_BS/D");

//   newtree->Branch("massZErrOld",&massZErrOld,"massZErrOld/D");
//  newtree->Branch("PvX1",&PvX1,"PvX1/D");
//  newtree->Branch("PvX2",&PvX2,"PvX2/D");
//  newtree->Branch("PvY1",&PvY1,"PvY1/D");
//  newtree->Branch("PvY2",&PvY2,"PvY2/D");
//  newtree->Branch("PvZ1",&PvZ1,"PvZ1/D");
//  newtree->Branch("PvZ2",&PvZ2,"PvZ2/D");
//  newtree->Branch("BsX1",&BsX1,"BsX1/D");
//  newtree->Branch("BsX2",&BsX2,"BsX2/D");
//  newtree->Branch("BsY1",&BsY1,"BsY1/D");
//  newtree->Branch("BsY2",&BsY2,"BsY2/D");
//  newtree->Branch("BsZ1",&BsZ1,"BsZ1/D");
//  newtree->Branch("BsZ2",&BsZ2,"BsZ2/D");

  // newtree->Branch("pT1_FromMuonBestTrack",&pT1_FromMuonBestTrack,"pT1_FromMuonBestTrack/D");
  // newtree->Branch("pT2_FromMuonBestTrack",&pT2_FromMuonBestTrack,"pT2_FromMuonBestTrack/D");
  newtree->Branch("pT1",&pT1,"pT1/D");
  newtree->Branch("pT2",&pT2,"pT2/D");
  newtree->Branch("eta1",&eta1,"eta1/D");
  newtree->Branch("eta2",&eta2,"eta2/D");
  newtree->Branch("phi1",&phi1,"phi1/D");
  newtree->Branch("phi2",&phi2,"phi2/D");
  newtree->Branch("m1",&m1,"m1/D");
  newtree->Branch("m2",&m2,"m2/D");

  newtree->Branch("pT_FSR1",&pT_FSR1,"pT_FSR1/D");
  newtree->Branch("pT_FSR2",&pT_FSR2,"pT_FSR2/D");
  newtree->Branch("eta_FSR1",&eta_FSR1,"eta_FSR1/D");
  newtree->Branch("eta_FSR2",&eta_FSR2,"eta_FSR2/D");
  newtree->Branch("phi_FSR1",&phi_FSR1,"phi_FSR1/D");
  newtree->Branch("phi_FSR2",&phi_FSR2,"phi_FSR2/D");
  newtree->Branch("m_FSR1",&m_FSR1,"m_FSR1/D");
  newtree->Branch("m_FSR2",&m_FSR2,"m_FSR2/D");

//   newtree->Branch("Pixel1",&Pixel1,"Pixel1/D");
//   newtree->Branch("Pixel2",&Pixel2,"Pixel2/D");
//   newtree->Branch("Tracker1",&Tracker1,"Tracker1/D");
//   newtree->Branch("Tracker2",&Tracker2,"Tracker2/D");

  // newtree->Branch("vtx_pT1",&vtx_pT1,"vtx_pT1/D");
  // newtree->Branch("vtx_pT2",&vtx_pT2, "vtx_pT2/D");
  // newtree->Branch("vtx_eta1",&vtx_eta1, "vtx_eta1/D");
  // newtree->Branch("vtx_eta2",&vtx_eta2, "vtx_eta2/D");
  // newtree->Branch("vtx_phi1",&vtx_phi1, "vtx_phi1/D");
  // newtree->Branch("vtx_phi2",&vtx_phi2, "vtx_phi2/D");
  // newtree->Branch("vtx_m1",&vtx_m1, "vtx_m1/D");
  // newtree->Branch("vtx_m2",&vtx_m2, "vtx_m2/D");

  // newtree->Branch("vtx_BS_pT1",&vtx_BS_pT1,"vtx_BS_pT1/D");
  // newtree->Branch("vtx_BS_pT2",&vtx_BS_pT2, "vtx_BS_pT2/D");
  // newtree->Branch("vtx_BS_eta1",&vtx_BS_eta1, "vtx_BS_eta1/D");
  // newtree->Branch("vtx_BS_eta2",&vtx_BS_eta2, "vtx_BS_eta2/D");
  // newtree->Branch("vtx_BS_phi1",&vtx_BS_phi1, "vtx_BS_phi1/D");
  // newtree->Branch("vtx_BS_phi2",&vtx_BS_phi2, "vtx_BS_phi2/D");
  // newtree->Branch("vtx_BS_m1",&vtx_BS_m1, "vtx_BS_m1/D");
  // newtree->Branch("vtx_BS_m2",&vtx_BS_m2, "vtx_BS_m2/D");

//  newtree->Branch("vtx_pT_FSR1",&vtx_pT_FSR1,"vtx_pT_FSR1/D");
//   newtree->Branch("vtx_pT_FSR2",&vtx_pT_FSR2, "vtx_pT_FSR2/D");
//   newtree->Branch("vtx_eta_FSR1",&vtx_eta_FSR1, "vtx_eta_FSR1/D");
//   newtree->Branch("vtx_eta_FSR2",&vtx_eta_FSR2, "vtx_eta_FSR2/D");
//   newtree->Branch("vtx_phi_FSR1",&vtx_phi_FSR1, "vtx_phi_FSR1/D");
//   newtree->Branch("vtx_phi_FSR2",&vtx_phi_FSR2, "vtx_phi_FSR2/D");
//   newtree->Branch("vtx_m_FSR1",&vtx_m_FSR1, "vtx_m_FSR1/D");
//   newtree->Branch("vtx_m_FSR2",&vtx_m_FSR2, "vtx_m_FSR2/D");

//   newtree->Branch("vtx_BS_pT_FSR1",&vtx_BS_pT_FSR1,"vtx_BS_pT_FSR1/D");
//   newtree->Branch("vtx_BS_pT_FSR2",&vtx_BS_pT_FSR2, "vtx_BS_pT_FSR2/D");
//   newtree->Branch("vtx_BS_eta_FSR1",&vtx_BS_eta_FSR1, "vtx_BS_eta_FSR1/D");
//   newtree->Branch("vtx_BS_eta_FSR2",&vtx_BS_eta_FSR2, "vtx_BS_eta_FSR2/D");
//   newtree->Branch("vtx_BS_phi_FSR1",&vtx_BS_phi_FSR1, "vtx_BS_phi_FSR1/D");
//   newtree->Branch("vtx_BS_phi_FSR2",&vtx_BS_phi_FSR2, "vtx_BS_phi_FSR2/D");
//   newtree->Branch("vtx_BS_m_FSR1",&vtx_BS_m_FSR1, "vtx_BS_m_FSR1/D");
//   newtree->Branch("vtx_BS_m_FSR2",&vtx_BS_m_FSR2, "vtx_BS_m_FSR2/D");

  newtree->Branch("d0BS1",&d0BS1,"d0BS1/D"); 
  newtree->Branch("d0BS2",&d0BS2,"d0BS2/D"); 

  // newtree->Branch("d0BS_vtx_BS1",&d0BS_vtx_BS1,"d0BS_vtx_BS1/D"); 
  // newtree->Branch("d0BS_vtx_BS2",&d0BS_vtx_BS2,"d0BS_vtx_BS2/D"); 

//   newtree->Branch("d0PV1",&d0PV1,"d0PV1/D"); 
//   newtree->Branch("d0PV2",&d0PV2,"d0PV2/D"); 

//   newtree->Branch("Iso1",&Iso1,"Iso1/D");
//   newtree->Branch("Iso2",&Iso2,"Iso2/D");
  newtree->Branch("Id1",&Id1,"Id1/I");
  newtree->Branch("Id2",&Id2,"Id2/I");
//   newtree->Branch("Tight1",&Tight1,"Tight1/I");
//   newtree->Branch("Tight2",&Tight2,"Tight2/I");

  newtree->Branch("pterr1",&pterr1,"pterr1/D");
  newtree->Branch("pterr2",&pterr2,"pterr2/D");
  // newtree->Branch("pterr1_VX",&pterr1_VX,"pterr1_VX/D");
  // newtree->Branch("pterr2_VX",&pterr2_VX,"pterr2_VX/D");
  // newtree->Branch("pterr1_VX_BS",&pterr1_VX_BS,"pterr1_VX_BS/D");
  // newtree->Branch("pterr2_VX_BS",&pterr2_VX_BS,"pterr2_VX_BS/D");
//   newtree->Branch("pterr1old",&pterr1old,"pterr1old/D");
//   newtree->Branch("pterr2old",&pterr2old,"pterr2old/D");
//   newtree->Branch("Met", &Met, "Met/D");
  newtree->Branch("weight",&weight,"weight/D");

  newtree->Branch("nFSRPhotons", &nFSRPhotons, "nFSRPhotons/I");
if(fs == "2e"){
  newtree->Branch("lep1_ecalDriven", &lep1_ecalDriven, "lep1_ecalDriven/I");
  newtree->Branch("lep2_ecalDriven", &lep2_ecalDriven, "lep2_ecalDriven/I");
}
  
  newtree->Branch("LepMomId1", &LepMomId1, "LepMomId1/I");
  newtree->Branch("LepMomId2", &LepMomId2, "LepMomId2/I");
  newtree->Branch("LepGrandMomId1", &LepGrandMomId1, "LepGrandMomId1/I");
  newtree->Branch("LepGrandMomId2", &LepGrandMomId2, "LepGrandMomId2/I");

  cout<<"start reading tree "<<endl;
        Long64_t nentries = tree->GetEntries();
// 
        tree->SetBranchStatus("*",0);//1);  
        tree->SetBranchStatus("nInt",1);
        tree->SetBranchStatus("nVtx",1);
//         tree->SetBranchStatus("PV_x",1);     
//         tree->SetBranchStatus("PV_y",1);     
//         tree->SetBranchStatus("PV_z",1);     
//         tree->SetBranchStatus("BS_x",1);     
//         tree->SetBranchStatus("BS_y",1);     
//         tree->SetBranchStatus("BS_z",1);                          
//         tree->SetBranchStatus("passedFullSelection",1);     
//         tree->SetBranchStatus("passedTrig",1);              
//         tree->SetBranchStatus("commonPV_x",1);                  
//         tree->SetBranchStatus("commonPV_y",1);             
//         tree->SetBranchStatus("commonPV_z",1);                  
//         tree->SetBranchStatus("commonBS_x",1);               
//         tree->SetBranchStatus("commonBS_y",1);                 
//         tree->SetBranchStatus("commonBS_z",1);                          
        tree->SetBranchStatus("lep_d0BS",1);
        tree->SetBranchStatus("lep_d0PV",1);
        tree->SetBranchStatus("lep_id",1);                  
        tree->SetBranchStatus("lep_tightId",1);             
//         tree->SetBranchStatus("lep_numberOfValidPixelHits",1);
//         tree->SetBranchStatus("lep_trackerLayersWithMeasurement",1);
        // tree->SetBranchStatus("lep_pt_FromMuonBestTrack",1);                  
        tree->SetBranchStatus("lep_pt",1);                  
        tree->SetBranchStatus("lep_eta",1);                 
        tree->SetBranchStatus("lep_phi",1);                 
        tree->SetBranchStatus("lep_mass",1);                
        tree->SetBranchStatus("lepFSR_pt",1);               
        tree->SetBranchStatus("lepFSR_eta",1);                 
        tree->SetBranchStatus("lepFSR_phi",1);                 
        tree->SetBranchStatus("lepFSR_mass",1);                
        tree->SetBranchStatus("lep_RelIso",1);              
        tree->SetBranchStatus("lep_pterr",1);               
//         tree->SetBranchStatus("lep_pterrold",1);            
//         tree->SetBranchStatus("lep_Sip",1);                 
        tree->SetBranchStatus("lep_dataMC",1);              
        tree->SetBranchStatus("lep_genindex",1);            
        tree->SetBranchStatus("lep_ecalDriven", 1);
//         tree->SetBranchStatus("Run",1);                           
//         tree->SetBranchStatus("LumiSect",1);                      
//         tree->SetBranchStatus("Event",1);                         
//         tree->SetBranchStatus("met",1);                           
        tree->SetBranchStatus("GENZ_mass",1);                     
        tree->SetBranchStatus("GENlep_pt",1);                     
        tree->SetBranchStatus("GENlep_eta",1);                    
        tree->SetBranchStatus("GENlep_phi",1);                    
        tree->SetBranchStatus("GENlep_mass",1);                   
        tree->SetBranchStatus("nFSRPhotons",1);                   
        tree->SetBranchStatus("fsrPhotons_pt", 1);                   
        tree->SetBranchStatus("fsrPhotons_eta", 1);                   
        tree->SetBranchStatus("fsrPhotons_phi", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_pt", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_ptError", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_eta", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_phi", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_mass", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_BS_pt", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_BS_d0", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_BS_ptError", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_BS_eta", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_BS_phi", 1);                   
        // tree->SetBranchStatus("vtxRecoLep_BS_mass", 1);                   
        // tree->SetBranchStatus("mass2l_vtx_BS", 1);                   
        // tree->SetBranchStatus("mass2l_vtx", 1);                   
        // tree->SetBranchStatus("massZ_vtx_chi2", 1);                   
        // tree->SetBranchStatus("massZ_vtx_chi2_BS", 1);                   
		
			tree->SetBranchStatus("lep_matchedR03_MomId", 1);
			tree->SetBranchStatus("lep_matchedR03_MomMomId", 1);
			tree->SetBranchStatus("eventWeight", 1);                                                   
			tree->SetBranchStatus("GENZ_mass",1);                     
			tree->SetBranchStatus("GENlep_pt",1);                     
			tree->SetBranchStatus("GENlep_eta",1);                    
			tree->SetBranchStatus("GENlep_phi",1);                    
			tree->SetBranchStatus("GENlep_mass",1);                   

//         tree->SetBranchAddress("PV_x", &PV_x);
//         tree->SetBranchAddress("PV_y", &PV_y);   
//         tree->SetBranchAddress("PV_z", &PV_z);      
//         tree->SetBranchAddress("BS_x", &BS_x);
//         tree->SetBranchAddress("BS_y", &BS_y);
//         tree->SetBranchAddress("BS_z", &BS_z);                
//         tree->SetBranchAddress("passedFullSelection",&passedFullSelection);
//         tree->SetBranchAddress("commonPV_x",&commonPV_x);
//         tree->SetBranchAddress("commonPV_y",&commonPV_y);   
//         tree->SetBranchAddress("commonPV_z", &commonPV_z);      
//         tree->SetBranchAddress("commonBS_x", &commonBS_x);
//         tree->SetBranchAddress("commonBS_y", &commonBS_y);
//         tree->SetBranchAddress("commonBS_z", &commonBS_z);                
//         tree->SetBranchAddress("passedTrig",&passedTrig);  
        tree->SetBranchAddress("nInt", &nInt);
        tree->SetBranchAddress("nVtx", &nVtx);
 
        tree->SetBranchAddress("lep_tightId", &lep_tightId);      
        tree->SetBranchAddress("lep_d0BS", &lep_d0BS);
        tree->SetBranchAddress("lep_d0PV", &lep_d0PV);
        tree->SetBranchAddress("lep_id", &lep_id);                
//         tree->SetBranchAddress("lep_numberOfValidPixelHits", &lep_numberOfValidPixelHits);
//         tree->SetBranchAddress("lep_trackerLayersWithMeasurement", &lep_trackerLayersWithMeasurement);                
        // tree->SetBranchAddress("lep_pt_FromMuonBestTrack", &lep_pt_FromMuonBestTrack);                 
        tree->SetBranchAddress("lep_pt", &lep_pt);                 
        tree->SetBranchAddress("lep_eta",&lep_eta);               
        tree->SetBranchAddress("lep_phi",&lep_phi);               
        tree->SetBranchAddress("lep_mass",&lep_mass);             
        tree->SetBranchAddress("lepFSR_pt",&lepFSR_pt);           
        tree->SetBranchAddress("lepFSR_eta",&lepFSR_eta);               
        tree->SetBranchAddress("lepFSR_phi",&lepFSR_phi);               
        tree->SetBranchAddress("lepFSR_mass",&lepFSR_mass);             
        tree->SetBranchAddress("lep_RelIso",&lep_RelIso);         
        tree->SetBranchAddress("lep_pterr",&lep_pterr);           
//         tree->SetBranchAddress("lep_pterrold",&lep_pterrold);     
//         tree->SetBranchAddress("lep_Sip", &lep_Sip);              
        tree->SetBranchAddress("lep_dataMC", &lep_dataMC);        
        tree->SetBranchAddress("lep_genindex", &lep_genindex);    
        tree->SetBranchAddress("lep_ecalDriven", &lep_ecalDriven);  
//         tree->SetBranchAddress("Run",&Run);                       
//         tree->SetBranchAddress("LumiSect",&LumiSect);             
//         tree->SetBranchAddress("Event",&Event);                   
//         tree->SetBranchAddress("met", &met);   
        
        tree->SetBranchAddress("fsrPhotons_pt", &fsrPhotons_pt);   
        tree->SetBranchAddress("fsrPhotons_eta", &fsrPhotons_eta);   
        tree->SetBranchAddress("fsrPhotons_phi", &fsrPhotons_phi);   

        // tree->SetBranchAddress("vtxRecoLep_pt", &vtxRecoLep_pt);   
        // tree->SetBranchAddress("vtxRecoLep_ptError", &vtxRecoLep_ptError);   
        // tree->SetBranchAddress("vtxRecoLep_eta", &vtxRecoLep_eta);   
        // tree->SetBranchAddress("vtxRecoLep_phi", &vtxRecoLep_phi);   
        // tree->SetBranchAddress("vtxRecoLep_mass", &vtxRecoLep_mass);   

        // tree->SetBranchAddress("mass2l_vtx", &mass2l_vtx);
        // tree->SetBranchAddress("massZ_vtx_chi2", &massZ_vtx_chi2);   
                           
        // tree->Draw("lep_ecalDriven");

        //tree->SetBranchAddress("massZ1_vtx", &massZ1_vtx);
        //tree->SetBranchAddress("massZErr1_vtx", &massZErr1_vtx);
        //tree->SetBranchAddress("massZ1_vtx_chi2", &massZ1_vtx_chi2);

			tree->SetBranchAddress("lep_matchedR03_MomId", &lep_matchedR03_MomId);
			tree->SetBranchAddress("lep_matchedR03_MomMomId", &lep_matchedR03_MomMomId);
			tree->SetBranchAddress("eventWeight", &eventWeight);
			tree->SetBranchAddress("GENZ_mass", &GENZ_mass);
			tree->SetBranchAddress("GENlep_pt", &GENlep_pt);
			tree->SetBranchAddress("GENlep_eta", &GENlep_eta);
			tree->SetBranchAddress("GENlep_phi", &GENlep_phi);
			tree->SetBranchAddress("GENlep_mass", &GENlep_mass);
			tree->SetBranchAddress("nFSRPhotons", &nFSRPhotons);

        // return;

        std::cout<<"after tree set, total entries in tree = "<<tree->GetEntries()<<std::endl;
        
        if (max_entry) {
          std::cout << "Running over all entries." << endl;
          i_max = nentries;
        } else {
          cout << "Only running over entries in range: " << set_low << "-" << set_max << endl;
          i_max = set_max;
        }
    
        int count_trig = 0;
        
         for(int mcfmEvt_HZZ = set_low; mcfmEvt_HZZ < i_max; mcfmEvt_HZZ++) { //event loop

            tree->GetEntry(mcfmEvt_HZZ);
            if(mcfmEvt_HZZ % print_every == 0)
              std::cout<<mcfmEvt_HZZ<<" --- Dentro il tree --- "<<std::endl;

//            if(!passedTrig){ 
//		count_trig++;
//		continue;
//	    }
            if((*lep_id).size()<2) continue;
            vector<int> passLepIndex;
            for(unsigned int il=0; il<(*lep_pt).size(); il++){
           
                if(!(*lep_tightId)[il]) continue;  // 
                //  if((*lep_RelIso)[il] > 0.35) continue;  // For HZZ.
                if((*lep_RelIso)[il] > 0.25) continue;  // For Hmumu.
//                  if((*lep_Sip)[il] > 4) continue;
                 passLepIndex.push_back(il);
            }
	
            if(passLepIndex.size()!=2){ 
			count_trig++;
			continue;
	   }
  
            unsigned int L1 = passLepIndex[0];
            unsigned int L2 = passLepIndex[1];
            
            int idL1 = (*lep_id)[L1];
            int idL2 = (*lep_id)[L2];
            if((idL1+idL2)!=0) continue;
            if(fs=="2e" && abs(idL1)!=11) continue;
            if(fs=="2mu" && abs(idL1)!=13) continue;
            
//             PvX = PV_x;
//             PvY = PV_y;
//             PvZ = PV_z;
//             BsX = BS_x;
//             BsY = BS_y;
//             BsZ = BS_z;

            if(run_on_gen == "no")
	            weight = eventWeight*(*lep_dataMC)[L1]*(*lep_dataMC)[L2];
	        else
	        	weight = 1;
            
            NInt = nInt;
            NVtx = nVtx;
                    
//             PvX1 = (*commonPV_x)[L1];
//             PvY1 = (*commonPV_y)[L1]; 
//             PvZ1 = (*commonPV_z)[L1];
//             PvX2 = (*commonPV_x)[L2];
//             PvY2 = (*commonPV_y)[L2]; 
//             PvZ2 = (*commonPV_z)[L2];
// 
//             BsX1 = (*commonBS_x)[L1];
//             BsY1 = (*commonBS_y)[L1]; 
//             BsZ1 = (*commonBS_z)[L1];
//             BsX2 = (*commonBS_x)[L2];
//             BsY2 = (*commonBS_y)[L2]; 
//             BsZ2 = (*commonBS_z)[L2];

            TLorentzVector lep1(0,0,0,0);
            TLorentzVector lep2(0,0,0,0);

            d0BS1 = (*lep_d0BS)[L1]; d0PV1 = (*lep_d0PV)[L1];
            d0BS2 = (*lep_d0BS)[L2]; d0PV2 = (*lep_d0PV)[L2];
            
            
            eta1 = (*lep_eta)[L1];         eta2 = (*lep_eta)[L2];
            phi1 = double((*lep_phi)[L1]); phi2 = double((*lep_phi)[L2]); 
            m1 = double((*lep_mass)[L1]);  m2 = double((*lep_mass)[L2]);
            pT1 = (*lep_pt)[L1];           pT2 = (*lep_pt)[L2];
            
            // pT1_FromMuonBestTrack = (*lep_pt_FromMuonBestTrack)[L1];
            // pT2_FromMuonBestTrack = (*lep_pt_FromMuonBestTrack)[L2];

            eta_FSR1 = (*lepFSR_eta)[L1];         eta_FSR2 = (*lepFSR_eta)[L2];
            phi_FSR1 = double((*lepFSR_phi)[L1]); phi_FSR2 = double((*lepFSR_phi)[L2]); 
            m_FSR1 = double((*lepFSR_mass)[L1]);  m_FSR2 = double((*lepFSR_mass)[L2]);
            pT_FSR1 = (*lepFSR_pt)[L1];           pT_FSR2 = (*lepFSR_pt)[L2];

            Iso1 = (*lep_RelIso)[L1];             Iso2 = (*lep_RelIso)[L2];
            Id1 = (*lep_id)[L1];                  Id2 = (*lep_id)[L2];
            Tight1 = (*lep_tightId)[L1];          Tight2 = (*lep_tightId)[L2];
            
//             Pixel1 = (*lep_numberOfValidPixelHits)[L1]; Pixel2 = (*lep_numberOfValidPixelHits)[L2];
//             Tracker1 = (*lep_trackerLayersWithMeasurement)[L1]; Tracker2 = (*lep_trackerLayersWithMeasurement)[L2];

            lep1.SetPtEtaPhiM(double((*lep_pt)[L1]),double((*lep_eta)[L1]),double((*lep_phi)[L1]),double((*lep_mass)[L1]));
            lep2.SetPtEtaPhiM(double((*lep_pt)[L2]),double((*lep_eta)[L2]),double((*lep_phi)[L2]),double((*lep_mass)[L2]));

            massZ = (lep1+lep2).M();
            pterr1 = double((*lep_pterr)[L1]); pterr2 = double((*lep_pterr)[L2]);            
//             pterr1old = double((*lep_pterrold)[L1]); pterr2old = double((*lep_pterrold)[L2]);

            TLorentzVector lep1p, lep2p;
            lep1p.SetPtEtaPhiM(double((*lep_pt)[L1]+pterr1),double((*lep_eta)[L1]),double((*lep_phi)[L1]),double((*lep_mass)[L1]));
            lep2p.SetPtEtaPhiM(double((*lep_pt)[L2]+pterr2),double((*lep_eta)[L2]),double((*lep_phi)[L2]),double((*lep_mass)[L2]));

            double dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
            double dm2 = (lep1+lep2p).M()-(lep1+lep2).M();
 
            massZErr = TMath::Sqrt(dm1*dm1+dm2*dm2);


            lep1.SetPtEtaPhiM(double((*lepFSR_pt)[L1]),double((*lepFSR_eta)[L1]),double((*lepFSR_phi)[L1]),double((*lepFSR_mass)[L1]));
            lep2.SetPtEtaPhiM(double((*lepFSR_pt)[L2]),double((*lepFSR_eta)[L2]),double((*lepFSR_phi)[L2]),double((*lepFSR_mass)[L2]));

            massZ_FSR = (lep1+lep2).M();
            pterr1 = double((*lep_pterr)[L1]); pterr2 = double((*lep_pterr)[L2]);
//             pterr1old = double((*lep_pterrold)[L1]); pterr2old = double((*lep_pterrold)[L2]);

            lep1p.SetPtEtaPhiM(double((*lepFSR_pt)[L1]+pterr1),double((*lepFSR_eta)[L1]),double((*lepFSR_phi)[L1]),double((*lepFSR_mass)[L1]));
            lep2p.SetPtEtaPhiM(double((*lepFSR_pt)[L2]+pterr2),double((*lepFSR_eta)[L2]),double((*lepFSR_phi)[L2]),double((*lepFSR_mass)[L2]));

            dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
            dm2 = (lep1+lep2p).M()-(lep1+lep2).M();
 
            massZErr_FSR = TMath::Sqrt(dm1*dm1+dm2*dm2);

//             lep1p.SetPtEtaPhiM(double((*lepFSR_pt)[L1]+pterr1old),double((*lep_eta)[L1]),double((*lep_phi)[L1]),double((*lep_mass)[L1]));
//             lep2p.SetPtEtaPhiM(double((*lepFSR_pt)[L2]+pterr2old),double((*lep_eta)[L2]),double((*lep_phi)[L2]),double((*lep_mass)[L2]));
// 
//             dm1 = (lep1p+lep2).M()-(lep1+lep2).M();
//             dm2 = (lep1+lep2p).M()-(lep1+lep2).M();
// 
//             massZErrOld = TMath::Sqrt(dm1*dm1+dm2*dm2);
//             Met = met; 
            
			TLorentzVector photon1, photon2;
			TLorentzVector vtx_Lep1, vtx_Lep2; 
			TLorentzVector vtx_Lep1_FSR, vtx_Lep2_FSR; 



/*

			/////// vtx constraint
			vtx_eta1 = double((*vtxRecoLep_eta)[L1]); vtx_eta2 = double((*vtxRecoLep_eta)[1]);
            vtx_phi1 = double((*vtxRecoLep_phi)[L1]); vtx_phi2 = double((*vtxRecoLep_phi)[1]); 
            vtx_m1 = double((*vtxRecoLep_mass)[L1]); vtx_m2 = double((*vtxRecoLep_mass)[1]);
            vtx_pT1 = double((*vtxRecoLep_pt)[L1]); vtx_pT2 = double((*vtxRecoLep_pt)[1]);
 
            pterr1_VX = double((*vtxRecoLep_ptError)[L1]); pterr2_VX = double((*vtxRecoLep_ptError)[1]);
           
            vtx_Lep1.SetPtEtaPhiM(vtx_pT1, vtx_eta1, vtx_phi1, vtx_m1);
            vtx_Lep2.SetPtEtaPhiM(vtx_pT2, vtx_eta2, vtx_phi2, vtx_m2);

            massZ_vtx = mass2l_vtx;
            massZ_vtxChi2 = massZ_vtx_chi2;

            lep1p.SetPtEtaPhiM(double(vtx_pT1 + pterr1_VX), vtx_eta1, vtx_phi1, vtx_m1);
            lep2p.SetPtEtaPhiM(double(vtx_pT2 + pterr2_VX), vtx_eta2, vtx_phi2, vtx_m2);

            dm1 = (lep1p+vtx_Lep2).M()-(vtx_Lep1+vtx_Lep2).M();
            dm2 = (vtx_Lep1+lep2p).M()-(vtx_Lep1+vtx_Lep2).M();
 
            massZErr_vtx = TMath::Sqrt(dm1*dm1+dm2*dm2);
                        
			if((*fsrPhotons_pt).size() > 0){
	            photon1.SetPtEtaPhiM((*fsrPhotons_pt)[0], (*fsrPhotons_eta)[0], (*fsrPhotons_phi)[0], 0.0);
	            float DeltaR1 = sqrt(pow((*fsrPhotons_eta)[0] - vtx_Lep1.Eta(),2) + pow((*fsrPhotons_phi)[0] - vtx_Lep1.Phi(),2));
	            if(DeltaR1 < 0.5)
		            vtx_Lep1_FSR = vtx_Lep1 + photon1;
		        else
		        	vtx_Lep1_FSR = vtx_Lep1;
		            
				if((*fsrPhotons_pt).size() > 1){
		            photon2.SetPtEtaPhiM((*fsrPhotons_pt)[1], (*fsrPhotons_eta)[1], (*fsrPhotons_phi)[1], 0.0);
	    	        float DeltaR2 = sqrt(pow((*fsrPhotons_eta)[1] - vtx_Lep2.Eta(),2) + pow((*fsrPhotons_phi)[1] - vtx_Lep2.Phi(),2));
	        	    if(DeltaR2 < 0.5)
	            		vtx_Lep2_FSR = vtx_Lep2 + photon2;
		            else
		    	        vtx_Lep2_FSR = vtx_Lep2;
	            }
	            else
	    	        vtx_Lep2_FSR = vtx_Lep2;
	            	
	        }
	        else{ 
	        	vtx_Lep1_FSR = vtx_Lep1;
    	        vtx_Lep2_FSR = vtx_Lep2;
    	    }

                        
			vtx_eta_FSR1 = vtx_Lep1_FSR.Eta(); vtx_eta_FSR2 = vtx_Lep2_FSR.Eta();
            vtx_phi_FSR1 = vtx_Lep1_FSR.Phi();  vtx_phi_FSR2 = vtx_Lep2_FSR.Phi();
            vtx_pT_FSR1 = vtx_Lep1_FSR.Pt(); vtx_pT_FSR2 =  vtx_Lep2_FSR.Pt();
            vtx_m_FSR1 = vtx_Lep1_FSR.M(); vtx_m_FSR2 = vtx_Lep2_FSR.M();
            
            massZ_vtx_FSR = (vtx_Lep1_FSR + vtx_Lep2_FSR).M();
                        
            lep1p.SetPtEtaPhiM(double(vtx_pT_FSR1 + pterr1_VX), vtx_eta_FSR1, vtx_phi_FSR1, vtx_m_FSR1);
            lep2p.SetPtEtaPhiM(double(vtx_pT_FSR2 + pterr2_VX), vtx_eta_FSR2, vtx_phi_FSR2, vtx_m_FSR2);

            dm1 = (lep1p+vtx_Lep2_FSR).M()-(vtx_Lep1_FSR+vtx_Lep2_FSR).M();
            dm2 = (vtx_Lep1_FSR+lep2p).M()-(vtx_Lep1_FSR+vtx_Lep2_FSR).M();
 
            massZErr_vtx_FSR = TMath::Sqrt(dm1*dm1+dm2*dm2);
			/////// vtx constraint
*/

            if (lep_ecalDriven->size() > 0) {

               lep1_ecalDriven = (*lep_ecalDriven)[L1];
               lep2_ecalDriven = (*lep_ecalDriven)[L2];
 
              }
              
              
            /////// GEN info for MC
				LepMomId1 = (*lep_matchedR03_MomId)[L1];   LepMomId2 = (*lep_matchedR03_MomId)[L2];
				LepGrandMomId1 = (*lep_matchedR03_MomMomId)[L1];   LepGrandMomId2 = (*lep_matchedR03_MomMomId)[L2];

				genzm=0;
        GENmass2l=0;
				if(GENZ_mass->size()>0) genzm = (*GENZ_mass)[0];

				TLorentzVector GENlep1p, GENlep2p;
            
				if((*lep_genindex)[L1] >= 0 && (*lep_genindex)[L2] >= 0) {
          int genindex1 = (*lep_genindex)[L1];
				  int genindex2 = (*lep_genindex)[L2];

				  genLep_pt1=(*GENlep_pt)[genindex1];    genLep_pt2=(*GENlep_pt)[genindex2];
				  genLep_eta1=(*GENlep_eta)[genindex1];  genLep_eta2=(*GENlep_eta)[genindex2];
				  genLep_phi1=(*GENlep_phi)[genindex1];  genLep_phi2=(*GENlep_phi)[genindex2];
          genLep_mass1 = (*GENlep_mass)[genindex1]; genLep_mass2 = (*GENlep_mass)[genindex2];

				  GENlep1p.SetPtEtaPhiM(double(genLep_pt1),double(genLep_eta1),double(genLep_phi1),double(genLep_mass1));
				  GENlep2p.SetPtEtaPhiM(double(genLep_pt2),double(genLep_eta2),double(genLep_phi2),double(genLep_mass1));
				  GENmass2l = (GENlep1p+GENlep2p).M();
				}

            newtree->Fill();

        }
              
     cout<<"end reading tree"<<endl;                       
                                                               
     tmpFile->cd();                                        
                                                                  
     newtree->Write("passedEvents",TObject::kOverwrite);   
                                                                     
     tmpFile->Write();                                     
     tmpFile->Close();  
     
     std::cout<< "Number of events with passLepIndex.size()!=2: " << count_trig << std::endl;
     
}