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

void MioSkim_Higgs(){

//  gSystem->Load("$ROOTSYS/test/libEvent");


    TString nome_file = "/eos/home-d/drosenzw/Samples/GluGluHToZZTo4L_M125_2018.root";
    TString nome_file_output = "/eos/home-d/drosenzw/Samples/LiteSkim/MC2018ggH_LiteSkim_test02";  // Don't include '.root'.

    std::vector<TString> ProdMode;
    ProdMode.clear();
    ProdMode.push_back("ggH");
//  ProdMode.push_back("VBF");
//  ProdMode.push_back("WplusH");
//  ProdMode.push_back("WminusH");
//  ProdMode.push_back("ZH");
//  ProdMode.push_back("ttH");
    
    std::vector<TString> ProdMode_File;
    ProdMode_File.clear();
//  ProdMode_File.push_back("GluGluHToZZTo4L");
//  ProdMode_File.push_back("VBF_HToZZTo4L");
//  ProdMode_File.push_back("WplusH_HToZZTo4L");
//  ProdMode_File.push_back("WminusH_HToZZTo4L");
//  ProdMode_File.push_back("ZH_HToZZ_4LFilter");
//  ProdMode_File.push_back("ttH_HToZZ_4LFilter");
    
    std::vector<TString> Background;
    Background.clear();    
/*    
    Background.push_back("GluGluToContinToZZTo4e");
    Background.push_back("GluGluToContinToZZTo4mu");
    Background.push_back("GluGluToContinToZZTo4tau");
    Background.push_back("GluGluToContinToZZTo2e2mu");
    Background.push_back("GluGluToContinToZZTo2e2tau");
    Background.push_back("GluGluToContinToZZTo2mu2tau");
    Background.push_back("ZZTo4L_powheg");
    Background.push_back("ZZTo4L_amcatnlo");
*/
//    Background.push_back("DYJetsToLL_M-10to50");
//    Background.push_back("DYJetsToLL_M-50");


    std::vector<int> massPoint;
    massPoint.clear();
 // massPoint.push_back(120);
 // massPoint.push_back(124);
    massPoint.push_back(125);
 // massPoint.push_back(126);
 // massPoint.push_back(130);
    
    std::vector<TString> year;
    year.clear();
/*    
    year.push_back("2016");
    year.push_back("2017");
*/
    year.push_back("2018");

    std::cout<<"OK\nsize of year vector = "<<year.size()<<std::endl;    
    
    for(int y = 0; y < year.size(); y++){
        std::cout<< "Running over year: " << year.at(y) <<std::endl;        
        for(int mP = 0; mP <  massPoint.size(); mP++){
            std::cout<<"Running over Higgs mass point:" << massPoint.at(mP)<<std::endl;
        for(int PM = 0; PM < ProdMode.size(); PM++){
                std::cout<< "Running over production mode:" << ProdMode.at(PM) <<std::endl;
/*
                TString nome_file = Form("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/%s/%s_M%i_%s", ProdMode[PM].Data(), ProdMode_File[PM].Data(), massPoint.at(mP), year[y].Data());
*/
                std::cout<< "Opening file:\n" << nome_file<<std::endl;
                TFile* infile = TFile::Open(nome_file);
                
                TTree* tree; 
                if(!infile){ 
                    std::cout<<"ERROR could not find the file"<<std::endl;
                    continue;
                }
                else{
                    std::cout<<"File trovato.OK"<<std::endl;
                    infile->cd("Ana");
                    tree = (TTree*)gDirectory->Get("passedEvents");
//                  Event *event   = new Event();
                    tree->SetBranchStatus("*",0);
//    tree->SetBranchStatus("Run",1);     
//    tree->SetBranchStatus("LumiSect",1);          
//    tree->SetBranchStatus("Event",1);              
                    tree->SetBranchStatus("GENMH",1);     
                    tree->SetBranchStatus("GENmass4l",1);          
                    tree->SetBranchStatus("GENmassZ1",1);              
                    tree->SetBranchStatus("GENmassZ2",1);     
                    tree->SetBranchStatus("GENlep_pt",1);          
                    tree->SetBranchStatus("GENlep_eta",1);              
                    tree->SetBranchStatus("GENlep_phi",1);     
                    tree->SetBranchStatus("GENlep_mass",1);          

                    tree->SetBranchStatus("nVtx",1);     
                    tree->SetBranchStatus("triggersPassed",1);  
                    tree->SetBranchStatus("passedFullSelection",1);     
            
                    tree->SetBranchStatus("passedFiducialSelection",1);          
                    tree->SetBranchStatus("passedZ4lSelection",1);                  
                    tree->SetBranchStatus("passedZXCRSelection",1);             
                    tree->SetBranchStatus("nZXCRFailedLeptons",1);                  
                    tree->SetBranchStatus("finalState",1);               
                    tree->SetBranchStatus("dataMCWeight",1);                 
                    tree->SetBranchStatus("k_qqZZ_qcd_M",1);                 
                    tree->SetBranchStatus("k_qqZZ_ewk",1);                
                    tree->SetBranchStatus("k_ggZZ",1);              
  
                    tree->SetBranchStatus("lep_pt",1);  
                    tree->SetBranchStatus("lep_pterr",1);                    
                    tree->SetBranchStatus("lep_id",1);                  
                    tree->SetBranchStatus("lep_eta",1);                 
                    tree->SetBranchStatus("lep_phi",1);                 
                    tree->SetBranchStatus("lep_mass",1);                
                    tree->SetBranchStatus("lep_ecalDriven", 1);
                    tree->SetBranchStatus("met",1);                           
                    tree->SetBranchStatus("mass4l",1);  
                    tree->SetBranchStatus("mass4lErr",1);                    
                    tree->SetBranchStatus("mass4lREFIT",1);                  
                    tree->SetBranchStatus("mass4lErrREFIT",1);                 
                    tree->SetBranchStatus("massZ1REFIT",1);                 
                    tree->SetBranchStatus("massZ1",1);                 
                    tree->SetBranchStatus("massZ2",1);                 
                    tree->SetBranchStatus("mass4mu",1);                
                    tree->SetBranchStatus("mass4e", 1);
                    tree->SetBranchStatus("mass2e2mu",1);                

                    tree->SetBranchStatus("EventCat",1);                
                    tree->SetBranchStatus("eventWeight", 1);
                    tree->SetBranchStatus("crossSection", 1);
                    tree->SetBranchStatus("genWeight", 1);
                    tree->SetBranchStatus("pileupWeight", 1);

                    tree->SetBranchStatus("D_bkg_kin",1);                
                    tree->SetBranchStatus("D_bkg", 1);
                    tree->SetBranchStatus("D_VBF", 1);


//     tree->Branch("me_qqZZ_MCFM", &me_qqZZ_MCFM, "me_qqZZ_MCFM/F");
//     tree->Branch("p0plus_m4l", &p0plus_m4l, "p0plus_m4l/F");
//     tree->Branch("bkg_m4l", &bkg_m4l, "bkg_m4l/F");
//     tree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/F");
//     tree->Branch("D_g4", &D_g4, "D_g4/F");
//     tree->Branch("D_g1g4", &D_g4, "D_g1g4/F");
//     tree->Branch("D_VBF1j",&D_VBF1j,"D_VBF1j/F");
//     tree->Branch("D_HadWH",&D_HadWH,"D_HadWH/F");
//     tree->Branch("D_HadZH",&D_HadZH,"D_HadZH/F");
//     tree->Branch("D_bkg_VBFdec",&D_bkg_VBFdec,"D_bkg_VBFdec/F");
//     tree->Branch("D_bkg_VHdec",&D_bkg_VHdec,"D_bkg_VHdec/F");
//     tree->Branch("D_VBF_QG",&D_VBF_QG,"D_VBF_QG/F");
//     tree->Branch("D_VBF1j_QG",&D_VBF1j_QG,"D_VBF1j_QG/F");
//     tree->Branch("D_HadWH_QG",&D_HadWH_QG,"D_HadWH_QG/F");
//     tree->Branch("D_HadZH_QG",&D_HadZH_QG,"D_HadZH_QG/F");
//     tree->Branch("me_0plus_JHU", &me_0plus_JHU, "me_0plus_JHU/F");

//   newtree->Branch("Dgg10_VAMCFM", &Dgg10_VAMCFM, "Dgg10_VAMCFM/D");
//   newtree->Branch("D_g4", &D_g4, "D_g4/D");
//   newtree->Branch("Djet_VAJHU", &Djet_VAJHU, "Djet_VAJHU/D");
//   newtree->Branch("D_VBF1j_VAJHU",&D_VBF1j_VAJHU,"D_VBF1j_VAJHU/D");
//   newtree->Branch("D_WHh_VAJHU",&D_WHh_VAJHU,"D_WHh_VAJHU/D");
//   newtree->Branch("D_ZHh_VAJHU",&D_ZHh_VAJHU,"D_ZHh_VAJHU/D");
//   newtree->Branch("D_VBF2j",&D_VBF2j,"D_VBF2j/D");
//   newtree->Branch("D_VBF1j",&D_VBF1j,"D_VBF1j/D");
//   newtree->Branch("D_WHh",&D_WHh,"D_WHh/D");
//   newtree->Branch("D_ZHh",&D_ZHh,"D_ZHh/D");

//                  tree->SetBranchAddress("pTL1", &pTL1);          
//                  tree->SetBranchAddress("etaL1", &etaL1);        
//                  tree->SetBranchAddress("phiL1", &phiL1);        
//                  tree->SetBranchAddress("mL1", &mL1);        
//                  tree->SetBranchAddress("idL1", &idL1);          
//                  tree->SetBranchAddress("pTErrL1", &pTErrL1);        
//                  tree->SetBranchAddress("pTL2", &pTL2);          
//                  tree->SetBranchAddress("etaL2", &etaL2);        
//                  tree->SetBranchAddress("phiL2", &phiL2);        
//                  tree->SetBranchAddress("mL2", &mL2);        
//                  tree->SetBranchAddress("idL2", &idL2);          
//                  tree->SetBranchAddress("pTErrL2", &pTErrL2);        
//                  tree->SetBranchAddress("pTL3", &pTL3);          
//                  tree->SetBranchAddress("etaL3", &etaL3);        
//                  tree->SetBranchAddress("phiL3", &phiL3);        
//                  tree->SetBranchAddress("mL3", &mL3);        
//                  tree->SetBranchAddress("idL3", &idL3);          
//                  tree->SetBranchAddress("pTErrL3", &pTErrL3);        
//                  tree->SetBranchAddress("pTL4", &pTL4);          
//                  tree->SetBranchAddress("etaL4", &etaL4);        
//                  tree->SetBranchAddress("phiL4", &phiL4);        
//                  tree->SetBranchAddress("mL4", &mL4);        
//                  tree->SetBranchAddress("idL4", &idL4);          
//                  tree->SetBranchAddress("pTErrL4", &pTErrL4); 
                
//   newtree->Branch("pTErrREFITL1",&pTErrREFITL1,"pTErrREFITL1/F");
//   newtree->Branch("pTErrREFITL2",&pTErrREFITL2,"pTErrREFITL2/F");
//   newtree->Branch("pTErrREFITL3",&pTErrREFITL3,"pTErrREFITL3/F");
//   newtree->Branch("pTErrREFITL4",&pTErrREFITL4,"pTErrREFITL4/F");
//   newtree->Branch("pTREFITL1",&pTREFITL1,"pTREFITL1/F");
//   newtree->Branch("pTREFITL2",&pTREFITL2,"pTREFITL2/F");
//   newtree->Branch("pTREFITL3",&pTREFITL3,"pTREFITL3/F");
//   newtree->Branch("pTREFITL4",&pTREFITL4,"pTREFITL4/F");
//   newtree->Branch("correlation",&correlation,"correlation/F");
//   newtree->Branch("nFSRPhotons", &nFSRPhotons, "nFSRPhotons/I");
//   newtree->Branch("mass4l_noFSR",&mass4l_noFSR,"mass4l_noFSR/F");
//   newtree->Branch("massZ1_noFSR",&massZ1_noFSR,"massZ1_noFSR/F");
//   newtree->Branch("massZ2_noFSR",&massZ2_noFSR,"massZ2_noFSR/F");
//   newtree->Branch("mass4l_up",&mass4l_up,"mass4l_up/F");
//   newtree->Branch("mass4l_dn",&mass4l_dn,"mass4l_dn/F");
//   newtree->Branch("mass4lErr_old",&mass4lErr_old,"mass4lErr_old/F");
//   newtree->Branch("pT4l",&pT4l,"pT4l/F");
//   newtree->Branch("massZ1",&massZ1,"massZ1/F");
//   newtree->Branch("massZ1Err",&massZ1Err,"massZ1Err/F");
//   newtree->Branch("massZ2",&massZ2,"massZ2/F"); 
//   newtree->Branch("njets_pt30_eta4p7",&njets_pt30_eta4p7,"njets_pt30_eta4p7/F");
//   newtree->Branch("pTj1",&pTj1,"pTj1/D");
//   newtree->Branch("etaj1",&etaj1,"etaj1/D");
//   newtree->Branch("pTj2",&pTj2,"pTj2/D");
//   newtree->Branch("etaj2",&etaj2,"etaj2/D");
                
/*                
                    TString nome_file_output = Form("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/%s/%s_M%i_%s", ProdMode[PM].Data(), ProdMode_File[PM].Data(), massPoint.at(mP), year[y].Data());
*/
                    TFile *newfile = new TFile(nome_file_output + "_skimmed.root","recreate");
                    TTree *newtree = tree->CloneTree();
                    newtree->Print();
                    newfile->Write();
                }
            }
        }

/*
        std::cout<<"Passo al backgorund"<<std::endl;
        for(int Bg = 0; Bg < Background.size(); Bg++){
            if(y != 0) continue;
           // if(Bg != 7 && Bg != 6) continue;
            TString nome_file;
            if(Bg < 6)
                nome_file = Form("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/ggZZ/%s_%s", Background[Bg].Data(), year[y].Data());
            else if(Bg == 6 || Bg == 7)
                nome_file = Form("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/qqZZ/%s_%s", Background[Bg].Data(), year[y].Data());
            else
                nome_file = Form("root://cmsio5.rc.ufl.edu//store/user/t2/users/ferrico/Full_RunII/ZplusJets/%s_%s", Background[Bg].Data(), year[y].Data());
                std::cout<<nome_file<<std::endl;
                TFile* infile = TFile::Open(nome_file + ".root");
                
                TTree* tree; 
                if(!infile){ 
                    std::cout<<"ERROR could not find the file"<<std::endl;
                    continue;
                }
                else{
                    std::cout<<"File trovato.OK"<<std::endl;
                    infile->cd("Ana");
                    tree = (TTree*)gDirectory->Get("passedEvents");

                    tree->SetBranchStatus("*",0);

                    tree->SetBranchStatus("GENMH",1);     
                    tree->SetBranchStatus("GENmass4l",1);          
                    tree->SetBranchStatus("GENmassZ1",1);              
                    tree->SetBranchStatus("GENmassZ2",1);     
                    tree->SetBranchStatus("GENlep_pt",1);          
                    tree->SetBranchStatus("GENlep_eta",1);              
                    tree->SetBranchStatus("GENlep_phi",1);     
                    tree->SetBranchStatus("GENlep_mass",1);          

                    tree->SetBranchStatus("nVtx",1);     
                    tree->SetBranchStatus("triggersPassed",1);  
                    tree->SetBranchStatus("passedFullSelection",1);     
            
                    tree->SetBranchStatus("passedFiducialSelection",1);          
                    tree->SetBranchStatus("passedZ4lSelection",1);                  
                    tree->SetBranchStatus("passedZXCRSelection",1);             
                    tree->SetBranchStatus("nZXCRFailedLeptons",1);                  
                    tree->SetBranchStatus("finalState",1);               
                    tree->SetBranchStatus("dataMCWeight",1);                 
                    tree->SetBranchStatus("k_qqZZ_qcd_M",1);                 
                    tree->SetBranchStatus("k_qqZZ_ewk",1);                
                    tree->SetBranchStatus("k_ggZZ",1);              
  
                    tree->SetBranchStatus("lep_pt",1);  
                    tree->SetBranchStatus("lep_pterr",1);                    
                    tree->SetBranchStatus("lep_id",1);                  
                    tree->SetBranchStatus("lep_eta",1);                 
                    tree->SetBranchStatus("lep_phi",1);                 
                    tree->SetBranchStatus("lep_mass",1);                
                    tree->SetBranchStatus("lep_ecalDriven", 1);
                    tree->SetBranchStatus("met",1);                           
                    tree->SetBranchStatus("mass4l",1);  
                    tree->SetBranchStatus("mass4lErr",1);                    
                    tree->SetBranchStatus("mass4lREFIT",1);                  
                    tree->SetBranchStatus("mass4lErrREFIT",1);                 
                    tree->SetBranchStatus("massZ1REFIT",1);                 
                    tree->SetBranchStatus("mass4mu",1);                
                    tree->SetBranchStatus("mass4e", 1);
                    tree->SetBranchStatus("mass2e2mu",1);                

                    tree->SetBranchStatus("EventCat",1);                
                    tree->SetBranchStatus("eventWeight", 1);
                    tree->SetBranchStatus("crossSection", 1);
                    tree->SetBranchStatus("genWeight", 1);
                    tree->SetBranchStatus("pileupWeight", 1);
                
                    TString nome_file_output;
                if(Bg < 6)
                    nome_file_output = Form("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/ggZZ/%s_%s", Background[Bg].Data(), year[y].Data());
                else if(Bg == 6 || Bg == 7)
                    nome_file_output = Form("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/qqZZ/%s_%s", Background[Bg].Data(), year[y].Data());
                else
                    nome_file_output = Form("/raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/Full_RunII/ZplusJets/%s_%s", Background[Bg].Data(), year[y].Data());


                    TFile *newfile = new TFile(nome_file_output + "_skimmed.root","recreate");
                    TTree *newtree = tree->CloneTree();
                    newtree->Print();
                    newfile->Write();
                }
            }
*/






    }
                

    

}

