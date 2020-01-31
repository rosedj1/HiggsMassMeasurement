//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec 11 18:58:47 2019 by ROOT version 6.12/07
// from TTree passedEvents/passedEvents
// found on file: /raid/raid8/ferrico/Useful_Code_HZZ/CMSSW_10_2_15/src/provad0_2mu.root
//////////////////////////////////////////////////////////

#ifndef MySelector_h
#define MySelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// Headers needed by this particular selector


class MySelector : public TSelector {
public :
   TTreeReader     fReader;  //!the tree reader
   TTree          *fChain = 0;   //!pointer to the analyzed TTree or TChain

   // Readers to access the data (delete the ones you do not need).
   TTreeReaderValue<Double_t> massZ = {fReader, "massZ"};
   TTreeReaderValue<Double_t> massZErr = {fReader, "massZErr"};
   TTreeReaderValue<Double_t> massZErrOld = {fReader, "massZErrOld"};
   TTreeReaderValue<Double_t> pT1 = {fReader, "pT1"};
   TTreeReaderValue<Double_t> pT2 = {fReader, "pT2"};
   TTreeReaderValue<Double_t> eta1 = {fReader, "eta1"};
   TTreeReaderValue<Double_t> eta2 = {fReader, "eta2"};
   TTreeReaderValue<Double_t> phi1 = {fReader, "phi1"};
   TTreeReaderValue<Double_t> phi2 = {fReader, "phi2"};
   TTreeReaderValue<Double_t> m1 = {fReader, "m1"};
   TTreeReaderValue<Double_t> m2 = {fReader, "m2"};
   TTreeReaderValue<Double_t> d0BS1 = {fReader, "d0BS1"};
   TTreeReaderValue<Double_t> d0BS2 = {fReader, "d0BS2"};
   TTreeReaderValue<Double_t> d0PV1 = {fReader, "d0PV1"};
   TTreeReaderValue<Double_t> d0PV2 = {fReader, "d0PV2"};
   TTreeReaderValue<Double_t> Iso1 = {fReader, "Iso1"};
   TTreeReaderValue<Double_t> Iso2 = {fReader, "Iso2"};
   TTreeReaderValue<Int_t> Id1 = {fReader, "Id1"};
   TTreeReaderValue<Int_t> Id2 = {fReader, "Id2"};
   TTreeReaderValue<Int_t> Tight1 = {fReader, "Tight1"};
   TTreeReaderValue<Int_t> Tight2 = {fReader, "Tight2"};
   TTreeReaderValue<Double_t> pterr1 = {fReader, "pterr1"};
   TTreeReaderValue<Double_t> pterr2 = {fReader, "pterr2"};
   TTreeReaderValue<Double_t> pterr1old = {fReader, "pterr1old"};
   TTreeReaderValue<Double_t> pterr2old = {fReader, "pterr2old"};
   TTreeReaderValue<Double_t> Met = {fReader, "Met"};
   TTreeReaderValue<Double_t> weight = {fReader, "weight"};
   TTreeReaderValue<Double_t> genzm = {fReader, "genzm"};
   TTreeReaderValue<Double_t> GENmass2l = {fReader, "GENmass2l"};
   TTreeReaderValue<Double_t> genLep_pt1 = {fReader, "genLep_pt1"};
   TTreeReaderValue<Double_t> genLep_pt2 = {fReader, "genLep_pt2"};
   TTreeReaderValue<Double_t> genLep_eta1 = {fReader, "genLep_eta1"};
   TTreeReaderValue<Double_t> genLep_eta2 = {fReader, "genLep_eta2"};
   TTreeReaderValue<Double_t> genLep_phi1 = {fReader, "genLep_phi1"};
   TTreeReaderValue<Double_t> genLep_phi2 = {fReader, "genLep_phi2"};
   TTreeReaderValue<Int_t> nFSRPhotons = {fReader, "nFSRPhotons"};
   TTreeReaderValue<Int_t> lep1_ecalDriven = {fReader, "lep1_ecalDriven"};
   TTreeReaderValue<Int_t> lep2_ecalDriven = {fReader, "lep2_ecalDriven"};


   MySelector(TTree * /*tree*/ =0) { }
   virtual ~MySelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   ClassDef(MySelector,0);

};

#endif

#ifdef MySelector_cxx
void MySelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the reader is initialized.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   fReader.SetTree(tree);
}

Bool_t MySelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}


#endif // #ifdef MySelector_cxx
