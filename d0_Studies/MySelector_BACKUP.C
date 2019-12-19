#define MySelector_cxx
// The class definition in MySelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.


// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("MySelector.C")
// root> T->Process("MySelector.C","some options")
// root> T->Process("MySelector.C+")
//


#include "MySelector.h"
#include <TH2.h>
#include <TStyle.h>

//----- Defining Functions -----//
void MySelector::SetNumBins(int n_bins) {nb = n_bins;}
void MySelector::SetNumBins_ZErr(int n_bins_ZErr) {nbZErr = n_bins_ZErr;}
void MySelector::SetRange_massZ(double low, double high) { massZ_lo = low; massZ_hi = high;}
void MySelector::SetRange_massZErr(double low, double high) { massZErr_lo = low; massZErr_hi = high;}
// void MySelector::SetLambda(int doLambda1_, double lambda1_, double lambda2_) {
// 
//      doLambda1 = doLambda1_; 
//      lambda1 = lambda1_; 
//      lambda2 = lambda2_;

void MySelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

   //----- User-defined stuff -----//
   rv_weight   = new RooRealVar("weight","weight", 0.00001, 10000);
   rv_massZ    = new RooRealVar("massZ","massZ", massZ_lo, massZ_hi);
   rv_massZErr = new RooRealVar("massZErr","massZErr", massZErr_lo, massZErr_hi);
   rastmp      = new RooArgSet(*rv_massZ, *rv_massZErr, *rv_weight);
   Data_Zlls   = new RooDataSet("Zlls","Zlls", *rastmp);

//   massZ_massZErr = new TH2F("massZ_massZErr", "", 300, massZ_lo, massZ_hi, 100, massZErr_lo, massZErr_hi);

//   rv_massZ->setBins(300,"fft");
//   rv_massZErr->setBins(50,"fft");

    rv_massZ->setBins(nb,"cache");
    rv_massZErr->setBins(nbZErr,"cache");
}

void MySelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t MySelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // When processing keyed objects with PROOF, the object is already loaded
   // and is available via the fObject pointer.
   //
   // This function should contain the \"body\" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

   fReader.SetEntry(entry);

      // I think rv stands for "return value"
   rv_massZ->setVal(*massZ);
   rv_weight->setVal(*weight);
   bool DEBUG  = false;

   if (DEBUG && (entry%50000==0)) {
      std::cout << "entry " << entry << ": massZ " << *massZ << "; weight " << *weight << endl;
   }

/*
   if (doLambda1 == 1) {
      // massZErr is the branch in the tree 
      rv_massZErr->setVal(*massZErr);
      if (DEBUG == 1 && (entry%50000==0)) {
      std::cout << "entry " << entry << ": massZErr " << *massZErr << endl;
      }
//      massZ_massZErr->Fill(*massZ, *massZErr);

      } else {

             TLorentzVector lep1 = TLorentzVector(0,0,0,0);
             TLorentzVector lep2 = TLorentzVector(0,0,0,0);
             TLorentzVector lep1p = TLorentzVector(0,0,0,0);
             TLorentzVector lep2p = TLorentzVector(0,0,0,0);
             
             // initialize 4-vectors with NO pT-err applied
             lep1.SetPtEtaPhiM(*pT1,*eta1,*phi1,*m1);
             lep2.SetPtEtaPhiM(*pT2,*eta2,*phi2,*m2);

             double nominalmZ = (lep1+lep2).M();

             // initialize 4-vectors WITH pT-err applied
             lep1p.SetPtEtaPhiM( (*pT1)+(*pterr1)*lambda1, *eta1, *phi1, *m1);// Initial lambdas are: lambda1=1  and lambda2=1
             lep2p.SetPtEtaPhiM( (*pT2)+(*pterr2)*lambda2, *eta2, *phi2, *m2);// Initial lambdas are: lambda1=1  and lambda2=1

             double dm1       = (lep1p+lep2).M() - nominalmZ; // dm1 tells us how much lep1pTerr affects the calculated mZ mass
             double dm2       = (lep1+lep2p).M() - nominalmZ; // dm2 tells us how much lep2pTerr affects the calculated mZ mass
//             double dm1 = (lep1p+lep2).M() - (lep1+lep2).M();
//             double dm2 = (lep1+lep2p).M() - (lep1+lep2).M();

             double new_massZErr = TMath::Sqrt(dm1*dm1+dm2*dm2);
             rv_massZErr->setVal(new_massZErr); // Update RooRealVar with a new massZErr val
             if (DEBUG == 1 && (entry%50000==0)) {
                std::cout << "entry " << entry << ": new_massZErr " << new_massZErr << endl;
        
//             massZ_massZErr->Fill(*massZ, new_massZErr);

             }
        }
*/
//   rv_massZ->setBins(40,"fft");
//   rv_massZErr->setBins(40,"fft");

   // Recall that rastmp is a RooArgSet(*rv_massZ, *rv_massZErr, *rv_weight).
   Data_Zlls->add(*rastmp);

   return kTRUE;
}

void MySelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void MySelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

}
