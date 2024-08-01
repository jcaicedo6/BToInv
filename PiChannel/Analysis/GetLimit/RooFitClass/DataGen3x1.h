//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon May 18 12:14:36 2020 by ROOT version 6.20/04
// from TTree tree/This is a tree
// found on file: pairPYTHIA8_3x1_nolepton_r3_1000000_pID_15.root
//////////////////////////////////////////////////////////

#ifndef DataGen3x1_h
#define DataGen3x1_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DataGen3x1 {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           nchparts;
   Double_t        px[4];   //[nchparts]
   Double_t        py[4];   //[nchparts]
   Double_t        pz[4];   //[nchparts]
   Double_t        E[4];   //[nchparts]
   Int_t           pdgID[4];   //[nchparts]
   Int_t           pdgIDMom[4];   //[nchparts]
   Char_t          IsSignal;
   Int_t           parent;
   Int_t           ndaughtersTauN;
   Int_t           ndaughtersTauP;
   Double_t        daughterTauN[10];   //[ndaughtersTauN]
   Double_t        daughterTauP[9];   //[ndaughtersTauP]

   // List of branches
   TBranch        *b_nchparts;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_pdgID;   //!
   TBranch        *b_pdgIDMom;   //!
   TBranch        *b_IsSignal;   //!
   TBranch        *b_parent;   //!
   TBranch        *b_ndaughtersTauN;   //!
   TBranch        *b_ndaughtersTauP;   //!
   TBranch        *b_daughterTauN;   //!
   TBranch        *b_daughterTauP;   //!

   DataGen3x1(TTree *tree=0);
   virtual ~DataGen3x1();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
  //virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

//#ifdef DataGen3x1_cxx
DataGen3x1::DataGen3x1(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pairPYTHIA8_3x1_nolepton_r3_1000000_pID_15.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pairPYTHIA8_3x1_nolepton_r3_1000000_pID_15.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

DataGen3x1::~DataGen3x1()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DataGen3x1::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DataGen3x1::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void DataGen3x1::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("nchparts", &nchparts, &b_nchparts);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("pdgID", pdgID, &b_pdgID);
   fChain->SetBranchAddress("pdgIDMom", pdgIDMom, &b_pdgIDMom);
   fChain->SetBranchAddress("IsSignal", &IsSignal, &b_IsSignal);
   fChain->SetBranchAddress("parent", &parent, &b_parent);
   fChain->SetBranchAddress("ndaughtersTauN", &ndaughtersTauN, &b_ndaughtersTauN);
   fChain->SetBranchAddress("ndaughtersTauP", &ndaughtersTauP, &b_ndaughtersTauP);
   fChain->SetBranchAddress("daughterTauN", daughterTauN, &b_daughterTauN);
   fChain->SetBranchAddress("daughterTauP", daughterTauP, &b_daughterTauP);
   Notify();
}

Bool_t DataGen3x1::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DataGen3x1::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DataGen3x1::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
//#endif // #ifdef DataGen3x1_cxx
