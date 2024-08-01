//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep  3 22:47:42 2019 by ROOT version 6.14/06
// from TTree tree/
// found on file: myoutput1CMS.root
//////////////////////////////////////////////////////////

#ifndef DataGen_h
#define DataGen_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class DataGen {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMax__event_ = 1;
   static constexpr Int_t kMax__run_ = 1;
   static constexpr Int_t kMax__experiment_ = 1;
   static constexpr Int_t kMax__ncandidates_ = 1;
   static constexpr Int_t kMax__weight_ = 1;

   // Declaration of leaf types
   Int_t           __event__;
   Int_t           __run__;
   Int_t           __experiment__;
   Int_t           __ncandidates__;
   Float_t         __weight__;
   Double_t        M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        ErrM[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        SigM[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        InvM[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        px[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        py[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        pz[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        pt[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        p[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        E[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_0_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_0_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_0_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_0_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_0_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_0_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_1_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_1_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_1_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_1_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_1_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_1_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_e_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_e_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_e_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_e_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_e_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_e_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_invi_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_invi_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_invi_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_invi_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_invi_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_invi_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_0_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_0_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_0_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_0_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_0_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_0_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_1_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_1_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_1_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_1_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_1_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_1_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_2_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_2_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_2_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_2_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_2_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_pi_2_M[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_nu_tau_px_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_nu_tau_py_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_nu_tau_pz_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_nu_tau_p_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_nu_tau_E_CMS[kMax__ncandidates_];   //[__ncandidates__]
   Double_t        tau_nu_tau_M[kMax__ncandidates_];   //[__ncandidates__]

   // List of branches
   TBranch        *b___event__;   //!
   TBranch        *b___run__;   //!
   TBranch        *b___experiment__;   //!
   TBranch        *b___ncandidates__;   //!
   TBranch        *b___weight__;   //!
   TBranch        *b_M;   //!
   TBranch        *b_ErrM;   //!
   TBranch        *b_SigM;   //!
   TBranch        *b_InvM;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_pt;   //!
   TBranch        *b_p;   //!
   TBranch        *b_E;   //!
   TBranch        *b_tau_0_px_CMS;   //!
   TBranch        *b_tau_0_py_CMS;   //!
   TBranch        *b_tau_0_pz_CMS;   //!
   TBranch        *b_tau_0_p_CMS;   //!
   TBranch        *b_tau_0_E_CMS;   //!
   TBranch        *b_tau_0_M;   //!
   TBranch        *b_tau_1_px_CMS;   //!
   TBranch        *b_tau_1_py_CMS;   //!
   TBranch        *b_tau_1_pz_CMS;   //!
   TBranch        *b_tau_1_p_CMS;   //!
   TBranch        *b_tau_1_E_CMS;   //!
   TBranch        *b_tau_1_M;   //!
   TBranch        *b_tau_e_px_CMS;   //!
   TBranch        *b_tau_e_py_CMS;   //!
   TBranch        *b_tau_e_pz_CMS;   //!
   TBranch        *b_tau_e_p_CMS;   //!
   TBranch        *b_tau_e_E_CMS;   //!
   TBranch        *b_tau_e_M;   //!
   TBranch        *b_tau_invi_px_CMS;   //!
   TBranch        *b_tau_invi_py_CMS;   //!
   TBranch        *b_tau_invi_pz_CMS;   //!
   TBranch        *b_tau_invi_p_CMS;   //!
   TBranch        *b_tau_invi_E_CMS;   //!
   TBranch        *b_tau_invi_M;   //!
   TBranch        *b_tau_pi_0_px_CMS;   //!
   TBranch        *b_tau_pi_0_py_CMS;   //!
   TBranch        *b_tau_pi_0_pz_CMS;   //!
   TBranch        *b_tau_pi_0_p_CMS;   //!
   TBranch        *b_tau_pi_0_E_CMS;   //!
   TBranch        *b_tau_pi_0_M;   //!
   TBranch        *b_tau_pi_1_px_CMS;   //!
   TBranch        *b_tau_pi_1_py_CMS;   //!
   TBranch        *b_tau_pi_1_pz_CMS;   //!
   TBranch        *b_tau_pi_1_p_CMS;   //!
   TBranch        *b_tau_pi_1_E_CMS;   //!
   TBranch        *b_tau_pi_1_M;   //!
   TBranch        *b_tau_pi_2_px_CMS;   //!
   TBranch        *b_tau_pi_2_py_CMS;   //!
   TBranch        *b_tau_pi_2_pz_CMS;   //!
   TBranch        *b_tau_pi_2_p_CMS;   //!
   TBranch        *b_tau_pi_2_E_CMS;   //!
   TBranch        *b_tau_pi_2_M;   //!
   TBranch        *b_tau_nu_tau_px_CMS;   //!
   TBranch        *b_tau_nu_tau_py_CMS;   //!
   TBranch        *b_tau_nu_tau_pz_CMS;   //!
   TBranch        *b_tau_nu_tau_p_CMS;   //!
   TBranch        *b_tau_nu_tau_E_CMS;   //!
   TBranch        *b_tau_nu_tau_M;   //!

   DataGen(TTree *tree=0);
   virtual ~DataGen();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef DataGen_cxx
DataGen::DataGen(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("myoutput1CMS.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("myoutput1CMS.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

DataGen::~DataGen()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t DataGen::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t DataGen::LoadTree(Long64_t entry)
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

void DataGen::Init(TTree *tree)
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

   fChain->SetBranchAddress("__event__", &__event__, &b___event__);
   fChain->SetBranchAddress("__run__", &__run__, &b___run__);
   fChain->SetBranchAddress("__experiment__", &__experiment__, &b___experiment__);
   fChain->SetBranchAddress("__ncandidates__", &__ncandidates__, &b___ncandidates__);
   fChain->SetBranchAddress("__weight__", &__weight__, &b___weight__);
   fChain->SetBranchAddress("M", M, &b_M);
   fChain->SetBranchAddress("ErrM", ErrM, &b_ErrM);
   fChain->SetBranchAddress("SigM", SigM, &b_SigM);
   fChain->SetBranchAddress("InvM", InvM, &b_InvM);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("pt", pt, &b_pt);
   fChain->SetBranchAddress("p", p, &b_p);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("tau_0_px_CMS", tau_0_px_CMS, &b_tau_0_px_CMS);
   fChain->SetBranchAddress("tau_0_py_CMS", tau_0_py_CMS, &b_tau_0_py_CMS);
   fChain->SetBranchAddress("tau_0_pz_CMS", tau_0_pz_CMS, &b_tau_0_pz_CMS);
   fChain->SetBranchAddress("tau_0_p_CMS", tau_0_p_CMS, &b_tau_0_p_CMS);
   fChain->SetBranchAddress("tau_0_E_CMS", tau_0_E_CMS, &b_tau_0_E_CMS);
   fChain->SetBranchAddress("tau_0_M", tau_0_M, &b_tau_0_M);
   fChain->SetBranchAddress("tau_1_px_CMS", tau_1_px_CMS, &b_tau_1_px_CMS);
   fChain->SetBranchAddress("tau_1_py_CMS", tau_1_py_CMS, &b_tau_1_py_CMS);
   fChain->SetBranchAddress("tau_1_pz_CMS", tau_1_pz_CMS, &b_tau_1_pz_CMS);
   fChain->SetBranchAddress("tau_1_p_CMS", tau_1_p_CMS, &b_tau_1_p_CMS);
   fChain->SetBranchAddress("tau_1_E_CMS", tau_1_E_CMS, &b_tau_1_E_CMS);
   fChain->SetBranchAddress("tau_1_M", tau_1_M, &b_tau_1_M);
   fChain->SetBranchAddress("tau_e_px_CMS", tau_e_px_CMS, &b_tau_e_px_CMS);
   fChain->SetBranchAddress("tau_e_py_CMS", tau_e_py_CMS, &b_tau_e_py_CMS);
   fChain->SetBranchAddress("tau_e_pz_CMS", tau_e_pz_CMS, &b_tau_e_pz_CMS);
   fChain->SetBranchAddress("tau_e_p_CMS", tau_e_p_CMS, &b_tau_e_p_CMS);
   fChain->SetBranchAddress("tau_e_E_CMS", tau_e_E_CMS, &b_tau_e_E_CMS);
   fChain->SetBranchAddress("tau_e_M", tau_e_M, &b_tau_e_M);
   fChain->SetBranchAddress("tau_invi_px_CMS", tau_invi_px_CMS, &b_tau_invi_px_CMS);
   fChain->SetBranchAddress("tau_invi_py_CMS", tau_invi_py_CMS, &b_tau_invi_py_CMS);
   fChain->SetBranchAddress("tau_invi_pz_CMS", tau_invi_pz_CMS, &b_tau_invi_pz_CMS);
   fChain->SetBranchAddress("tau_invi_p_CMS", tau_invi_p_CMS, &b_tau_invi_p_CMS);
   fChain->SetBranchAddress("tau_invi_E_CMS", tau_invi_E_CMS, &b_tau_invi_E_CMS);
   fChain->SetBranchAddress("tau_invi_M", tau_invi_M, &b_tau_invi_M);
   fChain->SetBranchAddress("tau_pi_0_px_CMS", tau_pi_0_px_CMS, &b_tau_pi_0_px_CMS);
   fChain->SetBranchAddress("tau_pi_0_py_CMS", tau_pi_0_py_CMS, &b_tau_pi_0_py_CMS);
   fChain->SetBranchAddress("tau_pi_0_pz_CMS", tau_pi_0_pz_CMS, &b_tau_pi_0_pz_CMS);
   fChain->SetBranchAddress("tau_pi_0_p_CMS", tau_pi_0_p_CMS, &b_tau_pi_0_p_CMS);
   fChain->SetBranchAddress("tau_pi_0_E_CMS", tau_pi_0_E_CMS, &b_tau_pi_0_E_CMS);
   fChain->SetBranchAddress("tau_pi_0_M", tau_pi_0_M, &b_tau_pi_0_M);
   fChain->SetBranchAddress("tau_pi_1_px_CMS", tau_pi_1_px_CMS, &b_tau_pi_1_px_CMS);
   fChain->SetBranchAddress("tau_pi_1_py_CMS", tau_pi_1_py_CMS, &b_tau_pi_1_py_CMS);
   fChain->SetBranchAddress("tau_pi_1_pz_CMS", tau_pi_1_pz_CMS, &b_tau_pi_1_pz_CMS);
   fChain->SetBranchAddress("tau_pi_1_p_CMS", tau_pi_1_p_CMS, &b_tau_pi_1_p_CMS);
   fChain->SetBranchAddress("tau_pi_1_E_CMS", tau_pi_1_E_CMS, &b_tau_pi_1_E_CMS);
   fChain->SetBranchAddress("tau_pi_1_M", tau_pi_1_M, &b_tau_pi_1_M);
   fChain->SetBranchAddress("tau_pi_2_px_CMS", tau_pi_2_px_CMS, &b_tau_pi_2_px_CMS);
   fChain->SetBranchAddress("tau_pi_2_py_CMS", tau_pi_2_py_CMS, &b_tau_pi_2_py_CMS);
   fChain->SetBranchAddress("tau_pi_2_pz_CMS", tau_pi_2_pz_CMS, &b_tau_pi_2_pz_CMS);
   fChain->SetBranchAddress("tau_pi_2_p_CMS", tau_pi_2_p_CMS, &b_tau_pi_2_p_CMS);
   fChain->SetBranchAddress("tau_pi_2_E_CMS", tau_pi_2_E_CMS, &b_tau_pi_2_E_CMS);
   fChain->SetBranchAddress("tau_pi_2_M", tau_pi_2_M, &b_tau_pi_2_M);
   fChain->SetBranchAddress("tau_nu_tau_px_CMS", tau_nu_tau_px_CMS, &b_tau_nu_tau_px_CMS);
   fChain->SetBranchAddress("tau_nu_tau_py_CMS", tau_nu_tau_py_CMS, &b_tau_nu_tau_py_CMS);
   fChain->SetBranchAddress("tau_nu_tau_pz_CMS", tau_nu_tau_pz_CMS, &b_tau_nu_tau_pz_CMS);
   fChain->SetBranchAddress("tau_nu_tau_p_CMS", tau_nu_tau_p_CMS, &b_tau_nu_tau_p_CMS);
   fChain->SetBranchAddress("tau_nu_tau_E_CMS", tau_nu_tau_E_CMS, &b_tau_nu_tau_E_CMS);
   fChain->SetBranchAddress("tau_nu_tau_M", tau_nu_tau_M, &b_tau_nu_tau_M);
   Notify();
}

Bool_t DataGen::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DataGen::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t DataGen::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef DataGen_cxx
