//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Oct  6 14:26:56 2018 by ROOT version 6.15/01
// from TTree ntuple/data from ascii file
// found on file: SRIM_100K.root
//////////////////////////////////////////////////////////

#ifndef Make_SRIM_h
#define Make_SRIM_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class Make_SRIM {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         ev;
   Float_t         kene;
   Float_t         x;
   Float_t         y;
   Float_t         z;
   Float_t         rene;
   Float_t         dene;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_kene;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_rene;   //!
   TBranch        *b_dene;   //!

   Make_SRIM(TTree *tree=0);
   virtual ~Make_SRIM();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Make_SRIM_cxx
Make_SRIM::Make_SRIM(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SRIM_100K.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SRIM_100K.root");
      }
      f->GetObject("ntuple",tree);

   }
   Init(tree);
}

Make_SRIM::~Make_SRIM()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Make_SRIM::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Make_SRIM::LoadTree(Long64_t entry)
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

void Make_SRIM::Init(TTree *tree)
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

   fChain->SetBranchAddress("ev", &ev, &b_ev);
   fChain->SetBranchAddress("kene", &kene, &b_kene);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("rene", &rene, &b_rene);
   fChain->SetBranchAddress("dene", &dene, &b_dene);
   Notify();
}

Bool_t Make_SRIM::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Make_SRIM::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Make_SRIM::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Make_SRIM_cxx
