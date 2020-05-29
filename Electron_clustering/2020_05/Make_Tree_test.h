//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov 25 14:32:32 2019 by ROOT version 6.18/04
// from TTree tree1/Emu_Info
// found on file: Root_data_G4_sim_100eV_10k_0_1.root
//////////////////////////////////////////////////////////

#ifndef Make_Tree_test_h
#define Make_Tree_test_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class Make_Tree_test {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event;
   vector<int>     *vMC_TrackId;
   vector<int>     *vMC_StepId;
   vector<int>     *vMC_PartId;
   vector<int>     *vMC_ParentId;
   vector<int>     *vMC_CreatorId;
   vector<int>     *vMC_InteractionId;
   vector<int>     *vMC_VolumeId;
   vector<double>  *vMC_EnergyKin;
   vector<double>  *vMC_PX;
   vector<double>  *vMC_PY;
   vector<double>  *vMC_PZ;
   vector<int>     *vTrackId;
   vector<int>     *vStepId;
   vector<int>     *vPartId;
   vector<int>     *vParentId;
   vector<int>     *vCreatorId;
   vector<int>     *vInteractionId;
   vector<int>     *vVolumeId;
   vector<int>     *vCopyNo;
   vector<double>  *vEnergyKin;
   vector<double>  *vEnergyDep;
   vector<double>  *vPreHitX;
   vector<double>  *vPreHitY;
   vector<double>  *vPreHitZ;
   vector<double>  *vPostHitX;
   vector<double>  *vPostHitY;
   vector<double>  *vPostHitZ;
   vector<double>  *vHitPX;
   vector<double>  *vHitPY;
   vector<double>  *vHitPZ;
   vector<double>  *vTrackLength;
   vector<double>  *vStepLength;

   // List of branches
   TBranch        *b_Event;   //!
   TBranch        *b_vMC_TrackId;   //!
   TBranch        *b_vMC_StepId;   //!
   TBranch        *b_vMC_PartId;   //!
   TBranch        *b_vMC_ParentId;   //!
   TBranch        *b_vMC_CreatorId;   //!
   TBranch        *b_vMC_InteractionId;   //!
   TBranch        *b_vMC_VolumeId;   //!
   TBranch        *b_vMC_EnergyKin;   //!
   TBranch        *b_vMC_PX;   //!
   TBranch        *b_vMC_PY;   //!
   TBranch        *b_vMC_PZ;   //!
   TBranch        *b_vTrackId;   //!
   TBranch        *b_vStepId;   //!
   TBranch        *b_vPartId;   //!
   TBranch        *b_vParentId;   //!
   TBranch        *b_vCreatorId;   //!
   TBranch        *b_vInteractionId;   //!
   TBranch        *b_vVolumeId;   //!
   TBranch        *b_vCopyNo;   //!
   TBranch        *b_vEnergyKin;   //!
   TBranch        *b_vEnergyDep;   //!
   TBranch        *b_vPreHitX;   //!
   TBranch        *b_vPreHitY;   //!
   TBranch        *b_vPreHitZ;   //!
   TBranch        *b_vPostHitX;   //!
   TBranch        *b_vPostHitY;   //!
   TBranch        *b_vPostHitZ;   //!
   TBranch        *b_vHitPX;   //!
   TBranch        *b_vHitPY;   //!
   TBranch        *b_vHitPZ;   //!
   TBranch        *b_vTrackLength;   //!
   TBranch        *b_vStepLength;   //!

   Make_Tree_test(TTree *tree=0);
   virtual ~Make_Tree_test();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Make_Tree_test_cxx
Make_Tree_test::Make_Tree_test(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Root_data_G4_sim.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("Root_data_G4_sim.root");
      }
      f->GetObject("tree1",tree);

   }
   Init(tree);
}

Make_Tree_test::~Make_Tree_test()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Make_Tree_test::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Make_Tree_test::LoadTree(Long64_t entry)
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

void Make_Tree_test::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   vMC_TrackId = 0;
   vMC_StepId = 0;
   vMC_PartId = 0;
   vMC_ParentId = 0;
   vMC_CreatorId = 0;
   vMC_InteractionId = 0;
   vMC_VolumeId = 0;
   vMC_EnergyKin = 0;
   vMC_PX = 0;
   vMC_PY = 0;
   vMC_PZ = 0;
   vTrackId = 0;
   vStepId = 0;
   vPartId = 0;
   vParentId = 0;
   vCreatorId = 0;
   vInteractionId = 0;
   vVolumeId = 0;
   vCopyNo = 0;
   vEnergyKin = 0;
   vEnergyDep = 0;
   vPreHitX = 0;
   vPreHitY = 0;
   vPreHitZ = 0;
   vPostHitX = 0;
   vPostHitY = 0;
   vPostHitZ = 0;
   vHitPX = 0;
   vHitPY = 0;
   vHitPZ = 0;
   vTrackLength = 0;
   vStepLength = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("vMC_TrackId", &vMC_TrackId, &b_vMC_TrackId);
   fChain->SetBranchAddress("vMC_StepId", &vMC_StepId, &b_vMC_StepId);
   fChain->SetBranchAddress("vMC_PartId", &vMC_PartId, &b_vMC_PartId);
   fChain->SetBranchAddress("vMC_ParentId", &vMC_ParentId, &b_vMC_ParentId);
   fChain->SetBranchAddress("vMC_CreatorId", &vMC_CreatorId, &b_vMC_CreatorId);
   fChain->SetBranchAddress("vMC_InteractionId", &vMC_InteractionId, &b_vMC_InteractionId);
   fChain->SetBranchAddress("vMC_VolumeId", &vMC_VolumeId, &b_vMC_VolumeId);
   fChain->SetBranchAddress("vMC_EnergyKin", &vMC_EnergyKin, &b_vMC_EnergyKin);
   fChain->SetBranchAddress("vMC_PX", &vMC_PX, &b_vMC_PX);
   fChain->SetBranchAddress("vMC_PY", &vMC_PY, &b_vMC_PY);
   fChain->SetBranchAddress("vMC_PZ", &vMC_PZ, &b_vMC_PZ);
   fChain->SetBranchAddress("vTrackId", &vTrackId, &b_vTrackId);
   fChain->SetBranchAddress("vStepId", &vStepId, &b_vStepId);
   fChain->SetBranchAddress("vPartId", &vPartId, &b_vPartId);
   fChain->SetBranchAddress("vParentId", &vParentId, &b_vParentId);
   fChain->SetBranchAddress("vCreatorId", &vCreatorId, &b_vCreatorId);
   fChain->SetBranchAddress("vInteractionId", &vInteractionId, &b_vInteractionId);
   fChain->SetBranchAddress("vVolumeId", &vVolumeId, &b_vVolumeId);
   fChain->SetBranchAddress("vCopyNo", &vCopyNo, &b_vCopyNo);
   fChain->SetBranchAddress("vEnergyKin", &vEnergyKin, &b_vEnergyKin);
   fChain->SetBranchAddress("vEnergyDep", &vEnergyDep, &b_vEnergyDep);
   fChain->SetBranchAddress("vPreHitX", &vPreHitX, &b_vPreHitX);
   fChain->SetBranchAddress("vPreHitY", &vPreHitY, &b_vPreHitY);
   fChain->SetBranchAddress("vPreHitZ", &vPreHitZ, &b_vPreHitZ);
   fChain->SetBranchAddress("vPostHitX", &vPostHitX, &b_vPostHitX);
   fChain->SetBranchAddress("vPostHitY", &vPostHitY, &b_vPostHitY);
   fChain->SetBranchAddress("vPostHitZ", &vPostHitZ, &b_vPostHitZ);
   fChain->SetBranchAddress("vHitPX", &vHitPX, &b_vHitPX);
   fChain->SetBranchAddress("vHitPY", &vHitPY, &b_vHitPY);
   fChain->SetBranchAddress("vHitPZ", &vHitPZ, &b_vHitPZ);
   fChain->SetBranchAddress("vTrackLength", &vTrackLength, &b_vTrackLength);
   fChain->SetBranchAddress("vStepLength", &vStepLength, &b_vStepLength);
   Notify();
}

Bool_t Make_Tree_test::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Make_Tree_test::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Make_Tree_test::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Make_Tree_test_cxx
