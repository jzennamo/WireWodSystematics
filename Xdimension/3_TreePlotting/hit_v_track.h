//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Aug 28 09:53:06 2019 by ROOT version 6.12/06
// from TTree hit_v_track/HitPropertiesTrackProperties
// found on file: ../../outfile_mc.root
//////////////////////////////////////////////////////////

#ifndef hit_v_track_h
#define hit_v_track_h

#include <iostream>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <numeric>
#include <vector>
#include <algorithm>


// Header file for the classes stored in the TTree if any.

class hit_v_track {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Run;
   Int_t           Subrun;
   Int_t           Event;
   Int_t           Ntrks;
   Int_t           Nacpt;
   Int_t           n_acpt;
   Int_t           Nhits;
   Int_t           n_hits;
   Double_t        hit_Q;
   Double_t        hit_A;
   Double_t        hit_sigma;
   Double_t        hit_time;
   Int_t           hit_plane;
   Int_t           hit_isMC;
   Double_t        trk_x;
   Double_t        trk_y;
   Double_t        trk_z;
   Double_t        trk_L;
   Double_t        trk_StartX;
   Double_t        trk_StartY;
   Double_t        trk_StartZ;
   Double_t        trk_EndX;
   Double_t        trk_EndY;
   Double_t        trk_EndZ;
   Double_t        trk_cosTheta;
   Double_t        trk_phi;
   Double_t        trk_ThetaXZ;
   Double_t        trk_ThetaYZ;
   Double_t        trk_dEdx;
   Double_t        trk_fracMC;
   Double_t        clus_pand_width;
   Double_t        clus_pand_start_tick;
   Double_t        clus_pand_end_tick;
   Double_t        clus_patrec_width;
   Double_t        clus_patrec_start_tick;
   Double_t        clus_patrec_end_tick;
   Double_t        MC_optIntegral;
   Double_t        Data_optIntegral;
   Double_t        Overlay_optIntegral;

   // List of branches
   TBranch        *b_Run;   //!
   TBranch        *b_Subrun;   //!
   TBranch        *b_Event;   //!
   TBranch        *b_Ntrks;   //!
   TBranch        *b_Nacpt;   //!
   TBranch        *b_n_acpt;   //!
   TBranch        *b_Nhits;   //!
   TBranch        *b_n_hits;   //!
   TBranch        *b_hit_Q;   //!
   TBranch        *b_hit_A;   //!
   TBranch        *b_hit_sigma;   //!
   TBranch        *b_hit_time;   //!
   TBranch        *b_hit_plane;   //!
   TBranch        *b_hit_isMC;   //!
   TBranch        *b_trk_x;   //!
   TBranch        *b_trk_y;   //!
   TBranch        *b_trk_z;   //!
   TBranch        *b_trk_L;   //!
   TBranch        *b_trk_StartX;   //!
   TBranch        *b_trk_StartY;   //!
   TBranch        *b_trk_StartZ;   //!
   TBranch        *b_trk_EndX;   //!
   TBranch        *b_trk_EndY;   //!
   TBranch        *b_trk_EndZ;   //!
   TBranch        *b_trk_cosTheta;   //!
   TBranch        *b_trk_phi;   //!
   TBranch        *b_trk_ThetaXZ;   //!
   TBranch        *b_trk_ThetaYZ;   //!
   TBranch        *b_trk_dEdx;   //!
   TBranch        *b_trk_fracMC;   //!
   TBranch        *b_clus_pand_width;   //!
   TBranch        *b_clus_pand_start_tick;   //!
   TBranch        *b_clus_pand_end_tick;   //!
   TBranch        *b_clus_patrec_width;   //!
   TBranch        *b_clus_patrec_start_tick;   //!
   TBranch        *b_clus_patrec_end_tick;   //!
   TBranch        *b_MC_optIntegral;   //!
   TBranch        *b_Data_optIntegral;   //!
   TBranch        *b_Overlay_optIntegral;   //!

   hit_v_track(TTree *tree=0);
   virtual ~hit_v_track();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual double TruncMean(std::vector< double > v);
};

#endif

#ifdef hit_v_track_cxx
hit_v_track::hit_v_track(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/uboone/data/users/jaz8600/outfile_mc_CVhighstats.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/uboone/data/users/jaz8600/outfile_mc_CVhighstats.root");
      }
      f->GetObject("hit_v_track",tree);

   }
   Init(tree);
}

hit_v_track::~hit_v_track()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t hit_v_track::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t hit_v_track::LoadTree(Long64_t entry)
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

void hit_v_track::Init(TTree *tree)
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

   fChain->SetBranchAddress("Run", &Run, &b_Run);
   fChain->SetBranchAddress("Subrun", &Subrun, &b_Subrun);
   fChain->SetBranchAddress("Event", &Event, &b_Event);
   fChain->SetBranchAddress("Nhits", &Nhits, &b_Nhits);
   fChain->SetBranchAddress("hit_Q", &hit_Q, &b_hit_Q);
   fChain->SetBranchAddress("hit_A", &hit_A, &b_hit_A);
   fChain->SetBranchAddress("hit_sigma", &hit_sigma, &b_hit_sigma);
   fChain->SetBranchAddress("hit_time", &hit_time, &b_hit_time);
   fChain->SetBranchAddress("hit_plane", &hit_plane, &b_hit_plane);
   fChain->SetBranchAddress("hit_isMC", &hit_isMC, &b_hit_isMC);
   fChain->SetBranchAddress("trk_x", &trk_x, &b_trk_x);
   fChain->SetBranchAddress("trk_y", &trk_y, &b_trk_y);
   fChain->SetBranchAddress("trk_z", &trk_z, &b_trk_z);
   fChain->SetBranchAddress("trk_L", &trk_L, &b_trk_L);
   fChain->SetBranchAddress("trk_StartX", &trk_StartX, &b_trk_StartX);
   fChain->SetBranchAddress("trk_StartY", &trk_StartY, &b_trk_StartY);
   fChain->SetBranchAddress("trk_StartZ", &trk_StartZ, &b_trk_StartZ);
   fChain->SetBranchAddress("trk_EndX", &trk_EndX, &b_trk_EndX);
   fChain->SetBranchAddress("trk_EndY", &trk_EndY, &b_trk_EndY);
   fChain->SetBranchAddress("trk_EndZ", &trk_EndZ, &b_trk_EndZ);
   fChain->SetBranchAddress("trk_cosTheta", &trk_cosTheta, &b_trk_cosTheta);
   fChain->SetBranchAddress("trk_phi", &trk_phi, &b_trk_phi);
   fChain->SetBranchAddress("trk_ThetaXZ", &trk_ThetaXZ, &b_trk_ThetaXZ);
   fChain->SetBranchAddress("trk_ThetaYZ", &trk_ThetaYZ, &b_trk_ThetaYZ);
   fChain->SetBranchAddress("trk_dEdx", &trk_dEdx, &b_trk_dEdx);
   fChain->SetBranchAddress("trk_fracMC", &trk_fracMC, &b_trk_fracMC);
   Notify();
}

Bool_t hit_v_track::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

double hit_v_track::TruncMean(std::vector< double > v)
{

  std::vector< double > temp;
  temp = v;
  double tm;

  double temp_med = 0;

  for(int n = 0; n < 1000; n++){
    v = temp;
    temp.clear();
    std::sort(v.begin(),v.end());
    
    double sum = std::accumulate(v.begin(), v.end(), 0.0);
    double mean = sum / double(v.size());
    double sq_sum = std::inner_product(v.begin(), v.end(), v.begin(), 0.0);
    double stdev = std::sqrt((sq_sum / double(v.size())) - (mean * mean));
    double med = 0;
    int N = v.size();
    
    if (N % 2 != 0) 
      med = v[N/2]; 
    else
      med = (v[(N-1)/2] + v[N/2])/2.0; 
    
    for(int i = 0; i < v.size(); i++){
      
      if( v[i] > (med - 2*stdev) && v[i] < (med + 1.75*stdev) ){
	
	temp.push_back(v[i]);
      } // clip at stdev    
    }// check all values

    std::sort(temp.begin(),temp.end());

    
    double tsum = std::accumulate(temp.begin(), temp.end(), 0.0);
    double tmean = tsum / double(temp.size());

    if( fabs(tmean - temp_med) < 0.0001){
      std::cout << "BREAK at " << n << " med " << tm << " > temp_med " << temp_med << " from " << "[ " << temp[0] << " , " <<  temp[temp.size()-1] << " ]" << " out of " << temp.size()<< std::endl;
      break;
    }

    temp_med = tmean;
    

    tm = tmean;
    /*
    int Nt = temp.size(); 

    if (Nt % 2 != 0) 
      tm = temp[Nt/2]; 
    else
      tm = (temp[(Nt-1)/2] + temp[Nt/2])/2.0; 
    
    if( fabs(tm - temp_med) < 0.000001){
      std::cout << "BREAK at " << n << " med " << tm << " > temp_med " << temp_med << " from " << "[ " << temp[0] << " , " <<  temp[temp.size()-1] << " ]" << " out of " << temp.size()<< std::endl;
      break;
    }

    temp_med = tm;
    */
    
  }

  return tm;

}



void hit_v_track::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t hit_v_track::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef hit_v_track_cxx
