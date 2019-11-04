//----------------------//
// ---- HISTOGRAMS ---- //
//----------------------//

#if !defined (__CINT__) || defined (__MAKECINT__)
#include "Rtypes.h"
#endif

#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TTree.h"
#include "TFile.h"
#include "TStyle.h"
#include <vector>
#include <iostream>

//Float_t pi = TMath::Pi();
const double pi=TMath::Pi();

Double_t first_x_pos=0;
Double_t first_y_pos=0;
Double_t first_z_pos=0;
Double_t last_x_pos=0;
Double_t last_y_pos=0;
Double_t last_z_pos=0;
Double_t trk_rec[3]={};

TH1F *h_Ene_MC = new TH1F("h_Ene_MC","; E_{MC} (MeV)",10000,0,10000); 
TH1F *h_Phi_MC = new TH1F("h_Phi_MC","; #phi_{MC} (rad)",100,-pi,pi);
TH1F *h_CosTheta_MC = new TH1F("h_CosTheta_MC","; cos#theta_{MC}",100,-1,-1);
TH1F *h_Theta_MC = new TH1F("h_Theta_MC","; #theta_{MC}",100,0,pi);

TH1F *h_PartName_EmuTop = new TH1F("h_PartName_EmuTop","; PartName",100,0,10);
TH1F *h_PartName_EmuDown = new TH1F("h_PartName_EmuDown","; PartName",100,0,10);
TH1F *h_num_trk_Emu = new TH1F("h_num_trk_Emu","; counts",100,0,10);

TH1F *h_trk_energy_Top = new TH1F("h_trk_energy_Top","; E_{MC} (keV)",100,0,100);
TH1F *h_trk_len_Top = new TH1F("h_trk_len_Top","; length(#mum)",100,0,10);
TH1F *h_trk_phi_Top = new TH1F("h_trk_phi_Top","; rad",100,-pi,pi);
TH1F *h_trk_theta_Top = new TH1F("h_trk_theta_Top","; rad",100,-pi,pi);

TH1F *h_trk_energy_Down = new TH1F("h_trk_energy_Down","; E_{MC} (keV)",100,0,100);
TH1F *h_trk_len_Down = new TH1F("h_trk_len_Down","; length(#mum)",100,0,10);
TH1F *h_trk_phi_Down = new TH1F("h_trk_phi_Down","; rad",100,-pi,pi);
TH1F *h_trk_theta_Down = new TH1F("h_trk_theta_Down","; rad",100,-pi,pi);

TH1F *h_ene_all_Top = new TH1F("h_ene_all_Top","; E_{MC} (keV)",100,0,100);
TH1F *h_len_all_Top = new TH1F("h_len_all_Top","; length(#mum)",100,0,10);
TH1F *h_phi_all_Top = new TH1F("h_phi_all_Top","; rad",100,-pi,pi);
TH1F *h_theta_all_Top = new TH1F("h_theta_all_Top","; rad",100,-pi,pi);

TH1F *h_ene_all_Down = new TH1F("h_ene_all_Down","; E_{MC} (keV)",100,0,100);
TH1F *h_len_all_Down = new TH1F("h_len_all_Down","; length(#mum)",100,0,10);
TH1F *h_phi_all_Down = new TH1F("h_phi_all_Down","; rad",100,-pi,pi);
TH1F *h_theta_all_Down = new TH1F("h_theta_all_Down","; rad",100,-pi,pi);



Int_t tmp_Event_Top;
Int_t tmp_Event_Down;
Int_t tmp_trk_Top=0;
Int_t tmp_trk_Down=0;
Int_t tmp_in_trk=0;
Int_t tmp_out_trk=0;
Int_t tmp_in_Event=0;
Int_t tmp_out_Event=0;
Int_t tmp_pdg_Top=0;
Int_t tmp_pdg_Down=0;
Int_t whichParticleTop=0;
Int_t whichParticleDown=0;
Int_t ntrack_contained_Top=0;
Int_t ntrack_contained_Down=0;
Double_t first_energy_step_Top=0;
Double_t first_energy_step_Down=0;
Double_t last_energy_step_Top=0;
Double_t last_energy_step_Down=0;
Double_t tot_ene_dep_Top=0;
Double_t tot_ene_dep_Down=0;
Double_t first_x_pos_Top=0;
Double_t first_y_pos_Top=0;
Double_t first_z_pos_Top=0;
Double_t last_x_pos_Top=0;
Double_t last_y_pos_Top=0;
Double_t last_z_pos_Top=0;
Double_t first_x_pos_Down=0;
Double_t first_y_pos_Down=0;
Double_t first_z_pos_Down=0;
Double_t last_x_pos_Down=0;
Double_t last_y_pos_Down=0;
Double_t last_z_pos_Down=0;
Bool_t start_Top=false;
Bool_t stop_Top=false;
Bool_t pass_Top=false;
Bool_t pass_Down=false;
Bool_t start_Down=false;
Bool_t stop_Down=false;
Bool_t out_trk=false;
Bool_t particle_in=false;
Bool_t particle_out=false;
Int_t ntrackEmuTop=0;
Int_t ntrackEmuDown=0;

double MC_angles[2]={};
double trk_rec_top[3]={};
double trk_rec_down[3]={};
double track_length_Top=0;
double phi_Top=0;
double theta_Top=0;
double track_length_Down=0;
double phi_Down=0;
double theta_Down=0;

double track_length_all_Top=0;
double phi_all_Top=0;
double theta_all_Top=0;
double track_length_all_Down=0;
double phi_all_Down=0;
double theta_all_Down=0;


Double_t hit_x_pre_last=0;
Double_t hit_y_pre_last=0;
Double_t hit_z_pre_last=0;

Bool_t pre_last=false;
Double_t step_len=0.0;
Double_t trk_poly_len=0.0;
Double_t tot_trk_poly_len=0.0;
Double_t dE_dx=0.;
Double_t first_ene=0;
Double_t sum_step_len=0.;

Int_t  is_a_father[10000];
Double_t final_range=0.;
Int_t final_step[10000];
Int_t trk_out[10000];
Int_t step_out=0;
Int_t end_trk=0;
Int_t count_down=0;
Int_t nhit=0;
Int_t ev_Id=0;
Int_t trk_Id=0;
Int_t parent_Id=0;
Int_t step_Id=0;
Int_t pdg_Id=0;
Int_t vol_Id=0;
Int_t copyNo=0;
Int_t emu_Id=0;
Int_t whichInt=0;
Int_t whichCre=0;
Bool_t outside=false;
Bool_t contained_Top=false;
Bool_t contained_Down=false;
Bool_t passed_Top=false;
Bool_t passed_Down=false;
Bool_t last_Top=false;
Bool_t last_Down=false;
Bool_t all_trk=false;
Bool_t trk_father=false;
Bool_t find_father=false;
Bool_t long_primary=false;
Bool_t long_secondary=false;
Int_t multiplicity=0;
Int_t tmp_primary=-1;
Int_t tmp_is_father=0;
Int_t father=0;
Double_t ene_dep_trk=0.0;
Double_t ene_dep_step=0.0;
Double_t trk_len=0.0;
Double_t trk_phi=0.0;
Double_t trk_the=0.0;
Double_t trk_ene=0.0;
Double_t MC_phi=0.0;
Double_t MC_the=0.0;
Double_t MC_ene=0.0;
Double_t x_pre=0.0;
Double_t y_pre=0.0;
Double_t z_pre=0.0;
Double_t x_pos=0.0;
Double_t y_pos=0.0;
Double_t z_pos=0.0;
Double_t delta_ene=0.;
///////////////////////
Int_t ev_id;
Int_t cl_id;
Int_t cl_id_big;
vector<int> cl_id_small;
vector<double> clx_small;
vector<double> cly_small;
vector<double> clz_small;
vector<double> cl_ene_small;
Int_t ntot;
Int_t Ntot;
Double_t cl_dep_ene=0;
Double_t big_dep_ene=0;

//vector<int>     *cl_el;
Int_t iel=0;
//Double_t ev_x_pos[1000]={};
//Double_t ev_y_pos[1000]={};
//Double_t ev_z_pos[1000]={};


//Bool_t visited[1000]={};
//Bool_t visited2[1000]={};
//Int_t cl_el[1000][1000]={};
//Float_t cl_x[1000][1000]={};
//Float_t cl_y[1000][1000]={};
//Float_t cl_z[1000][1000]={};
//Int_t n_el[1000]={};
//Int_t id_el2[1000]={};
Int_t tot_cl=0;

Double_t tmp_x_ev=0;
Double_t tmp_y_ev=0;
Double_t tmp_z_ev=0;

float thr_ene=3.5;
float epsilon1=0.1;
float epsilon2=0.45;
int min_points=1;
int icl=0;
int sphere_points=0;
int sphere_points2=0;
float dist=0;

double minx_cl; 
double miny_cl;
double minz_cl;
double maxx_cl;
double maxy_cl;
double maxz_cl;

char name_gr[100];
char name_grz[100];






//---------------------------
Int_t  GetChargeFromPdg(Int_t pdg){
  Int_t charge = -99;

  if(pdg==1000060120 || pdg==1000060130 || pdg==1000060140 || pdg==1000050110 || pdg==1000040090 || 1000060150) charge=6;
  if(pdg==1000070140 || pdg==1000070150) charge=7;
  if(pdg==1000080160 || pdg==1000080170 || pdg==1000080180) charge=8;
  if(pdg==1000160320 || pdg==1000160330 || pdg==1000160340 ||
     pdg==1000150320 || pdg==1000140290) charge=16;

  if(pdg==1000531270 || pdg==1000531280) charge=53;
  if(pdg==1000350790 || pdg==1000350810 ||  pdg==1000340790 ||
     pdg==1000350800 || pdg==1000350820 ||  pdg==1000340810) charge=35;
  if(pdg==1000471090 || pdg==1000471070 ||  pdg==1000461070 || 
     pdg==1000471080 || pdg==1000471100) charge=47;

  if(pdg==2212 || pdg==1000010020 || pdg==1000010030) charge = 1;

  if(pdg==1000020040) charge = 4;
    
  return charge;
}

//---------------------------
Int_t  GetWeightFromCharge(Int_t charge){
  Int_t weight = -99;

  // weight = 1 -> Proton
  if(charge==1)
    weight = 1;
  
  // weight = 2 -> Light Nuclei
  // C=6, N=7, O=8, S=16
  if(charge==6||charge==7||charge==8||charge==16)
    weight = 2;

  // weight = 3 -> Heavy Nuclei
  // Br=35, Ag=47, I=53
  if(charge==35||charge==47||charge==53)
    weight = 3;
  
  // weight = 4 -> Alfa
  if(charge==4)
    weight = 4;
  
  return weight;
}





void angles_MC(double px_MC,
	    double py_MC,
	    double pz_MC,  
	    double* return_F){
 
  Double_t moduloXY_MC = TMath::Sqrt(px_MC*px_MC+py_MC*py_MC/*+pz_MC*pz_MC*/);  
  //return_F[0] = TMath::ACos(pz_MC/moduloXYZ_MC);
  return_F[0] = TMath::ATan2(moduloXY_MC,pz_MC);
  //if(return_F[0]<=pi/2.) return_F[0] = pi-return_F[0];
  //else return_F[0] = return_F[0] - pi;
  return_F[1] = TMath::ATan2(py_MC,px_MC);  
}	

void trk_Rec(double x_first,
	       double y_first,
	       double z_first,
	       double x_last,
	       double y_last,
	       double z_last, 
	       double* return_F){

  Double_t delta_x = x_last - x_first;
  Double_t delta_y = y_last - y_first;
  Double_t delta_z = z_last - z_first;
  return_F[0] = TMath::Sqrt(TMath::Power(delta_x,2)+TMath::Power(delta_y,2)+TMath::Power(delta_z,2))*1000; // trk_len in nm  
  return_F[1] = TMath::ACos(delta_z/return_F[0]);
  return_F[2] = TMath::ATan2(delta_y,delta_x);
}
