#define Make_SRIM_cxx
#include "Make_SRIM.h"
//#include "Definitions.h"
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include <TF1.h>
#include <TStyle.h>
#include <fstream>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>
#include <TArc.h>

using namespace std;

Make_SRIM v;

void Make_SRIM::Loop()
{

  //TFile * f_out = new TFile("tree_SRIM.root","RECREATE");
  //TTree * Tree_out = new TTree("tree1","data");
  //ofstream log("log_data_105.txt");
  //log.is_open();

  //Tree_out->Branch("ev",&ev,"ev_Id/I");

  
  const int ncr_per_um =15*2;//15;// 20;//22;//0;  //5 per geant4
  const int nbin2=174;
  float lim_nbin1=ncr_per_um*100;//660;//502.5;
  float lim_nbin2=ncr_per_um*100;//660;//504.6;
  float zstep=75; //dist in z tra due layer 2(#crist_overlap)*cryst_diam*crist_sig*5;
  int nlayers=15;//15;  //5 per geant4
  int min_gap=0;
  int z_min_gap=0;
  int ncr_plane=ncr_per_um*ncr_per_um;
  float volume = 4*lim_nbin1*lim_nbin1*zstep*nlayers;
  float cut_eff=0.9;

  Double_t pi = TMath::Pi();
  ofstream log("sphere_coordinates.txt");
  ofstream log_trk("tracks.txt");
  ofstream log_ell("ell.txt");
  log.is_open();
  log_trk.is_open();
  log_ell.is_open();

  TCanvas *cc = new TCanvas("cc","cc",600,600);
  TCanvas *ccb = new TCanvas("ccb","ccb",600,600);
  TCanvas *ccc = new TCanvas("ccc","ccc",600,600);
  //TCanvas *ccd = new TCanvas("ccd","ccd",600,600);

  //TCanvas *c44[20];// = new TCanvas("c4","c4",600,600);
  TCanvas *c4[20];// = new TCanvas("c4","c4",600,600);
  //TGraph *gr[20];// = new TGraph();
  TGraph *gr2[20];// = new TGraph();

  TH1F *h_x_start = new TH1F("x_start","x_start",250,0,25);
  TH2F *h_yz_start = new TH2F("yz_start","yz_start",1000,-500,500,1000,-500,500);
  TH1F *h_phi = new TH1F("phi","phi",100,-pi/2.,pi/2.);
  TH1F *h_phi_bar = new TH1F("phi_bar","phi_bar",100,-pi/2.,pi/2.);
  TH1F *h_phi_maxbar = new TH1F("phi_maxbar","phi_maxbar",100,-pi/2.,pi/2.);
  TH1F *h_the = new TH1F("the","the",100,-pi/2.,pi/2.);
  TH1F *h_len = new TH1F("len","len",100,0.,1000);
  TH1F *h_len_2 = new TH1F("len2","len2",100,0.,1000);
  TH1F *h_len_3 = new TH1F("len3","len3",100,0.,1000);
  TH1F *h_nolen = new TH1F("no_len","no_len",100,0.,1000);
  TH1F *h_sens_len = new TH1F("sens_len","sens_len",100,0.,1000);

  TH1F *h_dist = new TH1F("h_dist","dist",100,0.,1000);

  TH2F *hgrid1 = new TH2F("grid1","grid1",ncr_per_um,-lim_nbin1,lim_nbin1,ncr_per_um,-lim_nbin1,lim_nbin1);
  TH2F *hgrid2 = new TH2F("grid2","grid2",nbin2,-lim_nbin2,lim_nbin2,nbin2,-lim_nbin2,lim_nbin2);
  //TH2F *hgrid2 = new TH2F("grid2","grid2",174,-nbin2,nbin2,174,-nbin2,nbin2);
  TH3F *hgrid3d = new TH3F("grid3d","grid3d",nbin2,-lim_nbin2,lim_nbin2,nbin2,-lim_nbin2,lim_nbin2,110,-100,1000);
  //TH2F *hgrid1 = new TH2F("grid1","grid1",ncr_per_um,-5040,5040,ncr_per_um,-5040,5040);
  //TH2F *hgrid2 = new TH2F("grid2","grid2",455,-5005,5005,455,-5005,5005);
  //TH3F *hgrid3d = new TH3F("grid3d","grid3d",455,-5005,5005,455,-5005,5005,110,-100,1000);
  TH1F *hzcr = new TH1F("zcr","zcr",100,-lim_nbin1*1.5,lim_nbin1*1.5);
  TH1F *hycr = new TH1F("ycr","ycr",100,-lim_nbin1*1.5,lim_nbin1*1.5);
  TH1F *hxcr = new TH1F("xcr","xcr",100,-100,nlayers*zstep*1.5);
  TH1F *hrcr = new TH1F("rcr","rcr",30,0,60);
  TH1F *hdist = new TH1F("dist","dist",250,0,250);
  TH1F *h_phi_nit = new TH1F("phi_nit","phi_nit",16,-pi/2.,pi/2.);
  TH1F *h_phi_nit_cut = new TH1F("phi_nit_cut","phi_nit_cut",16,-pi/2.,pi/2.);
  TH1F *h_phi_nit_gap = new TH1F("phi_nit_gap","phi_nit_gap",16,-pi/2.,pi/2.);
  TH1F *h_phi_nit_ell = new TH1F("phi_nit_ell","phi_nit_ell",16,-pi/2.,pi/2.);
  TH1F *h_phi_nit_gap_len = new TH1F("phi_nit_gap_len","phi_nit_gap_len",16,-pi/2.,pi/2.);
  TH1F *h_phi_nit_plasmon = new TH1F("phi_nit_plasmon","phi_nit_plasmon",16,-pi/2.,pi/2.);
  TH1F *h_len_nit = new TH1F("len_nit","len_nit",600,0,600);
  TH1F *h_len_nit_cut = new TH1F("len_nit_cut","len_nit_cut",600,0,600);

  TH1F *hell = new TH1F("hell","hell",100,0,10);
  TProfile *hle = new TProfile("hle","hle",100,0,10,0,1000);

  TH2F *hle2 = new TH2F("hle2","hle2",100,0,10,100,0,600);
  TH2F *hmajmin = new TH2F("hmm","hmm",100,0,1000,100,0,1000);
  TH2F *hmajphi = new TH2F("hmp","hmp",100,0,1000,16,-pi/2.,pi/2.);
  TH2F *h_len_ncr = new TH2F("len_ncr","len_ncr",50,0,50,100,0,1000);
  
  TProfile *h_length = new TProfile("len_prof","len_prof",100,0,600,0,600);
  TH1F *h_ncr = new TH1F("ncr","ncr",50,0,50);
  
  
  TRandom3 *t1 = new TRandom3();
  TRandom3 *t2 = new TRandom3();
  TRandom3 *t3 = new TRandom3();
  TRandom3 *t4 = new TRandom3();
  TRandom3 *t5 = new TRandom3();
  TRandom3 *t6 = new TRandom3();
  TRandom3 *t7 = new TRandom3();

  float cryst_diam=72.34;
  float cryst_sig=7.04;


  float z_rnd=0.;
  float y_rnd=0.;
  float x_rnd=0.;

  std::vector<float> z_pos;
  std::vector<float> y_pos;
  std::vector<float> x_pos;
  std::vector<int> cr_eff;
  std::vector<float> radius;

  //std::vector<float> zrnd;
  //std::vector<float> yrnd;

  
  float dist=0;
  float tmp_cr_dist=1000;


  float tot_vol=0;

  
  /// FRAME PER IL SINGOLO CRISTALLO + GEL
  for(int i=0;i<ncr_per_um;i++){
      for(int j=0;j<ncr_per_um;j++){
	hgrid1->SetBinContent(i+1,j+1,1);
      }
  }

  // hgrid1->Draw("colz");
  
  cc->cd();
  TGraph *gr = new TGraph();
  // SCORRELAMENTO DEI PIANI 
  for(int k=0;k<nlayers;k++){
    float xfr = t1->Uniform(-35,35);
    float yfr = t2->Uniform(-35,35);

    int icr=0;
    for(int i=0;i<ncr_per_um;i++){
      for(int j=0;j<ncr_per_um;j++){

	//radius.push_back(t5->Gaus(44,6.8)/2.);
	radius.push_back(t5->Gaus(cryst_diam,cryst_sig)/2.);
	//radius.push_back(44/2.);
	z_rnd = hgrid1->GetYaxis()->GetBinCenter(i+1) + xfr;
	y_rnd = hgrid1->GetXaxis()->GetBinCenter(j+1) + yfr;	
	x_rnd = k*zstep;

	if(k==0)gr->SetPoint(icr,z_rnd,y_rnd);
	//cout << z_rnd << " " << y_rnd << " " << x_rnd << endl;

	z_pos.push_back(z_rnd);
	y_pos.push_back(y_rnd);
	x_pos.push_back(x_rnd);
	
	//cout << icr << endl;
	icr++;
    }
  }        
    
    if(k==0){
      gr->Draw("AP");
      for(int s=0;s<z_pos.size();s++){
	TArc *arc = new TArc(z_pos.at(s),y_pos.at(s),radius.at(s));
	arc->SetLineColor(kRed);
	arc->SetLineWidth(4);
	arc->SetFillStyle(0);
	arc->Draw();
      }
      
    }
    
  }

  float delta_y=0;
  int count=0;
  TGraph *grb = new TGraph();
  for(int k=0;k<x_pos.size();k++){
    if(k%ncr_per_um==0){
      count=0;
      delta_y=0;
      y_pos.at(k) = y_pos.at(k) +  t2->Uniform(-35,0);
    }
    if(count>0){
      //cout << k << " " << count << " " << y_pos.at(k) << endl;
      delta_y = abs(y_pos.at(k)-y_pos.at(k-1)) - (radius.at(k)+radius.at(k-1));	
    }
    y_pos.at(k) = y_pos.at(k) -  delta_y + min_gap;
    //cout << k << " " << count << " " << y_pos.at(k) << " " << delta_y <<  " " << radius.at(k) <<  endl;
    if(k%ncr_plane<ncr_per_um)z_pos.at(k) = z_pos.at(k) +  t2->Uniform(-35,0);

    if(k<ncr_plane)grb->SetPoint(k,z_pos.at(k),y_pos.at(k));
    
    count++;
  }

  ccb->cd();
  grb->Draw("AP");
  for(int s=0;s<ncr_plane;s++){
    TArc *arc = new TArc(z_pos.at(s),y_pos.at(s),radius.at(s));
    arc->SetLineColor(kRed);
    arc->SetLineWidth(4);
    arc->SetFillStyle(0);
    arc->Draw();
  }
  

  TGraph *grc = new TGraph();
  float delta_z=0;
  int count2=0;
  float alpha=0;
  float zp=0;
  float yp=0;
  float tmp_gap[3]={};
  float tmp_min_dist=0;
 
  for(int k=0;k<x_pos.size();k++){
    if(k%ncr_plane>(ncr_per_um-1)){
      int ln= k-ncr_per_um;
      if(k%ncr_per_um!=0){
	tmp_gap[0] = sqrt(pow(z_pos.at(ln)-z_pos.at(k),2)+pow(y_pos.at(ln)-y_pos.at(k),2));
	tmp_gap[1] = sqrt(pow(z_pos.at(ln-1)-z_pos.at(k),2)+pow(y_pos.at(ln-1)-y_pos.at(k),2));
	tmp_gap[2] = sqrt(pow(z_pos.at(ln+1)-z_pos.at(k),2)+pow(y_pos.at(ln+1)-y_pos.at(k),2));      
	if(tmp_gap[0]<tmp_gap[1] && tmp_gap[0]<tmp_gap[2])ln=ln;
	if(tmp_gap[1]<tmp_gap[0] && tmp_gap[1]<tmp_gap[2])ln--;
	if(tmp_gap[2]<tmp_gap[0] && tmp_gap[2]<tmp_gap[1])ln++;
      }
      //cout << k << " " << ln << z_pos.at(k) << " " << z_pos.at(ln) << " " << y_pos.at(k) << " " << y_pos.at(ln) << endl;
      alpha = atan2((y_pos.at(k)-y_pos.at(ln)),(z_pos.at(k)-z_pos.at(ln)));
      zp = z_pos.at(ln) + radius.at(ln)*cos(alpha);
      yp = y_pos.at(ln) + radius.at(ln)*sin(alpha);
      delta_z = sqrt(pow(zp-z_pos.at(k),2)+pow(yp-y_pos.at(k),2)) - radius.at(k);
      z_pos.at(k) = z_pos.at(k) - delta_z /*- t4->Uniform(0,delta_z)*/;
      
      //cout << k << " " << ln << " " << z_pos.at(k) << " " << delta_z <<  " " << radius.at(k) <<  " " << alpha << " " << zp << " " << yp << endl;

      int start_id = ncr_plane*(k/ncr_plane);
      
      for(int j=start_id;j<k;j++){
	tmp_min_dist = sqrt(pow(z_pos.at(k)-z_pos.at(j),2)+pow(y_pos.at(k)-y_pos.at(j),2));
	if(tmp_min_dist < (radius.at(j)+radius.at(k)) && j!=k){
	  alpha = asin((y_pos.at(k)-y_pos.at(j))/(radius.at(j)+radius.at(k)));
	  zp = z_pos.at(j) + (radius.at(j)+radius.at(k))*cos(alpha);
	  //cout << k << " " << j << " " << zp << " " << z_pos.at(k) << endl;
	  z_pos.at(k) = zp + radius.at(k);
	  //cout << k << " " << j << " overlap" << endl;
	}
      }
    }
    if(k<ncr_plane)grc->SetPoint(k,z_pos.at(k),y_pos.at(k));
  }


  ccc->cd();
  grc->Draw("AP");
  for(int s=0;s<ncr_plane;s++){
    TArc *arc = new TArc(z_pos.at(s),y_pos.at(s),radius.at(s));
    arc->SetLineColor(kRed);
    arc->SetLineWidth(4);
    arc->SetFillStyle(0);
    arc->Draw();
  }

  //TGraph *grd = new TGraph();

   for(int i=0;i<x_pos.size();i++){
    tmp_min_dist=1000;
    for(int j=(i+1);j<x_pos.size();j++){
      tmp_min_dist = sqrt(pow(z_pos.at(j)-z_pos.at(i),2)+pow(y_pos.at(j)-y_pos.at(i),2)+pow(x_pos.at(j)-x_pos.at(i),2));
      //if(tmp_min_dist < (radius.at(i)+radius.at(j)))cout << i << " " << j << " " << (tmp_min_dist-radius.at(i)-radius.at(j)) << " overlap_zy " << endl;
      while(tmp_min_dist < (radius.at(i)+radius.at(j))){
	if(x_pos.at(i)>x_pos.at(j)){
	  x_pos.at(i)++;
	  x_pos.at(j)--;
	}
	else{
	  x_pos.at(i)--;
	  x_pos.at(j)++;
	}
	tmp_min_dist = sqrt(pow(z_pos.at(j)-z_pos.at(i),2)+pow(y_pos.at(j)-y_pos.at(i),2)+pow(x_pos.at(j)-x_pos.at(i),2));
      }
      if(tmp_min_dist < (radius.at(i)+radius.at(j)))cout << i << " " << j << " " << (tmp_min_dist-radius.at(i)-radius.at(j)) << " overlap_zy " << endl;
    }

    // cout << i << " hey " << x_pos.at(i) << radius.at(i) << endl;
    //cout << i << " " << z_pos.at(i) << " " << y_pos.at(i) << " " << x_pos.at(i) << endl;
    //grd->SetPoint(i,x_pos.at(i),y_pos.at(i));
    //hxcr->Fill(x_pos.at(i));
   }
   
   
  float delta_x=0;
  float tmp_delta_x=1000;
  bool not_overlap=false;
  
  for(int s=0;s<nlayers;s++){
    for(int j=0;j<z_pos.size();j++){
      delta_x=0;
      if(s==0 && j>=(ncr_plane*s) && j<(ncr_plane*(s+1))){
	x_pos.at(j) = x_pos.at(j) +  t2->Uniform(-30,0); //
	  // cout << "primo layer " << j << " " << x_pos.at(j) << endl;
      }
      if(s>0){
	//if(TMath::Abs(z_pos.at(j))<1000 && TMath::Abs(y_pos.at(j))<1000 && abs(x_pos.at(j)-s*zstep)<1){
	//	cout << "secondo layer " << j << " " << x_pos.at(j) << endl;
	if(j>=(ncr_plane*s) && j<(ncr_plane*(s+1))){ // layer attuale
	  int start_id = ncr_plane*(s-1);	  
	  //cout << start_id << " " << ncr_plane << " " << s << endl; 
	  tmp_delta_x=0;

	  float tmp_free_dist=100000;
	  float tmp_cr_x=-1;
	  float cr_dist=0;
	  float tmp_cr_dist=0;
	  not_overlap=false;

	  float tmp_free_dist_2=0;

	  
	  for(int v=0;v<100;v++){
	    for(int i=0;/*start_id*/i<(ncr_plane*s);i++){ // sui layer precedenti
	      cr_dist = sqrt(pow(z_pos.at(j)-z_pos.at(i),2)+pow(y_pos.at(j)-y_pos.at(i),2)+pow((x_pos.at(j)-v)-x_pos.at(i),2));// cr_dist
	      tmp_free_dist_2 = cr_dist - radius.at(i)-radius.at(j); // gap_disponibile

	      // se si intersecano, aggiungo salgo su della parte di intersezione
	      if(v==0 && tmp_free_dist_2<0){		
		delta_x=tmp_free_dist_2-2;
		not_overlap=true;
	      }
	      
	      if(tmp_free_dist_2<0 && !not_overlap){
		not_overlap=true;
		delta_x=v-1;
		//cout << "secondo layer " << j << " " << v << " " << x_pos.at(j) << " " << x_pos.at(i) << " " << tmp_free_dist_2 << endl;
		//cr_dist = sqrt(pow(z_pos.at(j)-z_pos.at(i),2)+pow(y_pos.at(j)-y_pos.at(i),2)+pow((x_pos.at(j)-v+1)-x_pos.at(i),2));
		//tmp_free_dist_2 = cr_dist - radius.at(i)-radius.at(j);
		
		//cout << "terzo layer " << j << " " << v << " " << x_pos.at(j) << " " << x_pos.at(i) << " " << tmp_free_dist_2 << endl;
		//cout << i << " " << j << " " << tmp_free_dist_2 << " overlap " << v << " " << delta_x <<  endl;
		break;
	      }
	    }
	    // if(tmp_free_dist_2<0) cout << "overlap " << tmp_free_dist_2 << " " << v << endl;
	  }
	  
	  //if(not_overlap)delta_x=0;
	  	  
	  //x_pos.at(j) = x_pos.at(j);// - t4->Uniform(0,zstep-40-radius.at(j));// + z_min_gap;
	  x_pos.at(j) = x_pos.at(j) - delta_x;//t4->Uniform(0,delta_x);// + z_min_gap;
	  
	  //cout << j  << " overlap_z " << delta_x << " " << tmp_delta_x << " " << x_pos.at(j) << endl;
	}	
      }
      //cout << s << " " << j << " " << z_pos.at(j) << " " << y_pos.at(j) << " " << x_pos.at(j) << endl;
      //grd->SetPoint(j,x_pos.at(j),y_pos.at(j));
      //hxcr->Fill(x_pos.at(j));
    }
     
  }

  cr_eff.resize(x_pos.size());
  
  for(int i=0;i<x_pos.size();i++){
    tmp_min_dist=1000;
    for(int j=(i+1);j<x_pos.size();j++){
      tmp_min_dist = sqrt(pow(z_pos.at(j)-z_pos.at(i),2)+pow(y_pos.at(j)-y_pos.at(i),2)+pow(x_pos.at(j)-x_pos.at(i),2));
      h_dist->Fill(tmp_min_dist);
      if(tmp_min_dist < (radius.at(i)+radius.at(j)))cout << i << " " << j << " " << (tmp_min_dist-(radius.at(i)+radius.at(j))) << " " << x_pos.at(j) << " " << x_pos.at(i) <<  " overlap_z " << " " << radius.at(i) << " " << radius.at(j) << endl;
    }
    if(i==1){
      for(int j=(i+1);j<x_pos.size();j++){
	if(TMath::Abs(z_pos.at(j))<1000 && y_pos.at(j)<-600){
	  
	  TArc *arc = new TArc(z_pos.at(j),x_pos.at(j),radius.at(j));
	  arc->SetLineColor(kRed);
	  arc->SetLineWidth(4);
	  arc->SetFillStyle(0);
	  arc->Draw();
	}
      }
    }
    tot_vol += (4./3.)*pi*pow(radius.at(i),3);
    	
    hgrid2->Fill(y_rnd,z_rnd);
    hgrid3d->Fill(y_rnd,z_rnd,x_rnd);
    hzcr->Fill(z_pos.at(i));
    hycr->Fill(y_pos.at(i));
    hxcr->Fill(x_pos.at(i));
    hrcr->Fill(radius.at(i));
    log << i << " " << x_pos.at(i) << " " <<  y_pos.at(i) << " " << z_pos.at(i) << " " << radius.at(i) << endl;
    if(t2->Uniform(0,1)>(1-cut_eff))cr_eff[i]=1;
    else cr_eff[i]=0;
    //cout << cr_eff.at(i) << endl;
  }

  // end of the crystal framework
 
  
  ///////// SRIM  //////////////
  
  Float_t x_start=0.;
  Float_t y_start=0.;
  Float_t z_start=0.;
  Float_t x_stop=0.;
  Float_t y_stop=0.;
  Float_t z_stop=0.;

  Float_t bar_y=0;
  Float_t bar_z=0;
  Int_t nstep=0;

  vector<float> x_hit;
  vector<float> y_hit;
  vector<float> z_hit;

 
  
  
  const int dimension=20000;
  std::vector <std::vector<float>> trk_hit;
  trk_hit.clear();
  trk_hit.resize(dimension, std::vector<float>(dimension));

   cout << "hello"<< endl;
  

  Float_t phi=0.;
  Float_t phi_maxbar=0.;
  Float_t the=0.;
  Float_t len=0.;
  Float_t len_2=0.;
  Float_t len_proj=0;

  Int_t tmp_ev1=-1;
  Int_t tmp_ev2=-1;
  Float_t tmp_dist=0;
  Float_t tmp_dist_2=0;
  Int_t tmp_first=-1;
  Int_t tmp_last=-1;
  int index=0;

  int not_sensitive=0;
  int sensitive=0;
  int one_grain_type_1=0;
  int one_grain_type_2=0;
  int two_grains=0;
  int track=0;

  int marginal_inside=0;
  int n_two_marginals_inside=0;
  int n_one_marginal_inside=0;
  int n_zero_marginals_inside=0;

  float phi_nit=0;
  float phi_nit_rnd=0;
  float phi_nit_ell=0;
  float len_nit=0;

  float minor=0;
  float major=0;

  Double_t gap=15; //nm

  int idg4, pdg;

  if (fChain == 0) return;
     
  Long64_t nentries = fChain->GetEntriesFast();
   //log << nentries << endl;
   
  //// FILE FROM SRIM
 for (Long64_t jentry=0; jentry<nentries;jentry++) {
     if(index==10000)break;
     GetEntry(jentry);
     Long64_t ientry = LoadTree(jentry);
     //log << jentry << endl;
     //cout << jentry << endl;

     x_hit.push_back(x);
     y_hit.push_back(y);
     z_hit.push_back(z);
     if(rene==0 && index<10000){   /// alla fine della traccia
       tmp_dist_2=0;
       tmp_first=-1;
       tmp_last=-1;
       x_hit.push_back(x);
       y_hit.push_back(y);
       z_hit.push_back(z);


       /// PUNTI MARGINALI
       for(int i=0;i<x_hit.size();i++){
	 for(int j=i;j<x_hit.size();j++){
	   len_2 = sqrt(pow(z_hit.at(j)-z_hit.at(i),2)+pow(y_hit.at(j)-y_hit.at(i),2));
	   if(len_2>tmp_dist_2){
	     tmp_first=i;
	     tmp_last=j;
	     tmp_dist_2=len_2;
	     phi_maxbar = atan2(z_hit.at(j)-z_hit.at(i),y_hit.at(j)-y_hit.at(i));
	   }
	 }
       }

       //cout << index << " " << z_hit.at(tmp_first) << " " << y_hit.at(tmp_first) << endl;
       //cout << tmp_first << endl;
       trk_hit[index][0]=z_hit.at(tmp_first);
       trk_hit[index][1]=y_hit.at(tmp_first);
       trk_hit[index][2]=x_hit.at(tmp_first);
       trk_hit[index][3]=z_hit.at(tmp_last);
       trk_hit[index][4]=y_hit.at(tmp_last);
       trk_hit[index][5]=x_hit.at(tmp_last);

       
       // cout << trk_hit[index][0] << " " << trk_hit[index][1] << " " << trk_hit[index][2] << " " << trk_hit[index][3] << endl;
       
       // h_yz_start->Fill(z_hit.at(tmp_first),y_hit.at(tmp_first));
       h_yz_start->Fill(z_hit.at(tmp_last),y_hit.at(tmp_last));
       h_len_2->Fill(tmp_dist_2);
       /*if(tmp_dist_2>100)*/h_phi_maxbar->Fill(phi_maxbar);
       x_hit.clear();
       y_hit.clear();
       z_hit.clear();
       index++;
     }
   }
   
  ////////////////////////////////////////////////

  TFile * f_in = new TFile("jp_minor_30keV_V.root","READ");
  TH1F * hjpmin = (TH1F*)f_in->Get("htemp");
   
   for(int i=0;i<index;i++){

     float dist_line=0;
     float z1,y1,z2,y2,x1,x2;
     bool is_inside=false;
     int n_inside=0;
     float dist_cr_hit=0;
     float tmp_dist_cr_hit=-1;
     vector <int> tmp_cr;
     tmp_cr.clear();
     len_nit=0;
     minor=0;
     major=0;
     phi_nit=0;
     phi_nit_rnd=0;
     phi_nit_ell=0;

     int ntimes=0;
     do{
       z1 = t3->Uniform(-2000,-1500); // distribuisco le tracce sul pattern dei cristalli
       y1 = t2->Uniform(-2000,-1500);
       x1 = t5->Uniform(250,750);
       z2 = z1 + trk_hit[i][3] - trk_hit[i][0];
       y2 = y1 + trk_hit[i][4] - trk_hit[i][1];
       x2 = x1 + trk_hit[i][5] - trk_hit[i][2];
       //cout << i << " " << x1 << " " << x2 << endl;
       ntimes++;
     }while((x2<0 || x2>1000) && ntimes<100);
     
     float zAB = (z2-z1);
     float yAB = (y2-y1);
     float xAB = (x2-x1);
     float distAB = sqrt(pow(zAB,2)+pow(yAB,2)+pow(xAB,2));
     len_proj = distAB;

     float len_srim = sqrt(pow(zAB,2)+pow(yAB,2));

     //cout << x1 << " " << y1 << " " << x2 << " " << y2 << endl;
     //cout << trk_hit[i][0] << " " << trk_hit[i][1] << " " << trk_hit[i][2] << " " << trk_hit[i][3] << endl;
     
    for(int j=0;j<z_pos.size();j++){
      if(cr_eff.at(j)>0){ // togli se non vuoi l'efficienza
       is_inside=false;
       dist_cr_hit=0;

       float zAC = z_pos.at(j)-z1;
       float yAC = y_pos.at(j)-y1;
       float xAC = x_pos.at(j)-x1;

       float zBC = z_pos.at(j)-z2;
       float yBC = y_pos.at(j)-y2;
       float xBC = x_pos.at(j)-x2;
                    
       float distAC = sqrt(pow(zAC,2)+pow(yAC,2)+pow(xAC,2));
       float distBC = sqrt(pow(zBC,2)+pow(yBC,2)+pow(xBC,2));

       float zABnorm = (z2-z1)/distAB;
       float yABnorm = (y2-y1)/distAB;
       float xABnorm = (x2-x1)/distAB;
       
       dist_cr_hit = sqrt(pow(yAC*xABnorm - xAC*yABnorm,2) + pow(-zAC*xABnorm + xAC*zABnorm,2) + pow(zAC*yABnorm - yAC*zABnorm,2));       
       float proj1 = sqrt(pow(distAC,2)-pow(dist_cr_hit,2));
       float proj2 = sqrt(pow(distBC,2)-pow(dist_cr_hit,2));
       
   
       //if((dist_cr_hit < radius.at(j)) && (proj1+proj2)<=distAB){
       if((distAC<radius.at(j) || distBC<radius.at(j)) || (dist_cr_hit < radius.at(j) && abs(proj1+proj2 - distAB)<1)){
	 //tmp_dist_cr_hit=dist_cr_hit;
	 is_inside=true;
	 //tmp_cr=j;
	 tmp_cr.push_back(j);
	 //cout << n_inside << " " << j << endl;
	 n_inside++;
       }
      }
     }

     //ecco

     
     if(n_inside==0) {
       h_nolen->Fill(len_proj);
       not_sensitive++;
     }
       if(n_inside>0) sensitive++;
       if(n_inside==1){
	 phi_nit = t6->Uniform(-1.57,1.57);
	 phi_nit_rnd = t6->Uniform(-1.57,1.57);
	 h_phi_nit->Fill(phi_nit);
	 //h_phi_nit_cut->Fill(phi_nit);
	 h_phi_nit_gap->Fill(phi_nit);
	 //h_len_nit->Fill(radius.at(tmp_cr.at(0))*2);
	 h_sens_len->Fill(len_proj);
	 //h_length->Fill(len_srim,len_nit);
	 one_grain_type_2++;

	 /*
	 float JP_pix=55; //nm
	 float ell_grain =  hjpmin->GetRandom();
	 float first_gr = ell_grain;
	 float second_gr = ell_grain + t3->Uniform(0,2*55*sqrt(2));
	 if(first_gr>second_gr){
	   major=first_gr;
	   minor=second_gr;
	 }
	 else{
	   minor=first_gr;
	   major=second_gr;
	 }
	 
	 if(major/minor>1.5)h_phi_nit_ell->Fill(phi_nit_rnd);
	 hell->Fill(major/minor);
	 hle->Fill(major/minor,len_nit);
	 hle2->Fill(major/minor,len_nit);
	 hmajmin->Fill(minor,major);
	 */
       }

       if(n_inside>1){

	 int id_tmp_cr1=-1;
	 int id_tmp_cr2=-1;
	 
	 
	 h_sens_len->Fill(len_proj);
	 float tmp_dist_cr=0;
	 for(int k=0;k<tmp_cr.size();k++){
	   for(int j=(k+1);j<tmp_cr.size();j++){
	     //if(t2->Uniform(0,1)>0.5)break;
	     float dist_cr = sqrt(pow(z_pos.at(tmp_cr.at(k))-z_pos.at(tmp_cr.at(j)),2)+pow(y_pos.at(tmp_cr.at(k))-y_pos.at(tmp_cr.at(j)),2)+pow(x_pos.at(tmp_cr.at(k))-x_pos.at(tmp_cr.at(j)),2));

	     //cout << k << " " << j << " " << dist_cr << " " << tmp_dist_cr << endl;
	     
	     if(dist_cr>tmp_dist_cr){
	       tmp_dist_cr=dist_cr;
	       float theta1 = TMath::ACos(t5->Uniform(-1,1));
	       float phi1 = t4->Uniform(0,2*pi);
	       float theta2 = TMath::ACos(t3->Uniform(-1,1));
	       float phi2 = t2->Uniform(0,2*pi);
	       float r1 = /*t1->Uniform(0,radius.at(tmp_cr.at(k)))*/radius.at(tmp_cr.at(k))+70;// + t1->Gaus(cryst_diam,cryst_sig);
	       float r2 = /*t1->Uniform(0,radius.at(tmp_cr.at(j)))*/radius.at(tmp_cr.at(j))+70;// + t1->Gaus(cryst_diam,cryst_sig);
	       float z1_cr_pos = z_pos.at(tmp_cr.at(k)) + r1*TMath::Sin(theta1)*TMath::Cos(phi1);
	       float y1_cr_pos = y_pos.at(tmp_cr.at(k)) + r1*TMath::Sin(theta1)*TMath::Sin(phi1);
	       float z2_cr_pos = z_pos.at(tmp_cr.at(j)) + r2*TMath::Sin(theta2)*TMath::Cos(phi2);
	       float y2_cr_pos = y_pos.at(tmp_cr.at(j)) + r2*TMath::Sin(theta2)*TMath::Sin(phi2);
	       len_nit = sqrt(pow(z2_cr_pos-z1_cr_pos,2)+pow(y2_cr_pos-y1_cr_pos,2));

	       
	       phi_nit_rnd = atan((z2_cr_pos-z1_cr_pos) / (y2_cr_pos-y1_cr_pos));
	       phi_nit_ell = t7->Gaus(phi_nit_rnd,0.420);
	       phi_nit_rnd = t6->Gaus(phi_nit_rnd,0.230);
	       //cout << len_nit << endl;
	       //cout << i << " " << tmp_dist_cr << " " << len_nit << " " << z1_cr_pos << " " << y1_cr_pos << " " << z2_cr_pos << " " << y2_cr_pos << endl;

	       float zAC1 = z_pos.at(tmp_cr.at(k))-z1;
	       float yAC1 = y_pos.at(tmp_cr.at(k))-y1;
	       float xAC1 = x_pos.at(tmp_cr.at(k))-x1;	       
	       float zBC1 = z_pos.at(tmp_cr.at(j))-z1;
	       float yBC1 = y_pos.at(tmp_cr.at(j))-y1;
	       float xBC1 = x_pos.at(tmp_cr.at(j))-x1;
	       float zAC2 = z_pos.at(tmp_cr.at(k))-z2;
	       float yAC2 = y_pos.at(tmp_cr.at(k))-y2;
	       float xAC2 = x_pos.at(tmp_cr.at(k))-x2;	       
	       float zBC2 = z_pos.at(tmp_cr.at(j))-z2;
	       float yBC2 = y_pos.at(tmp_cr.at(j))-y2;
	       float xBC2 = x_pos.at(tmp_cr.at(j))-x2;

	       float dist_pt1_cr1 = sqrt(pow(zAC1,2)+pow(yAC1,2)+pow(xAC1,2));
	       float dist_pt2_cr1 = sqrt(pow(zAC2,2)+pow(yAC2,2)+pow(xAC2,2));
	    
	       float dist_pt1_cr2 = sqrt(pow(zBC1,2)+pow(yBC1,2)+pow(xBC1,2));	       
	       float dist_pt2_cr2 = sqrt(pow(zBC2,2)+pow(yBC2,2)+pow(xBC2,2));
	       marginal_inside=0;

	       if(dist_pt1_cr1 < dist_pt2_cr1 && dist_pt1_cr1 < r1)marginal_inside++;
	       if(dist_pt1_cr1 > dist_pt2_cr1 && dist_pt2_cr1 < r1)marginal_inside++;

	       if(dist_pt1_cr2 < dist_pt2_cr2 && dist_pt1_cr2 < r2)marginal_inside++;
	       if(dist_pt1_cr2 > dist_pt2_cr2 && dist_pt2_cr2 < r2)marginal_inside++;


	       id_tmp_cr1=k;
	       id_tmp_cr2=j;
	       
	     }
	     
	   }
	 }

	 dist_line=0;
	 
	 for(int k=0;k<tmp_cr.size();k++){	   


	   /// redifinizione di minor
	   if(k!=id_tmp_cr1 && k!=id_tmp_cr2){

	     float a_coeff = (y_pos.at(tmp_cr.at(id_tmp_cr2)) - y_pos.at(tmp_cr.at(id_tmp_cr1)))/(z_pos.at(tmp_cr.at(id_tmp_cr2)) - z_pos.at(tmp_cr.at(id_tmp_cr1)));
	     float b_coeff = -1;
	     float c_coeff = y_pos.at(tmp_cr.at(id_tmp_cr1)) - z_pos.at(tmp_cr.at(id_tmp_cr1))*a_coeff;

	     float dist_from_line = (a_coeff*z_pos.at(tmp_cr.at(k)) + b_coeff*y_pos.at(tmp_cr.at(k)) + c_coeff)/sqrt(a_coeff*a_coeff + b_coeff*b_coeff);
	     if(dist_from_line>dist_line)dist_line=dist_from_line;
	   }

	   float JP_pix=55; //nm
	   float first_gr = hjpmin->GetRandom();
	   float second_gr = hjpmin->GetRandom();
	   float third_gr = hjpmin->GetRandom();
	   //float first_gr = t1->Gaus(208,12.97);
	   //float second_gr = t2->Gaus(208,12.97);
	   //float gr_overlap = len_nit_ell - first_gr - second_gr;
	   //minor = first_gr + second_gr;
	   float tmp_minor = (first_gr + second_gr)/2.;
	   if(dist_line>0)tmp_minor = (first_gr + second_gr)/4. + third_gr/2. + dist_line;
	   float tmp_major = (first_gr/2. + second_gr/2. + len_nit);
	   if(tmp_minor<tmp_major){
	     minor = tmp_minor;
	     major = tmp_major;
	   }
	   else{
	     minor = tmp_major;
	     major = tmp_minor;
	   }
	 }

	 phi_nit = atan((z2-z1)/(y2-y1));
	 phi_nit = t7->Gaus(phi_nit,0.230);
	 
	 h_phi_nit->Fill(phi_nit);
	 h_len_nit->Fill(len_nit);
	 h_len_nit_cut->Fill(len_nit);
	 h_length->Fill(len_srim,len_nit);
	 h_len_3->Fill(len_srim);
	 h_phi_nit_gap->Fill(phi_nit_rnd);

	  if(major/minor>1.5){
	   h_phi_nit_ell->Fill(phi_nit_ell);
	   hmajmin->Fill(minor,major);
	 }
	 

	 if(len_nit<185){
	   phi_nit = t6->Uniform(-1.57,1.57);
	   h_phi_nit_gap_len->Fill(phi_nit);
	 }
	 else h_phi_nit_gap_len->Fill(phi_nit_rnd);

	 if(len_nit>120 && len_nit<185){
	   phi_nit = t1->Uniform(-1.57,1.57);
	   h_phi_nit_plasmon->Fill(phi_nit);	   
	 }
	 if(len_nit>185)h_phi_nit_plasmon->Fill(phi_nit_rnd);

	 hell->Fill(major/minor);
	 hle->Fill(major/minor,len_nit);
	 hle2->Fill(major/minor,len_nit);	 
	 hmajphi->Fill(major,phi_nit);

	 //if((len_nit+gap)<100){
	 //phi_nit = t6->Uniform(-1.57,1.57);
	 //}
	 //h_phi_nit_gap->Fill(phi_nit);
	 
	 //if(len_nit<100){
	 // phi_nit = t6->Uniform(-1.57,1.57);
	 //}
	 if(n_inside>1) h_phi_nit_cut->Fill(phi_nit);
	 if(len_nit>100)track++;
	 two_grains++;	 
       }


       if(marginal_inside==2)n_two_marginals_inside++;
       if(marginal_inside==1)n_one_marginal_inside++; 
       if(marginal_inside==0)n_zero_marginals_inside++;

       h_ncr->Fill(n_inside);
       if(n_inside>1)h_len_ncr->Fill(n_inside, len_nit);
       //cout << i << " " << n_inside << " " << len_srim << " " << len_nit << endl;
       cout << i << " " << minor << " " << major << " " << len_nit << " " << dist_line << endl;
       log_trk << i << " " << z1 << " " << y1 << " " << x1 << " " << z2 << " " << y2 << " " << x2 << " " << minor << " " << major << " " << len_nit << " " << phi_nit_ell << endl;
       log_ell << i << " " << minor << " " << major << " " << len_nit << " " << phi_nit_ell << endl;
   }

   log.close();
   log_trk.close();
   log_ell.close();
   cout << "Volume totale dei cristalli " << tot_vol << endl;
   cout << "Volume a disposizione " << volume << endl;
   cout << "No crystal sensitized "<< not_sensitive << endl;
   cout << "Total number of visible events "<< sensitive << endl;
   //cout << "Inside one crystal "<< one_grain_type_1 << endl;
   cout << "One crystal sensitized, missed one point "<< one_grain_type_2 << endl;
   cout << "Two crystals sensitized "<< two_grains << endl;
   cout << "Track length larger than 100nm "<< track << endl;
   cout << "Marginal points inside crystals "<< n_two_marginals_inside << endl;
   cout << "One marginal point inside crystals "<< n_one_marginal_inside << endl;
   cout << "Zero marginal points inside crystals "<< n_zero_marginals_inside << endl;
   //log.close();
   //Tree_out->Write();


   TCanvas *c0 = new TCanvas("c0","c0",600,600);
   h_dist->Draw("");
   
   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   c1->Divide(2,2);
   c1->cd(1);
   hzcr->Draw("");
   c1->cd(2);
   hycr->Draw("");
   c1->cd(3);
   hxcr->Draw("");
   c1->cd(4);
   hrcr->Draw("");

   TCanvas *cc2 = new TCanvas("cc2","cc2",600,600);
   cc2->Divide(2,1);
   cc2->cd(1);
   h_nolen->Draw("");
   cc2->cd(2);
   h_sens_len->Draw("");
  
   TCanvas *c12 = new TCanvas("c12","c12",600,600);
   c12->Divide(2,1);
   c12->cd(1);
   h_len_3->Draw("HIST");
   c12->cd(2);
   h_len_nit->Draw("HIST");

    TCanvas *c15 = new TCanvas("c15","c15",600,600);
   //h_phi_nit_plasmon->Scale(1./h_phi_nit_plasmon->Integral());
   hell->Draw("HIST");

   TCanvas *c16 = new TCanvas("c16","c16",600,600);
   //h_phi_nit_plasmon->Scale(1./h_phi_nit_plasmon->Integral());
   c16->Divide(2,1);
   c16->cd(1);
   hle->Draw("");
   c16->cd(2);
   hle2->Draw("");

   TCanvas *c18 = new TCanvas("c18","c18",600,600);
   hmajmin->Draw("colz");

   TCanvas *c19 = new TCanvas("c19","c19",600,600);
   hmajphi->Draw("colz");

   TCanvas *c20 = new TCanvas("c20","c20",600,600);
   h_ncr->Draw("");
   TCanvas *c21 = new TCanvas("c21","c21",600,600);
   h_len_ncr->Draw("");

}

int myrun(){
  v.Loop();  
}
