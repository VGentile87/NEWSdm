#define Make_Tree_test_cxx
#include "Make_Tree_test.h"
#include "Definitions.h"
#include <TH2.h>
#include <TH3F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
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
#include <TMultiGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>

using namespace std;

Make_Tree_test v;

void Make_Tree_test::Loop()
{

  float thr_ene;

  cout << "Inserire threshod energy [keV] per la formazione di small clusters" << endl;
  cin >> thr_ene;
  cout << "La threshold utilizzata è " << thr_ene << " keV" << endl;

  
  std::vector<double> ev_x_pos;
  std::vector<double> ev_y_pos;
  std::vector<double> ev_z_pos;
  //std::vector<double> ev_x_pos2;
  //std::vector<double> ev_y_pos2;
  //std::vector<double> ev_z_pos2;
  std::vector<double> ev_ene_dep;
  //std::vector<double> ev_ene_dep2;
  std::vector<bool> visited;
  std::vector<bool> visited2;
  std::vector<int> n_el;
  std::vector<int> id_el2;
  std::vector <std::vector<double> > cl_el;
  std::vector <std::vector<double> > cl_x;
  std::vector <std::vector<double> > cl_y;
  std::vector <std::vector<double> > cl_z;
  std::vector <std::vector<double> > cl_ene_dep;
  std::vector <std::vector<double> > cl_x_abs;
  std::vector <std::vector<double> > cl_y_abs;
  std::vector <std::vector<double> > cl_z_abs;

  int width=400;
  int height=200;
  int depth=1;
  
  int xmin = -20000;
  int xmax = 20000;
  int ymin = -10000;
  int ymax = 10000;
  int zmin = 475;
  int zmax = 525;
  float xrange = (xmax-xmin)/width;
  float yrange = (ymax-ymin)/height;
  float zrange = (zmax-zmin)/depth;

  cout << xmin << " " << xmax << " " << xrange << endl;
  cout << ymin << " " << ymax << " " << yrange << endl;
  cout << zmin << " " << zmax << " " << zrange << endl;
  
  //std::vector<std::vector<std::vector<std::vector<unsigned int> > > > step_list;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > step_list2;
  //std::vector<std::vector<std::vector<std::vector<float> > > > ev_x_pos;
  //std::vector<std::vector<std::vector<std::vector<float> > > > ev_y_pos;
  //std::vector<std::vector<std::vector<std::vector<float> > > > ev_z_pos;
  //std::vector<std::vector<std::vector<std::vector<float> > > > ev_ene_dep;
  std::vector<std::vector<std::vector<std::vector<float> > > > ev_x_pos2;
  std::vector<std::vector<std::vector<std::vector<float> > > > ev_y_pos2;
  std::vector<std::vector<std::vector<std::vector<float> > > > ev_z_pos2;
  std::vector<std::vector<std::vector<std::vector<float> > > > ev_ene_dep2;
  //step_list.resize(width);
  step_list2.resize(width);
  //ev_x_pos.resize(width);
  //ev_y_pos.resize(width);
  //ev_z_pos.resize(width);
  //ev_ene_dep.resize(width);
  ev_x_pos2.resize(width);
  ev_y_pos2.resize(width);
  ev_z_pos2.resize(width);
  ev_ene_dep2.resize(width);
  for(int i=0;i<width;i++){
    //step_list[i].resize(height);
    step_list2[i].resize(height);
    //ev_x_pos[i].resize(height);
    //ev_y_pos[i].resize(height);
    //ev_z_pos[i].resize(height);
    //ev_ene_dep[i].resize(height);
    ev_x_pos2[i].resize(height);
    ev_y_pos2[i].resize(height);
    ev_z_pos2[i].resize(height);
    ev_ene_dep2[i].resize(height);
    for(int j=0;j<height;j++){
      //step_list[i][j].resize(depth);
      step_list2[i][j].resize(depth);
      //ev_x_pos[i][j].resize(depth);
      //ev_y_pos[i][j].resize(depth);
      //ev_z_pos[i][j].resize(depth);
      //ev_ene_dep[i][j].resize(depth);
      ev_x_pos2[i][j].resize(depth);
      ev_y_pos2[i][j].resize(depth);
      ev_z_pos2[i][j].resize(depth);
      ev_ene_dep2[i][j].resize(depth);
    }
  }
  

  Double_t clx_big, cly_big, clz_big;

  // ev_x_pos2.clear();
  //ev_y_pos2.clear();
  //ev_z_pos2.clear();
  //ev_ene_dep2.clear();

  int eps1 = epsilon1*1000;
  int eps2 = epsilon2*1000;

  TFile * f_out = new TFile(Form("data_eps1_%d_eps2_%d_thr_%.1fkeV.root",eps1,eps2,thr_ene),"RECREATE");
  TTree * Tree_cl_big = new TTree("tree_cl_big","data_cl_big");
  //Tree_cl_big->Branch("ev_id",&ev_id,"ev_id/I");
  Tree_cl_big->Branch("cl_id_big",&cl_id_big,"cl_id_big/I");
  Tree_cl_big->Branch("cl_id_small",&cl_id_small);
  Tree_cl_big->Branch("clx_small",&clx_small);
  Tree_cl_big->Branch("cly_small",&cly_small);
  Tree_cl_big->Branch("clz_small",&clz_small);
  Tree_cl_big->Branch("cl_ene_small",&cl_ene_small);
  Tree_cl_big->Branch("big_dep_ene",&big_dep_ene,"big_dep_ene/D");
  Tree_cl_big->Branch("clx_big",&clx_big,"clx_big/D");
  Tree_cl_big->Branch("cly_big",&cly_big,"cly_big/D");
  Tree_cl_big->Branch("clz_big",&clz_big,"clz_big/D");
  Tree_cl_big->Branch("Ntot",&Ntot,"Ntot/I");

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   // cout << nentries << endl;
   
   for (Long64_t jentry=0;jentry<nentries;jentry++){ // sul singolo gamma

     if(!(jentry>=8977600 && jentry<=8979000) && !(jentry>=9136341 && jentry<=9138000) && !(jentry>=26609622 && jentry<=26611000) && !(jentry>=11196720 && jentry<=11198000)  && !(jentry>=27203300 && jentry<=27205000)){  //4.5

     ev_x_pos.clear();
     ev_y_pos.clear();
     ev_z_pos.clear();
     ev_ene_dep.clear();
     visited.clear();
     visited2.clear();
     n_el.clear();
     id_el2.clear();
     cl_el.clear();
     cl_x.clear();
     cl_y.clear();
     cl_z.clear();
     cl_ene_dep.clear();
     cl_x_abs.clear();
     cl_y_abs.clear();
     cl_z_abs.clear();
     
     
     GetEntry(jentry);
     //Long64_t ientry = LoadTree(jentry);
     //log << jentry << endl;
     if(jentry%100000==0)cout << jentry << endl;
     iel=0;

     
     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     //// ALL STEPS
     for(Int_t in=0; in<vTrackId->size(); in++){

       stop_Top=false;
       if(vVolumeId->at(in)==1 && (vTrackId->at(in)!=tmp_trk_Top || Event!=tmp_Event_Top))stop_Top=false;
       if(vVolumeId->at(in)==1 && vEnergyKin->at(in)==vEnergyDep->at(in))stop_Top=true;
       
       
      if(stop_Top){  // gamma_env
	//cout << vPostHitX->at(in) << " " << vPostHitY->at(in) << " " << vPostHitZ->at(in) << endl;
	if(vPostHitX->at(in)>xmin && vPostHitX->at(in)<xmax && vPostHitY->at(in)>ymin && vPostHitY->at(in)<ymax && vPostHitZ->at(in)>zmin && vPostHitZ->at(in)<zmax ){


	  ev_x_pos.push_back(vPostHitX->at(in));
	  ev_y_pos.push_back(vPostHitY->at(in));
	  ev_z_pos.push_back(vPostHitZ->at(in));
	  ev_ene_dep.push_back(vEnergyDep->at(in));
	  

	  iel++;
	}
      }
     }//end loop rinculi



     //// CLUSTERING PER EVENTO
     icl=0;
     sphere_points=0;
     sphere_points2=0;
     
     visited.resize(ev_x_pos.size());
     visited2.resize(ev_x_pos.size());
     n_el.resize(ev_x_pos.size());
     id_el2.resize(ev_x_pos.size());
     cl_el.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_x.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_y.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_z.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_ene_dep.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_x_abs.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_y_abs.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     cl_z_abs.resize(ev_x_pos.size(), std::vector<double>(ev_x_pos.size()));
     
     //cout << "Evento " << jentry << " " << iel << endl;
     
     for(int j=0;j<ev_x_pos.size();j++){
       if(visited.at(j)==false){
	 sphere_points=0;
	 sphere_points2=0;
	 visited[j]=true;
	 cl_el[icl][sphere_points]=j;
	 if(j==0){
	   tmp_x_ev=ev_x_pos.at(j);
	   tmp_y_ev=ev_y_pos.at(j);
	   tmp_z_ev=ev_z_pos.at(j);
	 }
	 cl_x[icl][sphere_points]=ev_x_pos.at(j)-tmp_x_ev;
	 cl_y[icl][sphere_points]=ev_y_pos.at(j)-tmp_y_ev;
	 cl_z[icl][sphere_points]=ev_z_pos.at(j)-tmp_z_ev;
	 cl_ene_dep[icl][sphere_points]=ev_ene_dep.at(j);
	 cl_x_abs[icl][sphere_points]=ev_x_pos.at(j);
	 cl_y_abs[icl][sphere_points]=ev_y_pos.at(j);
	 cl_z_abs[icl][sphere_points]=ev_z_pos.at(j);
	 sphere_points++;
	 for(int k=0;k<ev_x_pos.size();k++){
	   //dist=TMath::Sqrt(TMath::Power(ev_x_pos[k]-ev_x_pos[j],2)+TMath::Power(ev_y_pos[k]-ev_y_pos[j],2)+TMath::Power(ev_z_pos[k]-ev_z_pos[j],2));
	   dist=TMath::Sqrt(TMath::Power(ev_x_pos.at(k)-ev_x_pos.at(j),2)+TMath::Power(ev_y_pos.at(k)-ev_y_pos.at(j),2)+TMath::Power(ev_z_pos.at(k)-ev_z_pos.at(j),2));
	   if(dist<epsilon1 && k!=j && !visited2.at(j)){
	     visited2[k]=true;
	     //cout <<"a "<< ev_ene_dep.at(k) << " " <<  j << " " << k << " " << sphere_points << " " << dist <<  endl;
	     cl_el[icl][sphere_points]=k;	    
	     cl_x[icl][sphere_points]=ev_x_pos.at(k)-tmp_x_ev;
	     cl_y[icl][sphere_points]=ev_y_pos.at(k)-tmp_y_ev;
	     cl_z[icl][sphere_points]=ev_z_pos.at(k)-tmp_z_ev;
	     cl_ene_dep[icl][sphere_points]=ev_ene_dep.at(k);
	     cl_x_abs[icl][sphere_points]=ev_x_pos.at(k);
	     cl_y_abs[icl][sphere_points]=ev_y_pos.at(k);
	     cl_z_abs[icl][sphere_points]=ev_z_pos.at(k);
	     sphere_points++;
	   }
	 }
	 
	 if(sphere_points < (min_points+1) && !visited2.at(j)) {
	   //cout <<"b "<< icl << " " <<  j << " " <<  sphere_points <<  endl;
	   cl_el[icl][sphere_points]=j;
	   cl_x[icl][sphere_points]=ev_x_pos.at(j)-tmp_x_ev;
	   cl_y[icl][sphere_points]=ev_y_pos.at(j)-tmp_y_ev;
	   cl_z[icl][sphere_points]=ev_z_pos.at(j)-tmp_z_ev;
	   cl_ene_dep[icl][sphere_points]=ev_ene_dep.at(j);
	   cl_x_abs[icl][sphere_points]=ev_x_pos.at(j);
	   cl_y_abs[icl][sphere_points]=ev_y_pos.at(j);
	   cl_z_abs[icl][sphere_points]=ev_z_pos.at(j);
	   n_el[icl]=1;
	   //cout <<"cl_el2 "<< icl << " " << n_el[icl] << " " << cl_el[icl][sphere_points] << endl;
	   icl++;
	   //cout << "e " << icl << " " << sphere_points << " " <<  cl_x[icl][sphere_points] << " " << cl_y[icl][sphere_points] << endl;
	   continue;
	 }
	 else {
	   for(int j=0;j<sphere_points;j++){
	     //sphere_points2=0;
	     if(visited.at(cl_el[icl][j])==false){
	       visited[cl_el[icl][j]]=true;
	       for(int k=0;k<iel;k++){
		 // dist = (TMath::Sqrt(TMath::Power(ev_x_pos[k]-ev_x_pos[cl_el[icl][j]],2)+TMath::Power(ev_y_pos[k]-ev_y_pos[cl_el[icl][j]],2)
		 //		     + TMath::Power(ev_z_pos[k]-ev_z_pos[cl_el[icl][j]],2)));
		 dist = (TMath::Sqrt(TMath::Power(ev_x_pos.at(k)-ev_x_pos.at(cl_el[icl][j]),2)+TMath::Power(ev_y_pos.at(k)-ev_y_pos.at(cl_el[icl][j]),2)
				     + TMath::Power(ev_z_pos.at(k)-ev_z_pos.at(cl_el[icl][j]),2)));
		 if(dist<epsilon1 && k!=cl_el[icl][j] && !visited2.at(k) && !visited.at(k)){
		   id_el2[sphere_points2]=k;
		   visited2[k]=true;
		   //cout <<"c "<< icl << " " <<  j << " " << cl_el[icl][j] << " " << k << " " << sphere_points2 <<  " " << dist << " " << id_el2[sphere_points2] << endl;	
		   sphere_points2++;
		 }
	       }
	       if(sphere_points2 >= min_points){
		 for(int i=0;i<sphere_points2;i++){
		   cl_el[icl][sphere_points+i]=id_el2[i];
		   cl_x[icl][sphere_points+i]=ev_x_pos.at(id_el2[i])-tmp_x_ev;
		   cl_y[icl][sphere_points+i]=ev_y_pos.at(id_el2[i])-tmp_y_ev;
		   cl_z[icl][sphere_points+i]=ev_z_pos.at(id_el2[i])-tmp_z_ev;
		   cl_ene_dep[icl][sphere_points+i]=ev_ene_dep.at(id_el2[i]);
		   cl_x_abs[icl][sphere_points+i]=ev_x_pos.at(id_el2[i]);
		   cl_y_abs[icl][sphere_points+i]=ev_y_pos.at(id_el2[i]);
		   cl_z_abs[icl][sphere_points+i]=ev_z_pos.at(id_el2[i]);
		   //cout <<"d "<< icl << " " << i << " " << id_el2[i] << endl;
		   //cout << "e " << icl << " " << (sphere_points+i) << " " <<  cl_x[icl][sphere_points+i] << " " << cl_y[icl][sphere_points+i] << endl;
		 }
	       }
	       //else 
	     }
	   }
	 }
	 //cout << "icl " << icl << " " << sphere_points << " " << sphere_points2 << endl;
	  if(!visited2[j]){
	   // cout << "cl_el ";
	   for(int i=0;i<(sphere_points+sphere_points2);i++){
	     //cout << cl_el[icl][i] << " ";
	     n_el[icl]++;
	   }
	   //cout << endl;
	   //cout <<"n_el "<< icl << " " << n_el[icl] << endl;
	 }
	 //cout <<"f " <<  j << " " << visited[j] << " " << visited2[j] << " " << sphere_points << endl;
	 if(!visited2.at(j) || sphere_points<(min_points))icl++;
	 
       }
       } /// end clustering


     
     //cout << "number of clusters " << icl << endl;
     
     for(int i=0;i<icl;i++){
       clx_big=0;
       cly_big=0;
       clz_big=0;
       cl_dep_ene=0;
       
       //cout <<"gr nel "<<jentry << " " <<  i << " " << n_el[i] << endl;
       for(int j=0;j<n_el[i];j++){
	 //cout << jentry << " " << icl << " " << n_el[i] << endl;
	 cl_dep_ene += cl_ene_dep[i][j];
	 //cout <<"gr pt " <<  i << " " << j << " " <<  cl_x[i][j] << " " << cl_y[i][j] << " " << cl_z[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 //cout <<"gr pt abs " <<  i << " " << j << " " <<  cl_x_abs[i][j] << " " << cl_y_abs[i][j] << " " << cl_z_abs[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 clx_big += cl_x_abs[i][j]*cl_ene_dep[i][j];
	 cly_big += cl_y_abs[i][j]*cl_ene_dep[i][j];
	 clz_big += cl_z_abs[i][j]*cl_ene_dep[i][j];
	 	
       }
       clx_big /= cl_dep_ene;
       cly_big /= cl_dep_ene;
       clz_big /= cl_dep_ene;
       ev_id = jentry;
       //cl_id = icl;
       ntot = n_el[i];
       if(cl_dep_ene>thr_ene){   // gamma env
	 //ev_x_pos2.push_back(clx_big);
	 //ev_y_pos2.push_back(cly_big);
	 //ev_z_pos2.push_back(clz_big);
	 //ev_ene_dep2.push_back(cl_dep_ene);
	 //cout << "eccomi " << jentry << endl;
	 //cl_id = tot_cl;
	 
	 int ix = floor((clx_big-xmin)/xrange);
	 int iy = floor((cly_big-ymin)/yrange);
	 int iz = floor((clz_big-zmin)/zrange);
	 //cout << vx << " " << vy <<" " << vz<< endl;
	 //cout << ix << " " << iy << " " << iz << endl;
	  
	 step_list2.at(ix).at(iy).at(iz).push_back(tot_cl);
	 ev_x_pos2.at(ix).at(iy).at(iz).push_back(clx_big);
	 ev_y_pos2.at(ix).at(iy).at(iz).push_back(cly_big);
	 ev_z_pos2.at(ix).at(iy).at(iz).push_back(clz_big);
	 ev_ene_dep2.at(ix).at(iy).at(iz).push_back(cl_dep_ene);
	 
	 tot_cl++;
       }

       }

   } // avoid crash
     
   } /// end jentry
   cout << "small clusters: " << tot_cl <<", ncl_ev: " << icl << ", epsilon1: " << epsilon1 << ", thr_ene: "<< thr_ene << endl;
   //////////////////////////////////////////////

   //////////// BIG CLUSTERS //////////////////////////////////
   

   for(int i=0;i<width;i++){
     for(int j=0;j<height;j++){
       for(int k=0;k<depth;k++){
	 //cout << "size " << vtx_list.at(i).at(j).at(k).size() << endl;
	 int dim=step_list2.at(i).at(j).at(k).size();
	 if(dim>0){
	   //cout << jentry << " " << i << " " << j << " " << k << " " << dim << endl;
   
	   icl=0;
	   sphere_points=0;
	   sphere_points2=0;
	   visited.clear();
	   visited2.clear();
	   n_el.clear();
	   id_el2.clear();
	   cl_el.clear();
	   cl_x.clear();
	   cl_y.clear();
	   cl_z.clear();
	   cl_ene_dep.clear();
	   cl_x_abs.clear();
	   cl_y_abs.clear();
	   cl_z_abs.clear();
	   visited.clear();
	   visited2.clear();
	   visited.resize(dim);
	   visited2.resize(dim);
	   n_el.resize(dim);
	   id_el2.resize(dim);
	   
	   cl_el.resize(dim);
	   cl_x.resize(dim);
	   cl_y.resize(dim);
	   cl_z.resize(dim);
	   cl_ene_dep.resize(dim);
	   cl_x_abs.resize(dim);
	   cl_y_abs.resize(dim);
	   cl_z_abs.resize(dim);




     //cout << visited.size() << " " << visited2.size() << endl;
     
     
     for(int l=0;l<dim;l++){

       if(l%10000==0)cout <<"SMALL "<< l << " " << dim << endl;
       //cout << j << endl;
       if(visited.at(l)==false){
	 sphere_points=0;
	 sphere_points2=0;
	 visited[l]=true;
	 cl_el[icl].push_back(l);
	 if(j==0){
	   tmp_x_ev=ev_x_pos2.at(i).at(j).at(k).at(l);
	   tmp_y_ev=ev_y_pos2.at(i).at(j).at(k).at(l);
	   tmp_z_ev=ev_z_pos2.at(i).at(j).at(k).at(l);
	 }
	 //cout << "hello " << endl;
	 cl_x[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(l)-tmp_x_ev);
	 cl_y[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(l)-tmp_y_ev);
	 cl_z[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(l)-tmp_z_ev);
	 cl_ene_dep[icl].push_back(ev_ene_dep2.at(i).at(j).at(k).at(l));
	 cl_x_abs[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(l));
	 cl_y_abs[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(l));
	 cl_z_abs[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(l));
	 //cout << "hello1 " << endl;
	 sphere_points++;
	 for(int m=0;m<dim;m++){
	   //dist=TMath::Sqrt(TMath::Power(ev_x_pos2[k]-ev_x_pos2[j],2)+TMath::Power(ev_y_pos2[k]-ev_y_pos2[j],2)+TMath::Power(ev_z_pos2[k]-ev_z_pos2[j],2));
	   dist=TMath::Sqrt(TMath::Power(ev_x_pos2.at(i).at(j).at(k).at(m)-ev_x_pos2.at(i).at(j).at(k).at(l),2)+TMath::Power(ev_y_pos2.at(i).at(j).at(k).at(m)-ev_y_pos2.at(i).at(j).at(k).at(l),2)+TMath::Power(ev_z_pos2.at(i).at(j).at(k).at(m)-ev_z_pos2.at(i).at(j).at(k).at(l),2));
	   if(dist<epsilon2 && m!=l && !visited2.at(l)){
	     visited2[m]=true;
	     //cout <<"a "<< icl << " " <<  j << " " << k << " " << sphere_points << " " << dist <<  endl;
	     //cout << "hello2 " << endl;
	     cl_el[icl].push_back(m);
	     cl_x[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(m)-tmp_x_ev);
	     cl_y[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(m)-tmp_y_ev);
	     cl_z[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(m)-tmp_z_ev);
	     cl_ene_dep[icl].push_back(ev_ene_dep2.at(i).at(j).at(k).at(m));
	     cl_x_abs[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(m));
	     cl_y_abs[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(m));
	     cl_z_abs[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(m));
	     
	     sphere_points++;
	   }
	 }
	 
	 if(sphere_points < (min_points+1) && !visited2.at(l)) {
	   //cout <<"b "<< icl << " " <<  j << " " <<  sphere_points <<  endl;
	   //cout << "hell3 " << endl;
	   cl_el[icl].push_back(l);
	   cl_x[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(l)-tmp_x_ev);
	   cl_y[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(l)-tmp_y_ev);
	   cl_z[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(l)-tmp_z_ev);
	   cl_ene_dep[icl].push_back(ev_ene_dep2.at(i).at(j).at(k).at(l));
	   cl_x_abs[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(l));
	   cl_y_abs[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(l));
	   cl_z_abs[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(l));
	   n_el[icl]=1;
	   //cout <<"cl_el2 "<< icl << " " << n_el[icl] << " " << cl_el[icl][sphere_points] << endl;
	   //cout << "c " << icl << " " << sphere_points << " " <<  cl_x[icl][sphere_points] << " " << cl_y[icl][sphere_points] << endl;
	   //cout << "big " << icl << " " << cl_el[icl][0] << " " << sphere_points << " " << sphere_points2 << endl;
	   icl++;	   
	   continue;
	 }
	 else {
	   //cout << "hello4 " << endl;
	   for(int m=0;m<sphere_points;m++){
	     //sphere_points2=0;
	     if(visited.at(cl_el[icl][m])==false){
	       visited[cl_el[icl][m]]=true;
	       for(int n=0;n<dim;n++){
		 
		 dist = (TMath::Sqrt(TMath::Power(ev_x_pos2.at(i).at(j).at(k).at(n)-ev_x_pos2.at(i).at(j).at(k).at(cl_el[icl][m]),2)+TMath::Power(ev_y_pos2.at(i).at(j).at(k).at(n)-ev_y_pos2.at(i).at(j).at(k).at(cl_el[icl][m]),2) + TMath::Power(ev_z_pos2.at(i).at(j).at(k).at(n)-ev_z_pos2.at(i).at(j).at(k).at(cl_el[icl][m]),2)));
		 if(dist<epsilon2 && n!=cl_el[icl][m] && !visited2.at(n) && !visited.at(n)){
		   id_el2[sphere_points2]=n;
		   visited2[n]=true;
		   //cout <<"c1 "<< icl << " " <<  j << " " << cl_el[icl][j] << " " << k << " " << sphere_points2 <<  " " << dist << " " << id_el2[sphere_points2] << endl;	
		   sphere_points2++;
		 }
	       }
	       if(sphere_points2 >= min_points){
		 for(int p=0;p<sphere_points2;p++){

		   cl_el[icl].push_back(id_el2[p]);
		   cl_x[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(id_el2[p])-tmp_x_ev);
		   cl_y[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(id_el2[p])-tmp_y_ev);
		   cl_z[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(id_el2[p])-tmp_z_ev);
		   cl_ene_dep[icl].push_back(ev_ene_dep2.at(i).at(j).at(k).at(id_el2[p]));
		   cl_x_abs[icl].push_back(ev_x_pos2.at(i).at(j).at(k).at(id_el2[p]));
		   cl_y_abs[icl].push_back(ev_y_pos2.at(i).at(j).at(k).at(id_el2[p]));
		   cl_z_abs[icl].push_back(ev_z_pos2.at(i).at(j).at(k).at(id_el2[p]));
		   //cout <<"d "<< icl << " " << i << " " << id_el2[i] << endl;
		   //cout << "e " << icl << " " << (sphere_points+i) << " " <<  cl_x[icl][sphere_points+i] << " " << cl_y[icl][sphere_points+i] << endl;
		 }
	       }
	       //else 
	     }
	   }
	 }
	 //cout << "icl " << icl << " " << sphere_points << " " << sphere_points2 << endl;
	 if(!visited2[l]){
	   // cout << "cl_el ";
	   for(int q=0;q<(sphere_points+sphere_points2);q++){
	     //cout << "big " << icl << " " << i << " " << cl_el[icl][i] << " " << sphere_points << " " << sphere_points2 << endl;	     
	     n_el[icl]++;
	   }
	   //cout << endl;
	   //cout <<"n_el "<< icl << " " << n_el[icl] << endl;
	 }
	 //cout <<"f " <<  j << " " << visited[j] << " " << visited2[j] << " " << sphere_points << endl;
	 if(!visited2.at(l) || sphere_points<(min_points)){
	   //for(int i=0;i<(sphere_points+sphere_points2);i++){
	   //cout << "big " << icl << " " << i << " " << cl_el[icl][i] << " " << sphere_points << " " << sphere_points2 << endl;	     }
	   icl++;
	 }
       }
       // cout << j << " " << icl << " " << ev_x_pos2[j] << " " << ev_y_pos2[j] << " " << ev_z_pos2[j] << endl;
     } /// end clustering
     //cout << "ntot cl2 " << icl << endl;


     /////// 
     
     for(int r=0;r<icl;r++){
       clx_big=0;
       cly_big=0;
       clz_big=0;
       big_dep_ene=0;
       cl_id_small.clear();
       clx_small.clear();
       cly_small.clear();
       clz_small.clear();
       cl_ene_small.clear();
       //cout <<"gr nel "<< i << " " << n_el[i] << endl;
       for(int s=0;s<n_el[r];s++){
	 big_dep_ene += cl_ene_dep[r][s];
	 //cout <<"gr pt " <<  i << " " << j << " " <<  cl_x[i][j] << " " << cl_y[i][j] << " " << cl_z[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 //cout <<"gr pt abs " <<  i << " " << j << " " <<  cl_x_abs[i][j] << " " << cl_y_abs[i][j] << " " << cl_z_abs[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 clx_big += cl_x_abs[r][s]*cl_ene_dep[r][s];
	 cly_big += cl_y_abs[r][s]*cl_ene_dep[r][s];
	 clz_big += cl_z_abs[r][s]*cl_ene_dep[r][s];
	 cl_id_small.push_back(cl_el[r][s]);
	 clx_small.push_back(cl_x_abs[r][s]);
	 cly_small.push_back(cl_y_abs[r][s]);
	 clz_small.push_back(cl_z_abs[r][s]);
	 cl_ene_small.push_back(cl_ene_dep[r][s]);
       }
       clx_big /= big_dep_ene;
       cly_big /= big_dep_ene;
       clz_big /= big_dep_ene;
       
       cl_id_big = r;
       Ntot = n_el[r];
       Tree_cl_big->Fill();       
     }     
     //////////////////////////////////////////////

     

	 }
       }
     }
   }
     // mettere parentesi

     
     //cout << "big clusters: " << icl <<  ", epsilon2: " << epsilon2 <<  endl;
     //////////////////////////////////////////////
     
    /////// 
   
     for(int i=0;i<icl;i++){
       clx_big=0;
       cly_big=0;
       clz_big=0;
       big_dep_ene=0;
       cl_id_small.clear();
       clx_small.clear();
       cly_small.clear();
       clz_small.clear();
       cl_ene_small.clear();
       //cout <<"gr nel "<< i << " " << n_el[i] << endl;
       for(int j=0;j<n_el[i];j++){
	 big_dep_ene += cl_ene_dep[i][j];
	 //cout <<"gr pt " <<  i << " " << j << " " <<  cl_x[i][j] << " " << cl_y[i][j] << " " << cl_z[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 //cout <<"gr pt abs " <<  i << " " << j << " " <<  cl_x_abs[i][j] << " " << cl_y_abs[i][j] << " " << cl_z_abs[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	clx_big += cl_x_abs[i][j]*cl_ene_dep[i][j];
	cly_big += cl_y_abs[i][j]*cl_ene_dep[i][j];
	clz_big += cl_z_abs[i][j]*cl_ene_dep[i][j];
	 cl_id_small.push_back(cl_el[i][j]);
	 clx_small.push_back(cl_x_abs[i][j]);
	 cly_small.push_back(cl_y_abs[i][j]);
	 clz_small.push_back(cl_z_abs[i][j]);
	 cl_ene_small.push_back(cl_ene_dep[i][j]);
       }
       clx_big /= big_dep_ene;
       cly_big /= big_dep_ene;
       clz_big /= big_dep_ene;
       
       cl_id_big = i;
       Ntot = n_el[i];
       Tree_cl_big->Fill();       
       }     
   //////////////////////////////////////////////
   


   //log.close();
   //Tree_out->Write();
   //Tree_cl_small->Write();
   Tree_cl_big->Write();
      
}

int myrun(){
  v.Loop();  
}



     /*
     DBSCAN(D, epsilon, min_points):
      C = 0
      for each unvisited point P in dataset
            mark P as visited
            sphere_points = regionQuery(P, epsilon)
            if sizeof(sphere_points) < min_points
                  ignore P
            else
                  C = next cluster
                  expandCluster(P, sphere_points, C, epsilon, min_points)

expandCluster(P, sphere_points, C, epsilon, min_points):
      add P to cluster C
      for each point P’ in sphere_points
            if P’ is not visited
                  mark P’ as visited
                  sphere_points’ = regionQuery(P’, epsilon)
                  if sizeof(sphere_points’) >= min_points
                        sphere_points = sphere_points joined with sphere_points’
                  if P’ is not yet member of any cluster
                        add P’ to cluster C

regionQuery(P, epsilon):
      return all points within the n-dimensional sphere centered at P with radius epsilon (including P)

     */



     
     //// CLUSTERING PER EVENTO
     /*
     for(int i=0;i<width;i++){
       for(int j=0;j<height;j++){
	 for(int k=0;k<depth;k++){
	   //cout << "size " << vtx_list.at(i).at(j).at(k).size() << endl;
	   int dim=step_list.at(i).at(j).at(k).size();
	   if(dim>0){
	     //cout << jentry << " " << i << " " << j << " " << k << " " << dim << endl;
	     icl=0;
	     sphere_points=0;
	     sphere_points2=0;
	     
	     visited.resize(dim);
	     visited2.resize(dim);
	     n_el.resize(dim);
	     id_el2.resize(dim);
	     cl_el.resize(dim, std::vector<double>(dim));
	     cl_x.resize(dim, std::vector<double>(dim));
	     cl_y.resize(dim, std::vector<double>(dim));
	     cl_z.resize(dim, std::vector<double>(dim));
	     cl_ene_dep.resize(dim, std::vector<double>(dim));
	     cl_x_abs.resize(dim, std::vector<double>(dim));
	     cl_y_abs.resize(dim, std::vector<double>(dim));
	     cl_z_abs.resize(dim, std::vector<double>(dim));
	     
	     //cout << "Evento " << jentry << " " << iel << endl;
	     
	     for(int l=0;l<dim;l++){	       
	       if(visited.at(l)==false){
		 sphere_points=0;
		 sphere_points2=0;
		 visited[l]=true;
		 cl_el[icl][sphere_points]=l;
		 if(l==0){
		   tmp_x_ev=ev_x_pos.at(i).at(j).at(k).at(l);
		   tmp_y_ev=ev_y_pos.at(i).at(j).at(k).at(l);
		   tmp_z_ev=ev_z_pos.at(i).at(j).at(k).at(l);
		 }
		 cl_x[icl][sphere_points]=ev_x_pos.at(i).at(j).at(k).at(l)-tmp_x_ev;
		 cl_y[icl][sphere_points]=ev_y_pos.at(i).at(j).at(k).at(l)-tmp_y_ev;
		 cl_z[icl][sphere_points]=ev_z_pos.at(i).at(j).at(k).at(l)-tmp_z_ev;
		 cl_ene_dep[icl][sphere_points]=ev_ene_dep.at(i).at(j).at(k).at(l);
		 cl_x_abs[icl][sphere_points]=ev_x_pos.at(i).at(j).at(k).at(l);
		 cl_y_abs[icl][sphere_points]=ev_y_pos.at(i).at(j).at(k).at(l);
		 cl_z_abs[icl][sphere_points]=ev_z_pos.at(i).at(j).at(k).at(l);
		 
		 sphere_points++;
		 
		 
		 for(int m=0;m<dim;m++){
		   
		   dist=TMath::Sqrt(TMath::Power(ev_x_pos.at(i).at(j).at(k).at(m)-ev_x_pos.at(i).at(j).at(k).at(l),2)+TMath::Power(ev_y_pos.at(i).at(j).at(k).at(m)-ev_y_pos.at(i).at(j).at(k).at(l),2)+TMath::Power(ev_z_pos.at(i).at(j).at(k).at(m)-ev_z_pos.at(i).at(j).at(k).at(l),2));
		   if(dist<epsilon1 && m!=l && !visited2.at(l)){
		     visited2[m]=true;
		     

		     if(jentry==65)cout << "mmm " <<jentry << " " << i << " " << j << " " << k << " " <<  dim <<  " " << icl << " " <<   sphere_points << endl;
		     //cl_el[icl][sphere_points]=m;	    
		     cl_x[icl][sphere_points]=ev_x_pos.at(i).at(j).at(k).at(m)-tmp_x_ev;
		     cl_y[icl][sphere_points]=ev_y_pos.at(i).at(j).at(k).at(m)-tmp_y_ev;
		     cl_z[icl][sphere_points]=ev_z_pos.at(i).at(j).at(k).at(m)-tmp_z_ev;
		     cl_ene_dep[icl][sphere_points]=ev_ene_dep.at(i).at(j).at(k).at(m);
		     cl_x_abs[icl][sphere_points]=ev_x_pos.at(i).at(j).at(k).at(m);
		     cl_y_abs[icl][sphere_points]=ev_y_pos.at(i).at(j).at(k).at(m);
		     cl_z_abs[icl][sphere_points]=ev_z_pos.at(i).at(j).at(k).at(m);
		     
		     sphere_points++;
		   }
		 }
		 
		 
		 if(sphere_points < (min_points+1) && !visited2.at(l)){
		   //cout <<"b "<< icl << " " <<  j << " " <<  sphere_points <<  endl;
		   
		   cl_el[icl][sphere_points]=l;
		   cl_x[icl][sphere_points]=ev_x_pos.at(i).at(j).at(k).at(l)-tmp_x_ev;
		   cl_y[icl][sphere_points]=ev_y_pos.at(i).at(j).at(k).at(l)-tmp_y_ev;
		   cl_z[icl][sphere_points]=ev_z_pos.at(i).at(j).at(k).at(l)-tmp_z_ev;
		   cl_ene_dep[icl][sphere_points]=ev_ene_dep.at(i).at(j).at(k).at(l);
		   cl_x_abs[icl][sphere_points]=ev_x_pos.at(i).at(j).at(k).at(l);
		   cl_y_abs[icl][sphere_points]=ev_y_pos.at(i).at(j).at(k).at(l);
		   cl_z_abs[icl][sphere_points]=ev_z_pos.at(i).at(j).at(k).at(l);
		   n_el[icl]=1;
		   
		   //cout <<"cl_el2 "<< icl << " " << n_el[icl] << " " << cl_el[icl][sphere_points] << endl;
		   icl++;
		   //cout << "e " << icl << " " << sphere_points << " " <<  cl_x[icl][sphere_points] << " " << cl_y[icl][sphere_points] << endl;
		   continue;
		 }
		 else {
		   for(int m=0;m<sphere_points;m++){
		   //sphere_points2=0;
		   if(visited.at(cl_el[icl][m])==false){
		   visited[cl_el[icl][m]]=true;
		   for(int n=0;n<iel;n++){
		   
		   dist = (TMath::Sqrt(TMath::Power(ev_x_pos.at(i).at(j).at(k).at(n)-ev_x_pos.at(i).at(j).at(k).at(cl_el[icl][m]),2)+TMath::Power(ev_y_pos.at(i).at(j).at(k).at(n)-ev_y_pos.at(i).at(j).at(k).at(cl_el[icl][m]),2) + TMath::Power(ev_z_pos.at(i).at(j).at(k).at(n)-ev_z_pos.at(i).at(j).at(k).at(cl_el[icl][m]),2)));
		   if(dist<epsilon1 && n!=cl_el[icl][m] && !visited2.at(n) && !visited.at(n)){
		   id_el2[sphere_points2]=n;
		   visited2[n]=true;
		   
		   
		   sphere_points2++;
			 }
			 }
			 if(sphere_points2 >= min_points){
			 for(int p=0;p<sphere_points2;p++){
			 
			 cl_el[icl][sphere_points+p]=id_el2[p];
			 cl_x[icl][sphere_points+p]=ev_x_pos.at(i).at(j).at(k).at(id_el2[p])-tmp_x_ev;
			 cl_y[icl][sphere_points+p]=ev_y_pos.at(i).at(j).at(k).at(id_el2[p])-tmp_y_ev;
			 cl_z[icl][sphere_points+p]=ev_z_pos.at(i).at(j).at(k).at(id_el2[p])-tmp_z_ev;
			 cl_ene_dep[icl][sphere_points+p]=ev_ene_dep.at(i).at(j).at(k).at(id_el2[p]);
			 cl_x_abs[icl][sphere_points+p]=ev_x_pos.at(i).at(j).at(k).at(id_el2[p]);
			 cl_y_abs[icl][sphere_points+p]=ev_y_pos.at(i).at(j).at(k).at(id_el2[p]);
			 cl_z_abs[icl][sphere_points+p]=ev_z_pos.at(i).at(j).at(k).at(id_el2[p]);
			 //cout <<"d "<< icl << " " << i << " " << id_el2[i] << endl;
			 
			 //cout << "e " << icl << " " << (sphere_points+i) << " " <<  cl_x[icl][sphere_points+i] << " " << cl_y[icl][sphere_points+i] << endl;
			 }
			 }
			 //else 
			 }
			 }
			 }
		 //cout << "icl " << icl << " " << sphere_points << " " << sphere_points2 << endl;
		 
		 if(!visited2[l]){
		   // cout << "cl_el ";
		   for(int q=0;q<(sphere_points+sphere_points2);q++){
		     //cout << cl_el[icl][i] << " ";
		     n_el[icl]++;
		   }
		   //cout << endl;
		   //cout <<"n_el "<< icl << " " << n_el[icl] << endl;
		 }
		 //cout <<"f " <<  j << " " << visited[j] << " " << visited2[j] << " " << sphere_points << endl;
		 if(!visited2.at(l) || sphere_points<(min_points))icl++;
		 
	       }
	     } /// end clustering


	     
	     cout << "number of clusters " << icl << endl;
	     
	     for(int r=0;r<icl;r++){
	       clx_big=0;
	       cly_big=0;
	       clz_big=0;
	       cl_dep_ene=0;
	       
	       //cout <<"gr nel "<<jentry << " " <<  i << " " << n_el[i] << endl;
	       for(int s=0;s<n_el[r];s++){
		 //cout << jentry << " " << icl << " " << n_el[i] << endl;
		 cl_dep_ene += cl_ene_dep[r][s];
		 //cout <<"gr pt " <<  i << " " << j << " " <<  cl_x[i][j] << " " << cl_y[i][j] << " " << cl_z[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
		 //cout <<"gr pt abs " <<  i << " " << j << " " <<  cl_x_abs[i][j] << " " << cl_y_abs[i][j] << " " << cl_z_abs[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
		 clx_big += cl_x_abs[r][s]*cl_ene_dep[r][s];
		 cly_big += cl_y_abs[r][s]*cl_ene_dep[r][s];
		 clz_big += cl_z_abs[r][s]*cl_ene_dep[r][s];
		 
	       }
	       clx_big /= cl_dep_ene;
	       cly_big /= cl_dep_ene;
	       clz_big /= cl_dep_ene;
	       ev_id = jentry;
	       //cl_id = icl;
	       ntot = n_el[r];
	       if(cl_dep_ene>thr_ene){   // gamma env
		 //ev_x_pos2.push_back(clx_big);
		 //ev_y_pos2.push_back(cly_big);
		 //ev_z_pos2.push_back(clz_big);
		 //ev_ene_dep2.push_back(cl_dep_ene);
		 //cout << "eccomi " << jentry << endl;
		 //cl_id = tot_cl;
		 
		 int ix = floor((clx_big-xmin)/xrange);
		 int iy = floor((cly_big-ymin)/yrange);
		 int iz = floor((clz_big-zmin)/zrange);
		 //cout << vx << " " << vy <<" " << vz<< endl;
		 //cout << ix << " " << iy << " " << iz << endl;
		 
		 step_list2.at(ix).at(iy).at(iz).push_back(tot_cl);
		 ev_x_pos2.at(ix).at(iy).at(iz).push_back(clx_big);
		 ev_y_pos2.at(ix).at(iy).at(iz).push_back(cly_big);
		 ev_z_pos2.at(ix).at(iy).at(iz).push_back(clz_big);
		 ev_ene_dep2.at(ix).at(iy).at(iz).push_back(cl_dep_ene);
	 
		 tot_cl++;
	       }
	     }
	   }
	 }
       }
     }*/
     // mettere parentesi
