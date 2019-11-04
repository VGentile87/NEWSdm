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
  std::vector<double> ev_x_pos;
  std::vector<double> ev_y_pos;
  std::vector<double> ev_z_pos;
  std::vector<double> ev_x_pos2;
  std::vector<double> ev_y_pos2;
  std::vector<double> ev_z_pos2;
  std::vector<double> ev_ene_dep;
  std::vector<double> ev_ene_dep2;
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

  Double_t clx_mean, cly_mean, clz_mean;

  ev_x_pos2.clear();
  ev_y_pos2.clear();
  ev_z_pos2.clear();
  ev_ene_dep2.clear();

  TFile * f_out = new TFile("data.root","RECREATE");
  //TTree * Tree_out = new TTree("tree1","data");
  ofstream log("log_data.txt");
  log.is_open();

  /*
  Tree_out->Branch("ev_Id",&ev_Id,"ev_Id/I");
  Tree_out->Branch("trk_Id",&trk_Id,"trk_Id/I");
  Tree_out->Branch("parent_Id",&parent_Id,"parent_Id/I");
  Tree_out->Branch("step_Id",&step_Id,"step_Id/I");
  Tree_out->Branch("emu_Id",&emu_Id,"emu_Id/I");
  Tree_out->Branch("pdg_Id",&pdg_Id,"pdg_Id/I");
  Tree_out->Branch("whichInt",&whichInt,"whichInt/I");
  Tree_out->Branch("whichCre",&whichCre,"whichCre/I");
  Tree_out->Branch("particle_in",&particle_in,"particle_in/B");
  Tree_out->Branch("particle_out",&particle_out,"particle_out/B");
  Tree_out->Branch("contained_Top",&contained_Top,"contained_Top/B");
  Tree_out->Branch("contained_Down",&contained_Down,"contained_Down/B");
  Tree_out->Branch("passed_Top",&passed_Top,"passed_Top/B");
  Tree_out->Branch("passed_Down",&passed_Down,"passed_Down/B");
  Tree_out->Branch("start_Top",&start_Top,"start_Top/B");
  Tree_out->Branch("start_Down",&start_Down,"start_Down/B");
  Tree_out->Branch("stop_Top",&stop_Top,"stop_Top/B");
  Tree_out->Branch("stop_Down",&stop_Down,"stop_Down/B");
  Tree_out->Branch("outside",&outside,"outside/B");
  Tree_out->Branch("trk_father",&trk_father,"trk_father/B");
  Tree_out->Branch("all_trk",&all_trk,"all_trk/B");
  Tree_out->Branch("ene_dep_trk",&ene_dep_trk,"ene_dep_trk/D");
  Tree_out->Branch("ene_dep_step",&ene_dep_step,"ene_dep_step/D");
  Tree_out->Branch("trk_len",&trk_len,"trk_len/D");
  Tree_out->Branch("trk_phi",&trk_phi,"trk_phi/D");
  Tree_out->Branch("trk_the",&trk_the,"trk_the/D");
  Tree_out->Branch("trk_ene",&trk_ene,"trk_ene/D");
  Tree_out->Branch("copyNo",&copyNo,"copyNo/I");
  Tree_out->Branch("MC_phi",&MC_phi,"MC_phi/D");
  Tree_out->Branch("MC_the",&MC_the,"MC_the/D");
  Tree_out->Branch("MC_ene",&MC_ene,"MC_ene/D");
  Tree_out->Branch("step_len",&step_len,"step_len/D");
  Tree_out->Branch("trk_poly_len",&trk_poly_len,"trk_poly_len/D");
  Tree_out->Branch("final_range",&final_range,"final_range/D");
  Tree_out->Branch("dE_dx",&dE_dx,"dE_dx/D");
  Tree_out->Branch("first_ene",&first_ene,"first_ene/D");
  Tree_out->Branch("delta_ene",&delta_ene,"delta_ene/D");
  Tree_out->Branch("sum_step_len",&sum_step_len,"sum_step_len/D");
  Tree_out->Branch("x_pre",&x_pre,"x_pre/D");
  Tree_out->Branch("y_pre",&y_pre,"y_pre/D");
  Tree_out->Branch("z_pre",&z_pre,"z_pre/D");
  Tree_out->Branch("x_pos",&x_pos,"x_pos/D");
  Tree_out->Branch("y_pos",&y_pos,"y_pos/D");
  Tree_out->Branch("z_pos",&z_pos,"z_pos/D");
  */
  TTree * Tree_cl_big = new TTree("tree_cl_big","data_cl_big");
  //Tree_cl_big->Branch("ev_id",&ev_id,"ev_id/I");
  Tree_cl_big->Branch("cl_id_big",&cl_id_big,"cl_id_big/I");
  Tree_cl_big->Branch("cl_id_small",&cl_id_small);
  Tree_cl_big->Branch("clx_small",&clx_small);
  Tree_cl_big->Branch("cly_small",&cly_small);
  Tree_cl_big->Branch("clz_small",&clz_small);
  Tree_cl_big->Branch("cl_ene_small",&cl_ene_small);
  Tree_cl_big->Branch("big_dep_ene",&big_dep_ene,"big_dep_ene/D");
  Tree_cl_big->Branch("Ntot",&Ntot,"Ntot/I");
  /*
  TTree * Tree_cl_small = new TTree("tree_cl_small","data_cl_small");
  Tree_cl_small->Branch("ev_id",&ev_id,"ev_id/I");
  Tree_cl_small->Branch("cl_id",&cl_id,"cl_id/I");
  Tree_cl_small->Branch("ntot",&ntot,"ntot/I");
  Tree_cl_small->Branch("cl_dep_ene",&cl_dep_ene,"cl_dep_ene/D");
  Tree_cl_small->Branch("clx_mean",&clx_mean,"clx_mean/D");
  Tree_cl_small->Branch("cly_mean",&cly_mean,"cly_mean/D");
  Tree_cl_small->Branch("clz_mean",&clz_mean,"clz_mean/D");
  //Tree_cl_small->AddFriend(Tree_cl_big);
  */


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   cout << nentries << endl;
   
   for (Long64_t jentry=0; jentry<100000/*nentries*/;jentry++) {

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
     Long64_t ientry = LoadTree(jentry);
     //log << jentry << endl;
     if(jentry%1000==0)cout << jentry << endl;
     iel=0;
     end_trk=0;
     out_trk=false;
     step_out=0;
     father=0;
     tmp_is_father=0;
     long_primary=false;
     long_secondary=false;
     //cout << jentry <<" " <<  step_out << endl;
     
     for(int i=0;i<10000;i++){
       final_step[i]=0;
       trk_out[i]=-1;
       is_a_father[i]=-1;
     }


     
     /////// MC    NOT USED YET
     for(Int_t in=0; in<vMC_EnergyKin->size(); in++){
       // i generati //// MC
       if(vMC_ParentId->at(in)==0 && vMC_StepId->at(in)==1){
	 h_Ene_MC->Fill(vMC_EnergyKin->at(in)/1000.); // MeV
	 ////log << vFirstHitPX->at(in)<<endl;
	 angles_MC(vMC_PX->at(in),vMC_PY->at(in),vMC_PZ->at(in),MC_angles);
	 MC_the = MC_angles[0];
	 MC_phi = MC_angles[1];
	 h_CosTheta_MC->Fill(TMath::Cos(MC_angles[0]));
	 h_Theta_MC->Fill(MC_angles[0]);
	 h_Phi_MC->Fill(MC_angles[1]);
	 ////log << vEnergy->at(in) << endl;
       }
     }
     ////////////

     for(Int_t in=0; in<vTrackId->size(); in++){
       if(in<(vTrackId->size()-1)){
	 //cout << in << " " << vTrackId->size() << endl;
	 if(vVolumeId->at(in)==1 && vTrackId->at(in)!=vTrackId->at(in+1)){
	   final_step[end_trk]=in;
	   end_trk++;
	 }
       }
       /*if(in>=(vTrackId->size()-1) && vTrackId->size()!=1){
 	 if(vVolumeId->at(in)==1 && vTrackId->at(in)!=vTrackId->at(in-1)){
	   final_step[end_trk]=in;
	   end_trk++;
	 }
       }*/
       if(in==(vTrackId->size()-1)){
	 if(vVolumeId->at(in)==1){
	   final_step[end_trk]=in;
	   end_trk++;
	 }
       }
       if(vTrackId->size()==1){
	 if(vVolumeId->at(in)==1){
	   final_step[end_trk]=in;
	   end_trk++;
	 }
       }
       //cout << Event << " " << vTrackId->at(in) << " " << vStepId->at(in) << " " << final_step[end_trk] << " " << end_trk-1 << endl;
       if(vInteractionId->at(in)==-1){
	 if(trk_out[step_out]!=vTrackId->at(in) && out_trk)step_out++;
	 out_trk=true;
	 trk_out[step_out]=vTrackId->at(in);
	 //if(Event==98350 || Event==1640)cout << Event << " " << vTrackId->at(in) << " " << vStepId->at(in) << " " << trk_out[step_out] << " " << step_out << endl;
       }
       if(tmp_is_father!=vParentId->at(in) && vParentId->at(in)!=1){
	 is_a_father[father]=vParentId->at(in);
	 tmp_is_father=vParentId->at(in);
	 father++;
       }
            
     }


     ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
     
     tmp_is_father=0;
     father=0;
     out_trk=false;
     end_trk=0;
     step_out=0;
     //// ALL STEPS
     for(Int_t in=0; in<vTrackId->size(); in++){
       nhit++;
       if(jentry<10)log << "index "<< nhit << endl;
    
       //ev_Id=0;
       //trk_Id=0;
       //step_Id=0;
       //emu_Id=0;
       //pdg_Id=0;
       //first_x_pos = 0.0;
       //first_y_pos = 0.0;
       //first_z_pos = 0.0;
       last_x_pos = 0.0;
       last_y_pos = 0.0;
       last_z_pos = 0.0;
       passed_Top=false;
       passed_Down=false;
       contained_Top=false;
       contained_Down=false;
       //outside=false;
       last_Top=false;
       last_Down=false;      
       stop_Top=false;
       stop_Down=false;
       all_trk=false;
       trk_father=false;
       particle_in=false;
       particle_out=false;
       //find_father=false;
       ene_dep_trk=0;
       ev_Id=Event;
       trk_Id=vTrackId->at(in);
       parent_Id=vParentId->at(in);
       step_Id=vStepId->at(in);
       pdg_Id=vPartId->at(in);
       copyNo=vCopyNo->at(in);
       if(vVolumeId->at(in)==1)emu_Id=1;
       if(vVolumeId->at(in)==2)emu_Id=2;
       whichInt = vInteractionId->at(in);
       whichCre = vCreatorId->at(in);
       ene_dep_step = vEnergyDep->at(in);
       step_len = vStepLength->at(in);
       trk_poly_len = vTrackLength->at(in);

       if(in <= final_step[end_trk]){
	 final_range = vTrackLength->at(final_step[end_trk]);
	 //cout << end_trk << " " << final_step[end_trk] << endl;
	 //cout << Event << " " << vTrackId->at(in) << " " << vStepId->at(in) << " " << trk_poly_len << " " << final_range << endl;
	 if(in == final_step[end_trk])end_trk++;
       }
       
       
       if(trk_Id!=trk_out[step_out]){
	 if(outside==true)step_out++;
	 outside=false;
       }
       if(trk_Id==trk_out[step_out])outside=true;

       if(trk_Id!=is_a_father[father]){
	 if(!find_father && in!=0)father++;
	 find_father=true;
	 trk_father=false;
       }
       if(trk_Id==is_a_father[father]){
	 find_father=false;
	 trk_father=true;
       }
       // if(jentry<10) cout <<"c "<< Event << " " << vParentId->at(in) << " " << vTrackId->at(in) << " " << vStepId->at(in) << " " << father << " " << trk_father << endl; 
       
       //if(Event==98350 || Event==5 || Event==1640)cout << Event << " " << trk_Id << " " << trk_out[step_out] << " " << outside << " " << step_out << endl; 
       
       if(!pre_last){
	 pre_last=false;
	 hit_x_pre_last = 0.0;
	 hit_y_pre_last = 0.0;
	 hit_z_pre_last = 0.0;
       }
       
       if(jentry<10)log << "Emulsione " << emu_Id << endl;
       if(jentry<10)log << in << endl;

       if(vVolumeId->at(in)==5 && (vTrackId->at(in)!=tmp_in_trk || Event!=tmp_in_Event) && vParentId->at(in)==0){
	 particle_in=true;
	 tmp_in_trk = vTrackId->at(in);
	 tmp_in_Event  = Event;
       }

       if(vVolumeId->at(in)==4 && (vTrackId->at(in)!=tmp_out_trk || Event!=tmp_out_Event)){
	 particle_out=true;
	 tmp_out_trk = vTrackId->at(in);
	 tmp_out_Event = Event;
       }
       
       ////// EMU TOP
       if(vVolumeId->at(in)==1 && (vTrackId->at(in)!=tmp_trk_Top || Event!=tmp_Event_Top)){
	 //// tracce totali nell'emulsione Top	    
	 start_Top=false;
	 stop_Top=false;
	 ene_dep_trk=0;
	 tot_ene_dep_Top = 0.;
	 sum_step_len=0.;
	 dE_dx=0;
	 pass_Top=false;
       }
       
      if(vVolumeId->at(in)==1){
	tmp_trk_Top = vTrackId->at(in);
	tmp_Event_Top = Event;
	if(!pass_Top){  // primo step in emu Top
	  first_x_pos = vPreHitX->at(in);
	  first_y_pos = vPreHitY->at(in);
	  first_z_pos = vPreHitZ->at(in);
	  first_ene = vEnergyKin->at(in);
	  if(jentry<10)log << "primo_step_Top " << vEnergyKin->at(in) << " " << vEnergyDep->at(in) << " " << vTrackId->at(in) << " " << vStepId->at(in) << endl;
	  pass_Top=true;
	}
	if(pass_Top){  // exception (mi assicuro che sia l'ultimo step)
	  tot_ene_dep_Top += vEnergyDep->at(in);
	  sum_step_len += vStepLength->at(in);
	  /////////////////////////////////////
	  x_pre = vPreHitX->at(in);
	  y_pre = vPreHitY->at(in);
	  z_pre = vPreHitZ->at(in);
	  x_pos = vPostHitX->at(in);
	  y_pos = vPostHitY->at(in);
	  z_pos = vPostHitZ->at(in);
	  ////////////////////////////////////
	  
	  if(jentry<10)log <<"energie "<< vEnergyKin->at(in) << " " << vEnergyDep->at(in) << endl;
	  //	}
	  
	  if((in==vTrackId->size()-1) || vTrackId->at(in+1)!=tmp_trk_Top){
	    last_x_pos = vPostHitX->at(in);
	    last_y_pos = vPostHitY->at(in);
	    last_z_pos = vPostHitZ->at(in);
	    delta_ene = -1;
	  }
	  else{
	    delta_ene= vEnergyKin->at(in)-(vEnergyKin->at(in+1)+vEnergyDep->at(in));
	  }
	}
	//if(in==0 || (vTrackId->at(count_down+1)!=vTrackId->at(count_down)))  cout << Event << " " << trk_Id << " " << step_Id << " " << in << " " << count_down << endl;
      }
                
      if(vVolumeId->at(in)==1 && vStepId->at(in)==1){
	start_Top=true; // tracce che partono in emulsione
	if(jentry<10)log << "start in emu top " << vEnergyKin->at(in) << " " << vEnergyDep->at(in) <<  " " << vTrackId->at(in) << " " << vStepId->at(in) << endl;  
	}
      
      if(vVolumeId->at(in)==1 && vEnergyKin->at(in)==vEnergyDep->at(in)){ // tracce che si fermano in emulsione
	stop_Top=true;
        //log << "stop in emu top  " << vEnergyKin->at(in) << " " << vEnergyDep->at(in) << " " << vTrackId->at(in) << " " << vStepId->at(in) << endl;
	if(start_Top==true && stop_Top==true){ // tracce contenute in emulsione
	  contained_Top=true;
	  passed_Top=false;
	  last_x_pos = vPostHitX->at(in);
	  last_y_pos = vPostHitY->at(in);
	  last_z_pos = vPostHitZ->at(in);
	  if(jentry<10)log << "start and stop in top  "<< Event << " " << vPartId->at(in) << " " << vTrackId->at(in) << " " << vStepId->at(in) << endl;
	}
      }

      if(in==0)tmp_primary=vTrackId->at(in);
      if(in>0){
	if(vTrackId->at(in)!=vTrackId->at(in-1) && vParentId->at(in)!=tmp_primary){
	   long_primary=false;
	   long_secondary=false;
	   if(sum_step_len>100)tmp_primary=vTrackId->at(in);
	}

	
	//if((vTrackId->at(in)==vTrackId->at(in-1) && sum_step_len>100) || vParentId->at(in)!=tmp_primary)long_primary=true;
	//if(long_primary && sum_step_len>100)long_secondary=true;
	//if(vParentId->at(in)==tmp_primary && (sum_step_len<100 || (sum_step_len-step_len)>100))long_secondary=false;
	//if(long_primary && long_secondary)multiplicity++;
	//if(jentry<10) cout << Event << " " << vParentId->at(in) << " " << vTrackId->at(in) << " " << vStepId->at(in) << " " << tmp_primary <<" " << sum_step_len << " " << long_primary << " " << long_secondary << " " << multiplicity <<  endl; 
      } 
      /////////////////////

      
      ////// FOR CLUSTERING   //////////////
      //if(stop_Top && TMath::Abs(vPostHitX->at(in)+13250)<2000 && vPostHitY->at(in)>-3022.5 && vPostHitY->at(in)<-1877.5 && vPostHitZ->at(in) <-488){ // PN1
      //if(stop_Top && TMath::Abs(vPostHitY->at(in)-1000)<2000 && vPostHitX->at(in)>-13272.5 && vPostHitX->at(in)<-11227.5 && vPostHitZ->at(in) <-488){  // PN2
      //if(stop_Top && TMath::Abs(vPostHitX->at(in))<7250 && TMath::Abs(vPostHitY->at(in))<7250 && TMath::Abs(vPostHitZ->at(in))<7250 ){  // C14 ten grams
      if(stop_Top /*&& TMath::Abs(vPostHitX->at(in))<7250 && TMath::Abs(vPostHitY->at(in))<7250 && TMath::Abs(vPostHitZ->at(in))<7250*/){  // gamma_env
	ev_x_pos.push_back(vPostHitX->at(in));
	ev_y_pos.push_back(vPostHitY->at(in));
	ev_z_pos.push_back(vPostHitZ->at(in));
	ev_ene_dep.push_back(vEnergyDep->at(in));
	//cout << jentry << " " << iel << " " << ev_x_pos[iel] << " " << ev_y_pos[iel] << " " << ev_y_pos[iel] << endl;
	iel++;
	}
     
      
      /////////////////////
      if(contained_Top || passed_Top)all_trk=true;
      if(vVolumeId->at(in)==1)ene_dep_trk = tot_ene_dep_Top;
      //if(vVolumeId->at(in)==2)ene_dep_trk = tot_ene_dep_Down;
      trk_Rec(first_x_pos,first_y_pos,first_z_pos,last_x_pos,last_y_pos,last_z_pos,trk_rec);
      trk_len = trk_rec[0];
      trk_the = trk_rec[1];
      trk_phi = trk_rec[2];
      trk_ene = vEnergyKin->at(in);

      //last_step_len = TMath::Sqrt(TMath::Power(last_x_pos - hit_x_pre_last,2)+TMath::Power(last_y_pos - hit_y_pre_last,2)+TMath::Power(last_z_pos - hit_z_pre_last,2));
      //step_len = vStepLength->at(in);
      //trk_poly_len = vTrackLength->at(in); 

      if(trk_len==0){
	trk_the=-5;
	trk_phi=-5;
	trk_ene=-5;
      }
      ///////////////////////////////
      if(jentry<10){
	log <<"Tree value " << ev_Id << " " << parent_Id << " " << trk_Id << " " << step_Id << " " << emu_Id << " " << pdg_Id << " " << ene_dep_trk << " " << trk_len << " " << trk_the << " " << trk_phi << endl;
      log <<"Tree bool "<<  passed_Top << " " << start_Top << " " << stop_Top << " "  << contained_Top << " " << trk_father << endl;
      log <<"first position "<< first_x_pos << " " << first_y_pos << " " << first_z_pos << endl;
      log <<"last position "<< last_x_pos << " " << last_y_pos << " " << last_z_pos << endl;
      log <<"last step length "<< step_len << " " << trk_poly_len << " " << first_ene << endl;
      log << "energies " << trk_ene << " " << ene_dep_step << " " << delta_ene << endl;
      log << endl;
      }
      //Tree_out->Fill();
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

     //cout << "icl " << jentry << " " << icl << endl;
     
     /////// GRAPHICS

     //if(jentry==70948)cout << "num cl " << jentry <<" " <<  icl << endl;

     
     ev_x_pos.clear();
     ev_y_pos.clear();
     ev_z_pos.clear();
     ev_ene_dep.clear();
     
     //_grafico_TCanvas *cc = new TCanvas();
     //_grafico_TH3F *gr[icl][jentry];

     for(int i=0;i<icl;i++){
       clx_mean=0;
       cly_mean=0;
       clz_mean=0;
       cl_dep_ene=0;

       //cout <<"gr nel "<<jentry << " " <<  i << " " << n_el[i] << endl;
       for(int j=0;j<n_el[i];j++){
	 //cout << jentry << " " << icl << " " << n_el[i] << endl;
	 cl_dep_ene += cl_ene_dep[i][j];
	 //cout <<"gr pt " <<  i << " " << j << " " <<  cl_x[i][j] << " " << cl_y[i][j] << " " << cl_z[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 //cout <<"gr pt abs " <<  i << " " << j << " " <<  cl_x_abs[i][j] << " " << cl_y_abs[i][j] << " " << cl_z_abs[i][j] << " " << cl_ene_dep[i][j] << " " << cl_dep_ene << endl;
	 clx_mean += cl_x_abs[i][j]*cl_ene_dep[i][j];
	 cly_mean += cl_y_abs[i][j]*cl_ene_dep[i][j];
	 clz_mean += cl_z_abs[i][j]*cl_ene_dep[i][j];
	 	
       }
       clx_mean /= cl_dep_ene;
       cly_mean /= cl_dep_ene;
       clz_mean /= cl_dep_ene;
       ev_id = jentry;
       //cl_id = icl;
       ntot = n_el[i];
       //cout << clx_mean << " " << cly_mean << " " << clz_mean << endl;
       //// Cut energy on small clusters 
       //if(TMath::Abs(clx_mean+13250)<1000 && cly_mean>-2022.5 && cly_mean<-1977.5 && clz_mean <-488 && cl_dep_ene>thr_ene){  // PN1
       //if(TMath::Abs(cly_mean-1000)<1000 && clx_mean>-12272.5 && clx_mean<-12227.5 && clz_mean <-488 && cl_dep_ene>thr_ene){   // PN2
       //if(TMath::Abs(clx_mean)<7150 && TMath::Abs(cly_mean)<7150 && TMath::Abs(clz_mean)<7150 && cl_dep_ene>thr_ene){   // C14 ten grams
       if(/*TMath::Abs(clx_mean)<7150 && TMath::Abs(cly_mean)<7150 && TMath::Abs(clz_mean)<7150 &&*/ cl_dep_ene>thr_ene){   // gamma env
	 ev_x_pos2.push_back(clx_mean);
	 ev_y_pos2.push_back(cly_mean);
	 ev_z_pos2.push_back(clz_mean);
	 ev_ene_dep2.push_back(cl_dep_ene);
	 //cout << "eccomi " << jentry << endl;
	 cl_id = tot_cl;
	 tot_cl++;
	 //Tree_cl_small->Fill();
       }
       //if(jentry==70948) cout << "mean " << clx_mean << " " << cly_mean << " " << clz_mean << " " << cl_dep_ene << endl;
      
       
       //_grafico_if(i==0)	 gr[i][jentry]->Draw("");       
       //_grafico_else gr[i][jentry]->Draw("SAMES");
       
     }     
   } /// end jentry
   cout << "small clusters: " << tot_cl <<", ncl_ev: " << icl << ", epsilon1: " << epsilon1 << ", thr_ene: "<< thr_ene << endl;
   //////////////////////////////////////////////

   //////////// BIG CLUSTERS //////////////////////////////////
   
   //// CLUSTERING PER EVENTO
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
     cl_x_abs.clear();
     cl_y_abs.clear();
     cl_z_abs.clear();
     visited.clear();
     visited2.clear();
     visited.resize(tot_cl);
     visited2.resize(tot_cl);
     n_el.resize(tot_cl);
     id_el2.resize(tot_cl);

     /* cl_el.resize(tot_cl, std::vector<double>(tot_cl));
     cl_x.resize(tot_cl, std::vector<double>(tot_cl));
     cl_y.resize(tot_cl, std::vector<double>(tot_cl));
     cl_z.resize(tot_cl, std::vector<double>(tot_cl));
     cl_ene_dep.resize(tot_cl, std::vector<double>(tot_cl));
     cl_x_abs.resize(tot_cl, std::vector<double>(tot_cl));
     cl_y_abs.resize(tot_cl, std::vector<double>(tot_cl));
     cl_z_abs.resize(tot_cl, std::vector<double>(tot_cl));
     */

     cl_el.resize(tot_cl);
     cl_x.resize(tot_cl);
     cl_y.resize(tot_cl);
     cl_z.resize(tot_cl);
     cl_ene_dep.resize(tot_cl);
     cl_x_abs.resize(tot_cl);
     cl_y_abs.resize(tot_cl);
     cl_z_abs.resize(tot_cl);




     //cout << visited.size() << " " << visited2.size() << endl;
		      
     
     for(int j=0;j<tot_cl;j++){
       //cout << j << endl;
       if(visited.at(j)==false){
	 sphere_points=0;
	 sphere_points2=0;
	 visited[j]=true;
	 cl_el[icl].push_back(j);
	 if(j==0){
	   tmp_x_ev=ev_x_pos2.at(j);
	   tmp_y_ev=ev_y_pos2.at(j);
	   tmp_z_ev=ev_z_pos2.at(j);
	 }
	 /*cl_x[icl][sphere_points]=ev_x_pos2.at(j)-tmp_x_ev;
	 cl_y[icl][sphere_points]=ev_y_pos2.at(j)-tmp_y_ev;
	 cl_z[icl][sphere_points]=ev_z_pos2.at(j)-tmp_z_ev;
	 cl_ene_dep[icl][sphere_points]=ev_ene_dep2.at(j);
	 cl_x_abs[icl][sphere_points]=ev_x_pos2.at(j);
	 cl_y_abs[icl][sphere_points]=ev_y_pos2.at(j);
	 cl_z_abs[icl][sphere_points]=ev_z_pos2.at(j);
	 */
	 cl_x[icl].push_back(ev_x_pos2.at(j)-tmp_x_ev);
	 cl_y[icl].push_back(ev_y_pos2.at(j)-tmp_y_ev);
	 cl_z[icl].push_back(ev_z_pos2.at(j)-tmp_z_ev);
	 cl_ene_dep[icl].push_back(ev_ene_dep2.at(j));
	 cl_x_abs[icl].push_back(ev_x_pos2.at(j));
	 cl_y_abs[icl].push_back(ev_y_pos2.at(j));
	 cl_z_abs[icl].push_back(ev_z_pos2.at(j));
	 
	 sphere_points++;
	 for(int k=0;k<tot_cl;k++){
	   //dist=TMath::Sqrt(TMath::Power(ev_x_pos2[k]-ev_x_pos2[j],2)+TMath::Power(ev_y_pos2[k]-ev_y_pos2[j],2)+TMath::Power(ev_z_pos2[k]-ev_z_pos2[j],2));
	   dist=TMath::Sqrt(TMath::Power(ev_x_pos2.at(k)-ev_x_pos2.at(j),2)+TMath::Power(ev_y_pos2.at(k)-ev_y_pos2.at(j),2)+TMath::Power(ev_z_pos2.at(k)-ev_z_pos2.at(j),2));
	   if(dist<epsilon2 && k!=j && !visited2.at(j)){
	     visited2[k]=true;
	     //cout <<"a "<< icl << " " <<  j << " " << k << " " << sphere_points << " " << dist <<  endl;
	     /*cl_el[icl][sphere_points]=k;	    
	     cl_x[icl][sphere_points]=ev_x_pos2.at(k)-tmp_x_ev;
	     cl_y[icl][sphere_points]=ev_y_pos2.at(k)-tmp_y_ev;
	     cl_z[icl][sphere_points]=ev_z_pos2.at(k)-tmp_z_ev;
	     cl_ene_dep[icl][sphere_points]=ev_ene_dep2.at(k);
	     cl_x_abs[icl][sphere_points]=ev_x_pos2.at(k);
	     cl_y_abs[icl][sphere_points]=ev_y_pos2.at(k);
	     cl_z_abs[icl][sphere_points]=ev_z_pos2.at(k);*/
	     cl_el[icl].push_back(k);
	     cl_x[icl].push_back(ev_x_pos2.at(k)-tmp_x_ev);
	     cl_y[icl].push_back(ev_y_pos2.at(k)-tmp_y_ev);
	     cl_z[icl].push_back(ev_z_pos2.at(k)-tmp_z_ev);
	     cl_ene_dep[icl].push_back(ev_ene_dep2.at(k));
	     cl_x_abs[icl].push_back(ev_x_pos2.at(k));
	     cl_y_abs[icl].push_back(ev_y_pos2.at(k));
	     cl_z_abs[icl].push_back(ev_z_pos2.at(k));
	     
	     sphere_points++;
	   }
	 }
	 
	 if(sphere_points < (min_points+1) && !visited2.at(j)) {
	   //cout <<"b "<< icl << " " <<  j << " " <<  sphere_points <<  endl;
	   /*cl_el[icl][sphere_points]=j;
	   cl_x[icl][sphere_points]=ev_x_pos2.at(j)-tmp_x_ev;
	   cl_y[icl][sphere_points]=ev_y_pos2.at(j)-tmp_y_ev;
	   cl_z[icl][sphere_points]=ev_z_pos2.at(j)-tmp_z_ev;
	   cl_ene_dep[icl][sphere_points]=ev_ene_dep2.at(j);
	   cl_x_abs[icl][sphere_points]=ev_x_pos2.at(j);
	   cl_y_abs[icl][sphere_points]=ev_y_pos2.at(j);
	   cl_z_abs[icl][sphere_points]=ev_z_pos2.at(j);*/
	   cl_el[icl].push_back(j);
	   cl_x[icl].push_back(ev_x_pos2.at(j)-tmp_x_ev);
	   cl_y[icl].push_back(ev_y_pos2.at(j)-tmp_y_ev);
	   cl_z[icl].push_back(ev_z_pos2.at(j)-tmp_z_ev);
	   cl_ene_dep[icl].push_back(ev_ene_dep2.at(j));
	   cl_x_abs[icl].push_back(ev_x_pos2.at(j));
	   cl_y_abs[icl].push_back(ev_y_pos2.at(j));
	   cl_z_abs[icl].push_back(ev_z_pos2.at(j));
	   n_el[icl]=1;
	   //cout <<"cl_el2 "<< icl << " " << n_el[icl] << " " << cl_el[icl][sphere_points] << endl;
	   //cout << "c " << icl << " " << sphere_points << " " <<  cl_x[icl][sphere_points] << " " << cl_y[icl][sphere_points] << endl;
	   //cout << "big " << icl << " " << cl_el[icl][0] << " " << sphere_points << " " << sphere_points2 << endl;
	   icl++;	   
	   continue;
	 }
	 else {
	   for(int j=0;j<sphere_points;j++){
	     //sphere_points2=0;
	     if(visited.at(cl_el[icl][j])==false){
	       visited[cl_el[icl][j]]=true;
	       for(int k=0;k<tot_cl;k++){
		 // dist = (TMath::Sqrt(TMath::Power(ev_x_pos2[k]-ev_x_pos2[cl_el[icl][j]],2)+TMath::Power(ev_y_pos2[k]-ev_y_pos2[cl_el[icl][j]],2)
		 //		     + TMath::Power(ev_z_pos2[k]-ev_z_pos2[cl_el[icl][j]],2)));
		 dist = (TMath::Sqrt(TMath::Power(ev_x_pos2.at(k)-ev_x_pos2.at(cl_el[icl][j]),2)+TMath::Power(ev_y_pos2.at(k)-ev_y_pos2.at(cl_el[icl][j]),2)
				     + TMath::Power(ev_z_pos2.at(k)-ev_z_pos2.at(cl_el[icl][j]),2)));
		 if(dist<epsilon2 && k!=cl_el[icl][j] && !visited2.at(k) && !visited.at(k)){
		   id_el2[sphere_points2]=k;
		   visited2[k]=true;
		   //cout <<"c1 "<< icl << " " <<  j << " " << cl_el[icl][j] << " " << k << " " << sphere_points2 <<  " " << dist << " " << id_el2[sphere_points2] << endl;	
		   sphere_points2++;
		 }
	       }
	       if(sphere_points2 >= min_points){
		 for(int i=0;i<sphere_points2;i++){
		   /*cl_el[icl][sphere_points+i]=id_el2[i];
		     cl_x[icl][sphere_points+i]=ev_x_pos2.at(id_el2[i])-tmp_x_ev;
		     cl_y[icl][sphere_points+i]=ev_y_pos2.at(id_el2[i])-tmp_y_ev;
		     cl_z[icl][sphere_points+i]=ev_z_pos2.at(id_el2[i])-tmp_z_ev;
		     cl_ene_dep[icl][sphere_points+i]=ev_ene_dep2.at(id_el2[i]);
		     cl_x_abs[icl][sphere_points+i]=ev_x_pos2.at(id_el2[i]);
		     cl_y_abs[icl][sphere_points+i]=ev_y_pos2.at(id_el2[i]);
		     cl_z_abs[icl][sphere_points+i]=ev_z_pos2.at(id_el2[i]);*/
		   cl_el[icl].push_back(id_el2[i]);
		   cl_x[icl].push_back(ev_x_pos2.at(id_el2[i])-tmp_x_ev);
		   cl_y[icl].push_back(ev_y_pos2.at(id_el2[i])-tmp_y_ev);
		   cl_z[icl].push_back(ev_z_pos2.at(id_el2[i])-tmp_z_ev);
		   cl_ene_dep[icl].push_back(ev_ene_dep2.at(id_el2[i]));
		   cl_x_abs[icl].push_back(ev_x_pos2.at(id_el2[i]));
		   cl_y_abs[icl].push_back(ev_y_pos2.at(id_el2[i]));
		   cl_z_abs[icl].push_back(ev_z_pos2.at(id_el2[i]));
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
	     //cout << "big " << icl << " " << i << " " << cl_el[icl][i] << " " << sphere_points << " " << sphere_points2 << endl;	     
	     n_el[icl]++;
	   }
	   //cout << endl;
	   //cout <<"n_el "<< icl << " " << n_el[icl] << endl;
	 }
	 //cout <<"f " <<  j << " " << visited[j] << " " << visited2[j] << " " << sphere_points << endl;
	 if(!visited2.at(j) || sphere_points<(min_points)){
	   //for(int i=0;i<(sphere_points+sphere_points2);i++){
	   //cout << "big " << icl << " " << i << " " << cl_el[icl][i] << " " << sphere_points << " " << sphere_points2 << endl;	     }
	   icl++;
	 }
       }
       // cout << j << " " << icl << " " << ev_x_pos2[j] << " " << ev_y_pos2[j] << " " << ev_z_pos2[j] << endl;
     } /// end clustering
     //cout << "ntot cl2 " << icl << endl;
     cout << "big clusters: " << icl <<  ", epsilon2: " << epsilon2 <<  endl;
     //////////////////////////////////////////////
     
    /////// 
   
     for(int i=0;i<icl;i++){
       //clx_mean=0;
       //cly_mean=0;
       //clz_mean=0;
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
	 //clx_mean += cl_x_abs[i][j]*cl_ene_dep[i][j];
	 //cly_mean += cl_y_abs[i][j]*cl_ene_dep[i][j];
	 //clz_mean += cl_z_abs[i][j]*cl_ene_dep[i][j];
	 cl_id_small.push_back(cl_el[i][j]);
	 clx_small.push_back(cl_x_abs[i][j]);
	 cly_small.push_back(cl_y_abs[i][j]);
	 clz_small.push_back(cl_z_abs[i][j]);
	 cl_ene_small.push_back(cl_ene_dep[i][j]);
       }
       //clx_mean /= cl_dep_ene;
       //cly_mean /= cl_dep_ene;
       //clz_mean /= cl_dep_ene;
       
       cl_id_big = i;
       Ntot = n_el[i];
       Tree_cl_big->Fill();       
       }     
   //////////////////////////////////////////////
   


   log.close();
   //Tree_out->Write();
   //Tree_cl_small->Write();
   Tree_cl_big->Write();
   
   TCanvas *CEne_MC = new TCanvas("CEne_MC","",600,600);
   h_Ene_MC->Draw();
   CEne_MC->SetLogy();
   h_Ene_MC->SetLineWidth(2);
   
   double max_phi_MC = h_Phi_MC->GetMaximum()*1.2;
   TCanvas *CPhi_MC = new TCanvas("CPhi_MC","",600,600);
   h_Phi_MC->Draw();
   h_Phi_MC->SetLineWidth(2);
   h_Phi_MC->GetYaxis()->SetRangeUser(0,max_phi_MC);
   
   TCanvas *CTheta_MC = new TCanvas("CTheta_MC","",600,600);
   h_Theta_MC->Draw();
   h_Theta_MC->SetLineWidth(2);
   
   TCanvas *CCosTheta_MC = new TCanvas("CCosTheta_MC","",600,600);
   h_CosTheta_MC->Draw();
   h_CosTheta_MC->SetLineWidth(2);
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
