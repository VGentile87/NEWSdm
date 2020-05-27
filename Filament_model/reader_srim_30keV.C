#define reader_srim_30keV_cxx
#include "reader_srim_30keV.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void reader_srim_30keV::Loop(){
  
  if (fChain == 0) return;
  
  Long64_t nentries = fChain->GetEntriesFast();

  float step;
  // track position in SRIM
  vector<float> xhit;
  vector<float> yhit;
  vector<float> zhit;
  vector<float> enedep;
  vector<float> stopping;
  vector<float> kinene;

  //TH2F *h2d = new TH2F("h2d","h2d",120,0,30,120,0,3);
  TProfile *h2d = new TProfile("h2d","h2d",300,0,30,0,3);
  
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
    GetEntry(jentry);
    Long64_t ientry = LoadTree(jentry);
    if(ientry < 0) break;
    if(jentry<10)cout << ev << " " << x << " " << rene << endl;
    
    xhit.push_back(z);
    yhit.push_back(y);
    zhit.push_back(x);
    enedep.push_back(dene);
    stopping.push_back(rene);
    kinene.push_back(kene);

    if(rene==0){
      
      step=0;
      xhit.push_back(x);
      yhit.push_back(y);
      zhit.push_back(z);
      enedep.push_back(dene);
      kinene.push_back(kene);
      
      
      for(int i=1;i<xhit.size();i++){
	step=0;
	TVector3 p1(xhit.at(i-1),yhit.at(i-1),zhit.at(i-1));
	TVector3 p2(xhit.at(i),yhit.at(i),zhit.at(i));
	TVector3 dp;
	dp = p2-p1;
	step = dp.Mag();
	h2d->Fill(kinene.at(i),enedep.at(i)/step);
	//h2d->Fill(2,2);
	cout << i << " " << enedep.at(i) << " " << kinene.at(i) << " " << step << endl;
      } // end of the track
      
      xhit.clear();
      yhit.clear();
      zhit.clear();
      enedep.clear();
      stopping.clear();
      kinene.clear();
      
    }
  }

  h2d->SaveAs("h2d.root");
  //TCanvas *c1 = new TCanvas();
  //h2d->Draw("colz");
}

void arun(){
  // gSystem->Load(".L filament_model.C");
  reader_srim_30keV srim;
  srim.Loop();
  // return 0;
}
