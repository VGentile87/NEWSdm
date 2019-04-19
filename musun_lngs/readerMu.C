#include "Riostream.h"
#include <math.h>
void readerMu() {
//  Read data from an ascii file and create a root file with an histogram and an ntuple.
//   see a variant of this macro in basic2.C
//Author: Rene Brun


// read file $ROOTSYS/tutorials/tree/basic.dat
// this file has 3 columns of float data
   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("readerMu.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   //in.open(Form("%smuons_gs_01.dat",dir.Data()));
   in.open(Form("%smuons_gs_seed_100.dat",dir.Data()));

   const int nbins=100;
   Double_t xbins[nbins+1];
   double dx = 6./nbins;
   double l10 = TMath::Log(10);
   for (int i=0;i<=nbins;i++) {
     xbins[i] = TMath::Exp(l10*i*dx);
   }
   
   Float_t pi = TMath::Pi();

   Float_t ev,id,ene,x,y,z,cx,cy,cz;
   Int_t nlines = 0;
   TFile *f = new TFile("muon_data.root","RECREATE");
   TH3F *h_pos = new TH3F("h_pos","x distribution",200,-4000,4000,200,-4000,4000,200,-4000,4000);
   //TH3F *h_pos = new TH3F("h_pos","x distribution",200,-0.5,0.5,200,-0.5,0.5,200,-0.5,0.5);
   TH1F *h_theta = new TH1F("h_theta","x distribution",100,0,pi);
   TH1F *h_cos_theta = new TH1F("h_cos_theta","x distribution",100,0,1);
   TH1F *h_domega = new TH1F("h_domega","x distribution",100,0,2*pi);
   TH1F *h_domega_2 = new TH1F("h_domega_2","x distribution",100,0,1);
   TH1F *h_phi = new TH1F("h_phi","x distribution",100,-pi,pi);
   TH1F *h_ene = new TH1F("h_ene","x distribution",nbins,xbins);
   TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","ev:id:ene:x:y:z:cx:cy:cz");

   Float_t phi=0;
   Float_t theta=0;
   Float_t c_mod=0;

   const Float_t glob_int = 0.5870;
   //const Float_t glob_int = 0.3256*TMath::Power(10,-7);


   do{
     in >> ev >> id >> ene >> x >> y >> z >> cx >> cy >> cz;
     if (!in.good()) break;
     cout << ev << endl;
     if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
     c_mod = TMath::Sqrt(cx*cx+cy*cy+cz*cz);
     phi = TMath::ATan2(cy,cx);
     //theta = TMath::ACos(-cz/c_mod);
     theta = TMath::ACos(-cz);
     int theta_bin = pi/100. + (pi/50.)*floor(theta/(pi/50.));
     //cout << theta << " " << theta_bin << endl;
     h_pos->Fill(x,y,z);
     h_phi->Fill(phi);
     h_theta->Fill(theta);
     h_cos_theta->Fill(TMath::Cos(theta),2*pi*TMath::Sin(theta));
     h_ene->Fill(ene);
     h_domega->Fill((2*pi*TMath::Sin(theta)));
     h_domega_2->Fill((TMath::Sin(theta)));
    
     ntuple->Fill(ev,id,ene,x,y,z,cx,cy,cz);
     nlines++;
   }while(nlines<1000000);
   printf(" found %d points\n",nlines);

   in.close();

   Double_t scale = glob_int/h_phi->Integral();
   //h_phi->Scale(scale);
   h_theta->Scale(scale);
   //h_cos_theta->Scale(scale);
   h_domega->Scale(scale);
   //h_ene->Scale(scale);

   //TH1F * h_cos_theta2;
   //h_cos_theta2->SetName("h_cos_theta2");
   //h_cos_theta2 = (TH1F*)h_cos_theta->Clone("h_cos_theta2");
   //h_cos_theta2->Divide(h_domega_2);
   //h_cos_theta2->Scale(2);
   
   TCanvas *c1 = new TCanvas("c1");
   h_cos_theta->Draw("L");

   TCanvas *c3 = new TCanvas("c3");
   h_phi->Draw("L");

   TCanvas *c2 = new TCanvas("c2");
   c2->SetLogx();
   h_ene->Draw("L");   

   f->Write();
}
