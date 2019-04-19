#include "Riostream.h"
void reader_SRIM() {
//  Read data from an ascii file and create a root file with an histogram and an ntuple.
//   see a variant of this macro in basic2.C
//Author: Rene Brun


// read file $ROOTSYS/tutorials/tree/basic.dat
// this file has 3 columns of float data
   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("reader_SRIM.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open("EXYZ.txt");

   Int_t ev;
   Float_t x,y,z, kene, dene, rene;
   Float_t tmp_x,tmp_y,tmp_z;
   Int_t nlines = 0;
   TFile *f = new TFile("SRIM.root","RECREATE");
   //TH1F *h1 = new TH1F("h1","x distribution",100,-4,4);
   TNtuple *ntuple = new TNtuple("ntuple","data from ascii file","ev:kene:x:y:z:rene:dene");

   while (1) {
     cout << "hello" << endl;
     in >> ev >> kene >> x >> y >> z >> rene >> dene;
     cout << ev <<" " << kene << " " <<  x << " " << y << " " << z << " " << rene << " " << dene <<   endl;
     if(kene!=0 && rene!=0){
       tmp_x = x;
       tmp_y = y;
       tmp_z = z;
     }
     if(kene!=0 && rene==0){
       x = tmp_x;
       y = tmp_y;
       z = tmp_z;
     }
     if (!in.good()) {
       cout << "break" << endl;
       break;
     }
     //if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
     //h1->Fill(x);
     //if(tmp_rene<rene)
     //tmp_rene = rene;
     ntuple->Fill(ev,kene,x/10.,y/10.,z/10.,rene,dene/1000.);
     nlines++;
   }
   printf(" found %d points\n",nlines);

   in.close();

   f->Write();
}
