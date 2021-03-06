#include "Riostream.h"
void MakeMuons() {
//  Read data from an MUSUN.dat file and create a hepMC file.


  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("MakeMuons.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%smuons_gs_seed_101.dat",dir.Data()));
  
  Float_t ev,id,ene,x,y,z,cx,cy,cz;
  //Float_t px=0, py=0, pz=0;
  Int_t nlines = 0;
  
  Int_t pdg=0;
  
  
    FILE *log = fopen("hepmuon_1M_seed_101.txt","w+");
    
    fprintf(log,"HepMC::Version 2.06.09\n");
    fprintf(log,"HepMC::IO_GenEvent-START_EVENT_LISTING\n");
    
    do{
      in >> ev >> id >> ene >> x >> y >> z >> cx >> cy >> cz;
      if (!in.good()) break;
      if(id==10)pdg=-13;
      if(id==11)pdg=13;
      // cout << ev << " " << pdg << endl;
      if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
      
      fprintf(log,"E %d  -1 -1.0 -1.0 -1.0 0 0 1 1 2 0 0\n",ev-1);
      fprintf(log,"U GEV MM\n");
      fprintf(log,"V -1 0 %f %f %f 0 0 1 0\n",x*10,y*10,z*10);
      fprintf(log,"P 1 %d %f %f %f %f %f 1 0 0 0 0\n",pdg,ene,ene,cx,cy,cz);
      
      nlines++;
    }while(nlines<100000)
      printf(" found %d points\n",nlines);

  in.close();

  fprintf(log,"HepMC::IO_GenEvent-END_EVENT_LISTING\n");
  fclose(log);
}


  

/*
  Int_t N = 100;
  Float_t px=0, py=0, pz=0;
  
  //valori da impostare
  Int_t pdg        =13;
  Float_t px_mean  = 0.;
  Float_t py_mean  = 0.;
  Float_t pz_mean  =10.;
  Float_t px_sigma = 0.;
  Float_t py_sigma = 0.;
  Float_t pz_sigma = 0.;
  
  TRandom3 r;
  
  fprintf(f,"HepMC::Version 2.06.09\n");
  fprintf(f,"HepMC::IO_GenEvent-START_EVENT_LISTING\n");

  for(Int_t i=0; i<N; i++){
    px=r.Gaus(px_mean,px_sigma);
    py=r.Gaus(py_mean,py_sigma);
    pz=r.Gaus(pz_mean,pz_sigma);
    fprintf(f,"E %d  -1 -1.0 -1.0 -1.0 0 0 1 1 2 0 0\n",i);
    fprintf(f,"U GEV MM\n");
    fprintf(f,"V -1 0 0.000000 0.000000 0.000000 0 0 1 0\n");
    fprintf(f,"P 1 %d %f %f %f 0.105658 0.105658 1 0 0 0 0\n",pdg,px,py,pz);
  }
  fprintf(f,"HepMC::IO_GenEvent-END_EVENT_LISTING\n");
  
  fclose(f);
  
}
*/
