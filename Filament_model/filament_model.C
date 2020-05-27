#include "Riostream.h"
#include <math.h>
//#include <tuple>
#define SRIM_data_cxx
#include "SRIM_data.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TVector3.h>


const float pi = TMath::Pi();
float  unit_cube = 2; //2nm

bool RayBox(TVector3 bc, TVector3 p1, TVector3 p2) {

  //bool inBox=false;

  float unit_length=2;
  
  float minx = bc.X()-unit_length;
  float maxx = bc.X()+unit_length;
  float miny = bc.Y()-unit_length;
  float maxy = bc.Y()+unit_length;
  float minz = bc.Z()-unit_length;
  float maxz = bc.Z()+unit_length;
  
  TVector3 dp = p2-p1;
  TVector3 dir = dp.Unit();
  
  float tmin = (minx - p1.X()) / dir.X(); 
  float tmax = (maxx - p1.X()) / dir.X(); 
  
  if (tmin > tmax) swap(tmin, tmax); 
  
  float tymin = (miny - p1.Y()) / dir.Y(); 
  float tymax = (maxy - p1.Y()) / dir.Y(); 
 
  if (tymin > tymax) swap(tymin, tymax); 
  
  if ((tmin > tymax) || (tymin > tmax)) 
  return false; 
  
  if (tymin > tmin) 
    tmin = tymin; 
 
  if (tymax < tmax) 
    tmax = tymax; 
  
  float tzmin = (minz - p1.Z()) / dir.Z(); 
  float tzmax = (maxz - p1.Z()) / dir.Z(); 
 
  if (tzmin > tzmax) swap(tzmin, tzmax); 
 
  if ((tmin > tzmax) || (tzmin > tmax)) 
  return false; 
  
  if (tzmin > tmin) 
    tmin = tzmin; 
  
  if (tzmax < tmax) 
    tmax = tzmax; 
  
  return true; 
}


tuple<vector<int>,vector<int>,vector<int>, vector<int>, vector<int>, vector<int>> MakeGrid(int event, float unit_length, vector<float> xcr, vector<float> ycr, vector<float> zcr, vector<float> rcr, vector<float> xfil, vector<float> yfil, vector<float> zfil, vector<float> slen){

  TH3F *hbox = new TH3F("box","box",500,-2200,-1200,500,-2200,-1200,500,0,1000);
  TH3F *hgrid = new TH3F("grid","grid",500,-2200,-1200,500,-2200,-1200,500,0,1000);

  vector<int> inx;
  vector<int> iny;
  vector<int> inz;
  vector<int> outx;
  vector<int> outy;
  vector<int> outz;
  
  for(int j=0;j<1/*rcr.size()*/;j++){    
    // cell
    int bx = hbox->GetXaxis()->FindBin(xcr.at(j));
    int by = hbox->GetYaxis()->FindBin(ycr.at(j));
    int bz = hbox->GetZaxis()->FindBin(zcr.at(j));
    int br = round(rcr.at(j));
    int ix = bx-br;
    int iy = by-br;
    int iz = bz-br;
    int fx = bx+br;
    int fy = by+br;
    int fz = bz+br;
    // cout << "ix " << ix << " " << iy << " " << iz << " " << fx << " " << fy << " " << fz << " " << br << endl;
    //cout << "ix2 " << bx << " " << by << " " << bz << " " << xcr.at(j) << " " << ycr.at(j) << " " << zcr.at(j) << endl;
    if(ix<0)ix=0;
    if(iy<0)iy=0;
    if(iz<0)iz=0;
    if(fx>hbox->GetNbinsX())fx=hbox->GetNbinsX();
    if(fy>hbox->GetNbinsY())fy=hbox->GetNbinsY();
    if(fz>hbox->GetNbinsZ())fz=hbox->GetNbinsZ();
    
    for(int ib=ix;ib<fx;ib++){
      for(int jb=iy;jb<fy;jb++){
	for(int kb=iz;kb<fz;kb++){
	  hbox->SetBinContent(ib,jb,kb,1);
	  //inx.push_back(ib);
	  //iny.push_back(jb);
	  //inz.push_back(kb);
	}
      }
    }  
  }

  

  // GRID MULTIPLE DIPOLE
  for(int s=0;s<9/*xfil.size()*/;s++){
    if(slen.at(s)!=0){
      int xin = hbox->GetXaxis()->FindBin(xfil.at(s-1));
      int xout = hbox->GetXaxis()->FindBin(xfil.at(s));
      int yin = hbox->GetYaxis()->FindBin(yfil.at(s-1));
      int yout = hbox->GetYaxis()->FindBin(yfil.at(s));
      int zin = hbox->GetZaxis()->FindBin(zfil.at(s-1));
      int zout = hbox->GetZaxis()->FindBin(zfil.at(s));

      int minx,miny,minz,maxx,maxy,maxz;

      TVector3 p1(xfil.at(s-1),yfil.at(s-1),zfil.at(s-1));
      TVector3 p2(xfil.at(s),yfil.at(s),zfil.at(s));

      if(xin<=xout){
	minx=xin;
	maxx=xout;
      }
      else {
	minx=xout;
	maxx=xin;
      }

      if(yin<=yout){
	miny=yin;
	maxy=yout;
      }
      else {
	miny=yout;
	maxy=yin;
      } 

      if(zin<=zout){
	minz=zin;
	maxz=zout;
      }
      else {
	minz=zout;
	maxz=zin;
      } 

      cout << "s " << s << endl;
      cout << "grid " << xin << " " << xout << " " << yin << " " << yout << " " << zin << " " << zout << endl;
      cout << "max " << minx << " " << maxx << " " << miny << " " << maxy << " " << minz << " " << maxz << endl;
      
      for(int i=minx;i<(maxx+1);i++){
	for(int j=miny;j<(maxy+1);j++){
	  for(int k=minz;k<(maxz+1);k++){
	    if(hbox->GetBinContent(i,j,k)==1){
	      
	      TVector3 sc(-2200+unit_length*i -unit_length/2., -2200+unit_length*j -unit_length/2., 0+unit_length*k -unit_length/2.);
	      //cout << "sc " << sc(0) << " " << sc(1) << " " << sc(2) << endl;	      
	      bool inBox = RayBox(sc,p1,p2);
	      
	      if(inBox){
		hgrid->SetBinContent(i,j,k,1);
		inx.push_back(i);
		iny.push_back(j);
		inz.push_back(k);
	      }
	    }
	    
	    if(hbox->GetBinContent(i,j,k)!=1){

	      TVector3 sc(-2200+unit_length*i -unit_length/2., -2200+unit_length*j -unit_length/2., 0+unit_length*k -unit_length/2.);
	      //cout << "sc " << sc(0) << " " << sc(1) << " " << sc(2) << endl;	      
	      bool inBox = RayBox(sc,p1,p2);
	      
	      if(inBox){
		outx.push_back(i);
		outy.push_back(j);
		outz.push_back(k);
	      }
	    }
	  }
	}
      }   
    }
  }
  
  hgrid->SaveAs(Form("hgrid_ev%d.root",event));
  hbox->SaveAs(Form("hbox_ev%d.root",event));

  return make_tuple(inx,iny,inz,outx,outy,outz);
  
}

void DrawEvent(int event, vector<float> xcrs, vector<float> ycrs, vector<float> zcrs, vector<float> rcrs, vector<float> xhit, vector<float> yhit, vector<float> zhit, vector<float> xlim, vector<float> ylim, vector<float> zlim, vector<float> nkinks, vector<float> xfil, vector<float> yfil, vector<float> zfil){

  TCanvas *c1 = new TCanvas("c1","c1",500,500);

  gSystem->Load("libGeom");
  //delete previous geometry objects in case this script is re-executed
  if (gGeoManager) {
      gGeoManager->GetListOfNodes()->Delete();
      gGeoManager->GetListOfShapes()->Delete();
  }
  
  TRandom3 rnd;
  //if(event==2){
  cout << "Drawing event " << event << endl;
  
  c1->SetFillColor(1);
  TView *view = TView::CreateView(1);
  view->SetRange(-200,-200,-200,200,200,200);
  TBRIK *brik  = new TBRIK("BRIK","BRIK","void",100,100,100);
  
  Double_t x0= xcrs.at(0);
  Double_t y0= ycrs.at(0);
  Double_t z0= zcrs.at(0);
  Double_t pos[3]={1,1,1};
  TNode *node  = new TNode("NODE1","NODE1","BRIK");
  node->cd();
  TNode *node1[rcrs.size()];
  for(int j=0;j<rcrs.size();j++){
    TSPHE *sph = new TSPHE("sph","sphere","air",0,rcrs.at(j),0,180,0,360);
    node1[j]  = new TNode("NODE1","NODE1","sph",xcrs.at(j)-x0, ycrs.at(j)-y0,zcrs.at(j)-z0);
    node1[j]->SetFillStyle(3001);
    node1[j]->SetLineColor(0);
    node->cd();
    node->Draw("gl");
    c1->Update();
  }
	 
  TPolyLine3D *line =new TPolyLine3D(xhit.size());
  for(int i=0;i<xhit.size();i++){
    line->SetPoint(i,xhit.at(i)-x0, yhit.at(i)-y0,zhit.at(i)-z0);
  }
  line->Draw("");
  line->SetLineWidth(2);
  line->SetLineColor(kYellow);
  line->SetLineStyle(2);
  c1->Update();
  
  TPolyMarker3D *pm3d[xlim.size()]; 
  for(int k=0;k<xlim.size();k++){
    pm3d[k] = new TPolyMarker3D(1);
    pm3d[k]->SetName(Form("pm3d%d",k+1));
    pm3d[k]->SetPoint(0,xlim.at(k)-x0, ylim.at(k)-y0,zlim.at(k)-z0);
    pm3d[k]->SetMarkerSize(1.5);
    pm3d[k]->SetMarkerStyle(8);
    pm3d[k]->SetMarkerColor(2);
    pm3d[k]->Draw();
    c1->Update();
  }
	   
  int istep=0;
  int iline=0;
  
  TPolyLine3D *l[nkinks.size()];
  for(int i=0;i<nkinks.size();i++){
    l[i] = new TPolyLine3D(nkinks.at(i)+1);
    for(int j=0;j<nkinks.at(i)+1;j++){
      l[i]->SetPoint(j,xfil.at(istep)-x0, yfil.at(istep)-y0,zfil.at(istep)-z0);
      istep++;
    }
    l[i]->Draw("PL");
    l[i]->SetLineWidth(2);
    int icolor = rnd.Uniform(2,10);
    l[i]->SetLineColor(icolor);
    c1->Update();
  }
  //c1->SaveAs("newtest.C");
}

// grid
void DrawGrid(int event, vector<float> xcrs, vector<float> ycrs, vector<float> zcrs, vector<float> rcrs, vector<float> nkinks, vector<float> xfil, vector<float> yfil, vector<float> zfil, vector<int> inx, vector<int> iny, vector<int> inz, vector<int> outx, vector<int> outy, vector<int> outz){

  TRandom3 rnd;
  //if(event==2){
  cout << "Drawing event " << event << endl;
  TCanvas *c2 = new TCanvas("c2","c2",500,500);
  c2->SetFillColor(1);
  TView *view = TView::CreateView(1);
  //view->SetRange(-200,-200,-200,200,200,200);
  TBRIK *brik  = new TBRIK("BRIK","BRIK","void",500,500,500);
  //const Int_t n = 500;
  //r = new TRandom();
  Double_t x0= xcrs.at(0);
  Double_t y0= ycrs.at(0);
  Double_t z0= zcrs.at(0);

  //Double_t unit_cube = 2; //2nm
  
  Double_t pos[3]={1,1,1};
  TNode *node  = new TNode("NODE1","NODE1","BRIK");
  node->cd();
  for(int j=0;j<rcrs.size();j++){
    TSPHE *sph = new TSPHE("sph","sphere","air",0,rcrs.at(j),0,180,0,360);
    TNode *node1  = new TNode("NODE1","NODE1","sph",xcrs.at(j)-x0, ycrs.at(j)-y0,zcrs.at(j)-z0);
    node1->SetFillStyle(3001);
    //sph->SetFillColor(0);
    node1->SetLineColor(0);
    node->cd();
    node->Draw("gl");
    c2->Update();
    //cout << "node " << xcrs.at(j)-x0 << " " << ycrs.at(j)-y0 << " " << zcrs.at(j)-z0 << endl;
  }

  // INSIDE Ag
  for(int j=0;j<inx.size();j++){
    if(j%1==0){
      //cout << "cubes " << j << " " << inx.at(j) << " " << iny.at(j) << " " << inz.at(j) << endl;
      Double_t cx = -2200 + (unit_cube*inx.at(j));
      Double_t cy = -2200 + (unit_cube*iny.at(j));
      Double_t cz = 0 + (unit_cube*inz.at(j));

      //cout << "cx " << cx << " " << cy << " " << cz << " " << inx.at(j) << " " << iny.at(j) << " " << inz.at(j) << endl;
      
      TBRIK *cube = new TBRIK("cube","cube","void",2,2,2);
      TNode *node1  = new TNode("NODE1","NODE1","cube",-cx+x0,-cy+y0,-cz+z0);
      //node1->SetFillStyle(3001);
      //sph->SetFillColor(0);
      node1->SetLineColor(2);
      node->cd();
      node->Draw("gl");
      c2->Update();
      //cout << "node " << xcrs.at(j)-x0 << " " << ycrs.at(j)-y0 << " " << zcrs.at(j)-z0 << endl;
    }
  }

  // OUTSIDE Ag
  for(int j=0;j<outx.size();j++){
    if(j%1==0){
      //cout << "cubes " << j << " " << outx.at(j) << " " << outy.at(j) << " " << outz.at(j) << endl;
      Double_t cx = -2200 + (unit_cube*outx.at(j));
      Double_t cy = -2200 + (unit_cube*outy.at(j));
      Double_t cz = 0 + (unit_cube*outz.at(j));

      //cout << "cx " << cx << " " << cy << " " << cz << " " << ox.at(j) << " " << oy.at(j) << " " << oz.at(j) << endl;
      
      TBRIK *cube = new TBRIK("cube","cube","void",2,2,2);
      TNode *node1  = new TNode("NODE1","NODE1","cube",-cx+x0,-cy+y0,-cz+z0);
      //node1->SetFillStyle(3001);
      //sph->SetFillColor(0);
      node1->SetLineColor(4);
      node->cd();
      node->Draw("gl");
      c2->Update();
      //cout << "node " << xcrs.at(j)-x0 << " " << ycrs.at(j)-y0 << " " << zcrs.at(j)-z0 << endl;
    }
  }
	   
  int istep=0;
  int iline=0;
  
  TPolyLine3D *l[nkinks.size()];
  for(int i=0;i<nkinks.size();i++){
    l[i] = new TPolyLine3D(nkinks.at(i)+1);
    for(int j=0;j<nkinks.at(i)+1;j++){
      Double_t xf0= xfil.at(istep)-x0;
      Double_t yf0= yfil.at(istep)-y0;
      Double_t zf0= zfil.at(istep)-z0;
      //if(j==0)cout << i << " " << j << " " << istep << " " << xf0 << " " <<  endl;
      //if(j>0)cout << i << " " << j << " " << istep << " " << xf0 << " " << (xfil.at(istep)-xfil.at(istep-1)*1000-xf0)  <<  endl;
      if(j==0)l[i]->SetPoint(j,xf0, yf0,zf0);
      else l[i]->SetPoint(j,(xfil.at(istep)-xfil.at(istep-1))-xf0, (yfil.at(istep)-yfil.at(istep-1))-yf0,(zfil.at(istep)-zfil.at(istep-1))-zf0);
      istep++;
    }
    l[i]->Draw("PL");
    l[i]->SetLineWidth(2);
    //l[i]->SetMarkerStyle(8);
    int icolor = rnd.Uniform(2,10);
    l[i]->SetLineColor(icolor);
    c2->Update();
  }
  c2->SaveAs("gridtest.C");
}


tuple<vector<float>,vector<float>,vector<float>,vector<float>,vector<float>> selected_crystal(vector<float> xhit, vector<float> yhit, vector<float> zhit, vector<float> xcr, vector<float> ycr, vector<float> zcr, vector<float> rcr){
  
  vector<float> xcr_sel;
  vector<float> ycr_sel;
  vector<float> zcr_sel;
  vector<float> rcr_sel;
  vector<float> edep_sel;
  
  float min_x=100000;
  float max_x=-100000;
  float min_y=100000;
  float max_y=-100000;
  float min_z=100000;
  float max_z=-100000;
  
  //---- define a smaller region where is located the track ----//
  for(int i=0;i<xhit.size();i++){
    if(xhit.at(i)<min_x)min_x=xhit.at(i);
    if(xhit.at(i)>max_x)max_x=xhit.at(i);
    if(yhit.at(i)<min_y)min_y=yhit.at(i);
    if(yhit.at(i)>max_y)max_y=yhit.at(i);
    if(zhit.at(i)<min_z)min_z=zhit.at(i);
    if(zhit.at(i)>max_z)max_z=zhit.at(i);
  }
  //cout <<"minmax " << min_x << " " << max_x << " " << min_y << " " << max_y << " " << min_z << " " << max_z << endl;
  //---------------------------------------------//
  
  //---- selection of a limited number of crystals ----// 
  for(int j=0;j<rcr.size();j++){
    if((xcr.at(j)+rcr.at(j))>min_x && (xcr.at(j)-rcr.at(j))<max_x && (ycr.at(j)+rcr.at(j))>min_y && (ycr.at(j)-rcr.at(j))<max_y && (zcr.at(j)+rcr.at(j))>min_z && (zcr.at(j)-rcr.at(j))<max_z){
      //cout << "nfound " << endl;
      xcr_sel.push_back(xcr.at(j));
      ycr_sel.push_back(ycr.at(j));
      zcr_sel.push_back(zcr.at(j));
      rcr_sel.push_back(rcr.at(j));
      edep_sel.push_back(0);
    }	 
  }
  //---------------------------------------//

  return make_tuple(xcr_sel,ycr_sel,zcr_sel,rcr_sel,edep_sel);
}

int InsideSphere(TVector3 p1, TVector3 p2, TVector3 sc, float r){

  int res;
  TVector3 dp1;
  TVector3 dp2;
  
  dp1 = p1-sc;
  dp2 = p2-sc;
  
  float dist1 = dp1.Mag();
  float dist2 = dp2.Mag();
  
  if(dist1<=r && dist2<=r) res=2;
  if(dist1<=r && dist2>r) res=-1;
  if(dist1>r && dist2<=r) res=1;
  if(dist1>r && dist2>r) res=0;

  //if(res==-1)cout <<"dist " <<  dist1 << " " << dist2 << endl;
  
  return res;
}

bool InSegment(TVector3 p1, TVector3 p2, TVector3 s1, TVector3 s2){

  bool inside=false;

  TVector3 dp11;
  TVector3 dp12;
  TVector3 dp21;
  TVector3 dp22;

  TVector3 dp;
  dp = p2-p1;
  
  dp11 = p1-s1;
  dp21 = p2-s1;
  dp12 = p1-s2;
  dp22 = p2-s2;

  float dist = dp.Mag();
  float dist11 = dp11.Mag();
  float dist12 = dp12.Mag();
  float dist21 = dp21.Mag();
  float dist22 = dp22.Mag();
  
  if(dist11<dist && dist12<dist && dist21<dist && dist22< dist) inside=true;
  else inside=false;
  
  return inside;
}

/*
   Calculate the intersection of a ray and 
   The line segment is defined from p1 to p2
   The sphere is of radius r and centered at sc
   There are potentially two points of intersection given by
   p = p1 + mu1 (p2 - p1)
   p = p1 + mu2 (p2 - p1)
   Return FALSE if the ray doesn't intersect the sphere.
*/
tuple<TVector3,TVector3> RaySphere(TVector3 p1, TVector3 p2, TVector3 sc, float r)
{
   double a,b,c;
   double bb4ac;
   TVector3 dp;
   float mu1, mu2;
   TVector3 s1,s2;

   dp.SetXYZ(p2.x() - p1.x(),p2.y() - p1.y(),p2.z() - p1.z());
   
   //dp.SetX(p2.x() - p1.x());
   //dp.SetY(p2.y() - p1.y());
   //dp.SetZ(p2.z() - p1.z());
   
   a = dp.x() * dp.x() + dp.y() * dp.y() + dp.z() * dp.z();
   b = 2 * (dp.x() * (p1.x() - sc.x()) + dp.y() * (p1.y() - sc.y()) + dp.z() * (p1.z() - sc.z()));
   c = sc.x() * sc.x() + sc.y() * sc.y() + sc.z() * sc.z();
   c += p1.x() * p1.x() + p1.y() * p1.y() + p1.z() * p1.z();
   c -= 2 * (sc.x() * p1.x() + sc.y() * p1.y() + sc.z() * p1.z());
   c -= r * r;
   bb4ac = b * b - 4 * a * c;
   if(bb4ac <= 0) {
     mu1 = 0;
     mu2 = 0;
     s1.SetXYZ(0,0,0);
     s2.SetXYZ(0,0,0);
   }
   else{
     //cout << "eccomi " << endl;
     mu1 = (-b + sqrt(bb4ac)) / (2 * a);
     mu2 = (-b - sqrt(bb4ac)) / (2 * a);

     s1 = p1 +  mu1*(p2-p1);
     s2 = p1 +  mu2*(p2-p1);
   }
   
   return make_tuple(s1,s2);
}


float StepInSphere(float dist, TVector3 invec, TVector3 outvec, TVector3 s1, TVector3 s2){

  TVector3 dp1;
  TVector3 dp2;
  
  dp1 = outvec-s1;
  dp2 = outvec-s2;
  
  float dist1 = dp1.Mag();
  float dist2 = dp2.Mag();

  float step_inside;
  TVector3 vec;

  if(dist1< dist && dist2< dist) {
    cout << "WARNING: " << endl;
    step_inside=0;
  }

  if(dist1<=dist){
    vec = s1-invec;
    step_inside = vec.Mag();
  }
  if(dist2<=dist){
    vec = s2-invec;
    step_inside = vec.Mag();
  }

  return step_inside;
}



tuple<vector<int>,vector<float>,vector<float>,vector<float>,vector<float>> crystal_coordinates() {

   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("filament_model.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("%scoord_crystals.txt",dir.Data()));

   Int_t ev;
   Float_t x,y,z,r;
   Int_t nlines = 0;
   vector<int> icr;
   vector<float> xcr;
   vector<float> ycr;
   vector<float> zcr;
   vector<float> rcr;
   
   do{
     in >> ev >> z >> x >> y >> r;
     if (!in.good()) break;
     icr.push_back(ev);
     xcr.push_back(x);
     ycr.push_back(y);
     zcr.push_back(z);
     rcr.push_back(r);
     nlines++;
   }while(!in.eof());
   printf(" found %d points\n",nlines);
   in.close();
   return make_tuple(icr,xcr,ycr,zcr,rcr);
}

vector<float>  crystal_energy_deposit(vector<float> xhit, vector<float> yhit, vector<float> zhit, vector<float> enedep, vector<float> xcr_sel, vector<float> ycr_sel, vector<float> zcr_sel, vector<float> rcr_sel, vector<float> edep_sel){
  
  // track energy deposit evaluation in crystals 
  for(int i=1;i<xhit.size();i++){    // for each step of the track (start from the second hit)
    
    TVector3 p1(xhit.at(i-1),yhit.at(i-1),zhit.at(i-1));   // pre step point
    TVector3 p2(xhit.at(i),yhit.at(i),zhit.at(i));         // post step point
    TVector3 dp;
    dp = p2-p1;
    float dist = dp.Mag();   // step length 
    
    // for each crystal selected (sphere)
    for(int j=0;j<rcr_sel.size();j++){
      
      float step_inside=0; // effective step length in a crystal
      TVector3 s1,s2;      // intersection points segment-sphere
      
      // sphere coordinates
      TVector3 sc(xcr_sel.at(j),ycr_sel.at(j),zcr_sel.at(j));
      
      // evaluate if step marginal points are inside the sphere
      int res = InsideSphere(p1,p2,sc,rcr_sel.at(j));
      
      if(res==2)step_inside=dist; // both are inside 
      else{
	
	//int res2 = InsideSphere(s1,s2,sc,rcr_sel.at(j)); // cross-check
	tie(s1,s2) = RaySphere(p1, p2, sc, rcr_sel.at(j)); // evaluate interesection points
	
	if(res==-1) step_inside = StepInSphere(dist,p1, p2, s1, s2);  // only the first is inside
	if(res==1) step_inside = StepInSphere(dist, p2, p1, s1, s2);  // only the second is inside
	if(res==0){                                                   // nobody is inside
	  if(s1.Mag()!=0 && s2.Mag()!=0){                             // the are intersection points
	    
	    // evaluate if the intersection points belong to the segment
	    bool inside = InSegment(p1,p2,s1,s2);        
	    if(inside){         // if inside evaluate step_inside
	      TVector3 vec;
	      vec = s2-s1;
	      step_inside = vec.Mag();
	    }
	    else {             // intersection points are outside therefore step_inside is 0 
	      step_inside=0;
	    }
	  }
	  if(s1.Mag()==0 && s2.Mag()==0){    // no intersection points are found therefore step_inside is 0
	    step_inside=0;
	  }
	}
      }
      
      // energy deposit in crystal evaluated as:  (dE)*(delta_step_inside)/(delta_step_length)
      float inside_fraction = step_inside/dist;
      edep_sel.at(j) += enedep.at(i-1)*inside_fraction;
    }
  } // end of the track

  return edep_sel;
}

float nlatent_image(float ene_thr, float sig_thr, float edep, TRandom3 rnd){

  float gaus[4];
  float min_dist=1000000;
  float nlatent_im=0;
  //float eff1=0.97;
  //float eff2=0.95;
  //float eff3=0.94;

  gaus[0] = rnd.Landau(ene_thr,sig_thr);
  gaus[1] = rnd.Landau(2*ene_thr,sig_thr);
  gaus[2] = rnd.Landau(3*ene_thr,sig_thr);
  gaus[3] = rnd.Landau(4*ene_thr,sig_thr);

  for(int s=0;s<4;s++){
    if(abs(edep-gaus[s])<min_dist){
      min_dist = abs(edep-gaus[s]);
      nlatent_im=s+1;
    }
  }

  return nlatent_im;
}

tuple<vector<vector<float>>,vector<vector<float>>,vector<vector<float>>,vector<vector<float>>,vector<vector<float>>> kink_weights() {

   TString dir = gSystem->UnixPathName(__FILE__);
   dir.ReplaceAll("filament_model.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("%snkinks_weights.txt",dir.Data()));

   Float_t w1,w2,w3,w4,w5;
   Int_t nlines = 0;
   //Int_t Nlines = 0;
   do{
     in >> w1 >> w2 >> w3 >> w4 >> w5;
     if (!in.good()) break;
     nlines++;
   }while(!in.eof());
   in.close();
   //cout << "nlines " << nlines << endl;
   vector<vector<float>> wk1;
   vector<vector<float>> wk2;
   vector<vector<float>> wk3;
   vector<vector<float>> wk4;
   vector<vector<float>> wk5;
   wk1.resize(nlines);
   wk2.resize(nlines);
   wk3.resize(nlines);
   wk4.resize(nlines);
   wk5.resize(nlines);

   in.open(Form("%snkinks_weights.txt",dir.Data()));
   nlines=0;
   do{
     //cout << "hello" << endl;
     in >> w1 >> w2 >> w3 >> w4 >> w5;
     if (!in.good()) break;
     wk1.at(nlines).push_back(0);
     wk1.at(nlines).push_back(w1);
     wk2.at(nlines).push_back(w1);
     wk2.at(nlines).push_back(w1+w2);
     wk3.at(nlines).push_back(w1+w2);
     wk3.at(nlines).push_back(w1+w2+w3);
     wk4.at(nlines).push_back(w1+w2+w3);
     wk4.at(nlines).push_back(w1+w2+w3+w4);
     wk5.at(nlines).push_back(w1+w2+w3+w4);
     wk5.at(nlines).push_back(w1+w2+w3+w4+w5);
     nlines++;
   }while(!in.eof());
   printf(" found %d kinks\n",nlines);
   in.close();
   //cout << "size " << wk1.size() << endl;
   return make_tuple(wk1,wk2,wk3,wk4,wk5);
}

int num_kinks(float flen, float fmax, int nbins, vector<vector<float>> wk1,vector<vector<float>> wk2,vector<vector<float>> wk3,vector<vector<float>> wk4,vector<vector<float>> wk5, TRandom3 rnd){

  //TRandom3 rnd;
  int ibin = flen/(fmax/nbins);
  float val = rnd.Uniform(0,1);
  int nkinks=0;
  if(val>wk1.at(ibin).at(0)&&val<wk1.at(ibin).at(1))nkinks=1;
  if(val>wk2.at(ibin).at(0)&&val<wk2.at(ibin).at(1))nkinks=2;
  if(val>wk3.at(ibin).at(0)&&val<wk3.at(ibin).at(1))nkinks=3;
  if(val>wk4.at(ibin).at(0)&&val<wk4.at(ibin).at(1))nkinks=4;
  if(val>wk5.at(ibin).at(0)&&val<wk5.at(ibin).at(1))nkinks=5;

  // no weigths avaiable
  if(nkinks==0 && ibin==0) nkinks=1;  
  if(nkinks==0 && ibin>=16) nkinks=rnd.Uniform(4,6);
  return nkinks;
}

TVector3 start_latent_image(TVector3 cr_pos, float cr_radius, TRandom3 rnd){

  TVector3 p0lim;
  float u = rnd.Uniform(0,1);
  float v = rnd.Uniform(0,1);
  float phi = 2*pi*u;
  float theta = acos(2*v-1);
  p0lim.SetX(cr_pos.X() + cr_radius*sin(theta)*cos(phi));
  p0lim.SetY(cr_pos.Y() + cr_radius*sin(theta)*sin(phi));
  p0lim.SetZ(cr_pos.Z() + cr_radius*cos(theta));

  return p0lim;
}

tuple<vector<float>,vector<float>,vector<float>,vector<float>,vector<float>,vector<float>> make_filament(TVector3 sc, TVector3 p0lim, vector<float> klen, TRandom3 rnd){

  TVector3 dir;
  dir = p0lim-sc;
  
  //cout << "sc " << sc(0) << " " << sc(1) << " " << sc(2) << endl;
  //cout << "dir " << dir(0) << " " << dir(1) << " " << dir(2) << endl;
  //cout << "p0lim " << p0lim(0) << " " << p0lim(1) << " " << p0lim(2) << endl;

  int npoints = klen.size()+1;
  vector<float> xf;
  vector<float> yf;
  vector<float> zf;
  vector<float> slen;
  vector<float> sthe;
  vector<float> sphi;
  xf.push_back(p0lim(0));
  yf.push_back(p0lim(1));
  zf.push_back(p0lim(2));
  slen.push_back(0);
  sthe.push_back(0);
  sphi.push_back(0);

  float phi,theta,step,dz,ang,u,v;
  int ntimes=0;
  
  for(int i=1;i<npoints;i++){
    if(i==1){
      xf.push_back(0);
      yf.push_back(0);
      zf.push_back(0);
      do{
	u = rnd.Uniform(0,1);
	v = rnd.Uniform(0,1);
	phi = 2*pi*u;
	theta = acos(2*v-1);
	step = 1000*klen.at(i-1); // nanometer
	dz = step/tan(theta);
	xf.at(i) = xf.at(i-1) + step*cos(phi);
	yf.at(i) = yf.at(i-1) + step*sin(phi);
	zf.at(i) = zf.at(i-1) + dz;
	TVector3 dp(xf.at(i)-xf.at(i-1),yf.at(i)-yf.at(i-1),zf.at(i)-zf.at(i-1));
	ang = dp.Angle(dir);
	ntimes++;
      }while(ang>pi/2.);
    }
    else{
      u = rnd.Uniform(0,1);
      v = rnd.Uniform(0,1);
      phi = 2*pi*u;
      theta = acos(2*v-1);
      step = 1000*klen.at(i-1); // nanometer
      dz = step/tan(theta);
      xf.push_back(xf.at(i-1) + step*cos(phi));
      yf.push_back(yf.at(i-1) + step*sin(phi));
      zf.push_back(zf.at(i-1) + dz);
    }
    slen.push_back(step);
    sthe.push_back(theta);
    sphi.push_back(phi);    
  }
  
  return make_tuple(xf,yf,zf,slen,sthe,sphi);
}

vector<float>  make_len_vector(int dim, float tot_len, TH1F *hist){
  
  float fLen=0;                            // length of single filament
  float diff_fil=1000000;                  // gap between total length and single lengths extracted
  vector<float> vflen(dim);                // vector with single filament lengths 
  if(dim==1)vflen.at(0)=tot_len;           // if there is only one filament (fLen=fTLen)
  else{                                    // if there are more filaments
    for(int s=0;s<500;s++){                // extract 500 times to minimize the gap (diff_fil)
      float sum_flen=0;                    // sum of filament lengths
      vector<float> tmp_len(dim);          // temporary set of single lengths
      // for each latent image
      for(int l=0;l<dim;l++){ 
	float rndlen= hist->GetRandom();
	tmp_len.at(l)=rndlen;
	sum_flen += rndlen;
      }
      float tmp_fil = abs(tot_len-sum_flen); // temporary gap between total length and single lengths extracted
      if(tmp_fil<diff_fil){                  // fill vector with single length if gap is at minimum
	vflen.clear();
	vflen.resize(dim);
	for(int l=0;l<dim;l++){
	  vflen.at(l)=tmp_len.at(l);
	}
	diff_fil=tmp_fil;
      }
    }
  }
  return vflen;
}


// SRIM
//SRIM_data srim;
void SRIM_data::Loop() {

  bool is_drawn=false;
  
  // crystal position
  vector<int> icr;
  vector<float> xcr;
  vector<float> ycr;
  vector<float> zcr;
  vector<float> rcr;

  // kink weights
  vector<vector<float>> wk1;
  vector<vector<float>> wk2;
  vector<vector<float>> wk3;
  vector<vector<float>> wk4;
  vector<vector<float>> wk5;
  //wk1.resize(25);
  // track position in SRIM
  vector<float> xhit;
  vector<float> yhit;
  vector<float> zhit;
  vector<float> enedep;

  vector<int> icrs;
  vector<int> nlim;
  vector<float> xcrs;
  vector<float> ycrs;
  vector<float> zcrs;
  vector<float> rcrs;
  vector<float> ecrs;
  
  vector<float> xcr_sel;
  vector<float> ycr_sel;
  vector<float> zcr_sel;
  vector<float> rcr_sel;
  vector<float> edep_sel;

  vector<int> idlim;
  vector<int> crlim;
  vector<float> xlim;
  vector<float> ylim;
  vector<float> zlim;
  vector<float> flen;
  vector<float> nkinks;

  vector<int> idfil;
  vector<int> crfil;
  vector<float> xfil;
  vector<float> yfil;
  vector<float> zfil;
  vector<float> slen;
  vector<float> thef;
  vector<float> phif;
  
  // crystal position
  tie(icr,xcr,ycr,zcr,rcr) = crystal_coordinates();
  
  int event=0;
  int ncrs=0;
  float dx,dy,dz;
  TRandom3 rnd_l;
  TRandom3 rnd;
  float_t ene_thr=4.5;  // arbitrario
  float_t sig_thr=0.45; // arbitrario

  //---- DATA HISTOGRAMS ----//
  TFile *fin = new TFile("data_plots.root","READ");
  TH1F *hndata = (TH1F*)fin->Get("HfilamentN");
  hndata->Scale(1./hndata->Integral());
  TH1F *hflen = (TH1F*)fin->Get("HfilamentLen");
  hflen->Scale(1./hflen->Integral());
  TH1F *hkink = (TH1F*)fin->Get("kink");
  TH1F *hklen = (TH1F*)fin->Get("slen");
  hklen->Scale(1./hklen->Integral());
  TH1F *htlen = (TH1F*)fin->Get("HtotalLen");
  htlen->Scale(1./htlen->Integral());
  TH1F *hfil1 = (TH1F*)fin->Get("hfil1");
  hfil1->Scale(1./hfil1->Integral());
  TH1F *hfil2 = (TH1F*)fin->Get("hfil2");
  hfil2->Scale(1./hfil2->Integral());

  //---- SIM HISTOGRAMS ----//
  TH1F *he0 = new TH1F("he0","he0",1000,0,100);
  TH1F *he = new TH1F("he","he",1000,0,100);
  TH1F *hn = new TH1F("hn","hn",6,0,6);

  TH1F *hrlen = new TH1F("hrlen","hrlen",25,0,0.5);
  TH1F *htotlen = new TH1F("htotlen","htotlen",20,0,1);
  TH1F *hrkink = new TH1F("hrkink","hrkink",6,0,6);
  TH1F *hdiff = new TH1F("hdiff","hdiff",250,0,0.05);
  TH1F *hslen = new TH1F("hslen","hslen",25,0,0.2);
  TH1F *hsthe = new TH1F("hsthe","hsthe",50,-1,1);
  TH1F *hsphi = new TH1F("hsphi","hsphi",100,0,2*pi);

  //---- OUTPUT TREE ----//
  TFile *fout = new TFile("filament.root","RECREATE");
  TTree *Tree = new TTree("data","Filament");
  Tree->Branch("event",&event);
  Tree->Branch("xhit",&xhit);
  Tree->Branch("yhit",&yhit);
  Tree->Branch("zhit",&zhit);
  Tree->Branch("enedep",&enedep);
  Tree->Branch("ncrs",&ncrs);
  Tree->Branch("icrs",&icrs);
  Tree->Branch("nlim",&nlim);
  Tree->Branch("xcrs",&xcrs);
  Tree->Branch("ycrs",&ycrs);
  Tree->Branch("zcrs",&zcrs);
  Tree->Branch("rcrs",&rcrs);
  Tree->Branch("ecrs",&ecrs);
  Tree->Branch("idlim",&idlim);
  Tree->Branch("crlim",&crlim);
  Tree->Branch("xlim",&xlim);
  Tree->Branch("ylim",&ylim);
  Tree->Branch("zlim",&zlim);
  Tree->Branch("flen",&flen);
  Tree->Branch("nkinks",&nkinks);
  Tree->Branch("idfil",&idfil);
  Tree->Branch("crfil",&crfil);
  Tree->Branch("xfil",&xfil);
  Tree->Branch("yfil",&yfil);
  Tree->Branch("zfil",&zfil);
  Tree->Branch("slen",&slen);

  
 

  // kink weights
  tie(wk1,wk2,wk3,wk4,wk5) = kink_weights();
  
  if (fChain == 0) return;

   // FILE SRIM CON LE TRACCE
   Long64_t nentries = fChain->GetEntriesFast();
   //cout << "nentries " << nentries << endl;
   
   for(Long64_t jentry=0; jentry<nentries;jentry++) {
     
     GetEntry(jentry);

     //----- tracks displaced on crystal pattern -----//
     if(kene==30){                    // to modify is the initial kinetic energy of simulated particle 
       dx = rnd.Uniform(-2000,-1500); // range depends on pattern  
       dy = rnd.Uniform(-2000,-1500);
       dz = rnd.Uniform(250,750);
     }
     xhit.push_back(z+dx);
     yhit.push_back(y+dy);
     zhit.push_back(x+dz);
     enedep.push_back(dene);
     //-----------------------------//

     //---- end of the  track (last hit) ----//
     if(rene==0){ 
       
       int cr_index=0;

       // selection of crystal in a smaller range where is located the track
       tie(xcr_sel,ycr_sel,zcr_sel,rcr_sel,edep_sel) = selected_crystal(xhit,yhit,zhit,xcr,ycr,zcr,rcr);

       // energy deposit evaluation in selected crystals
       edep_sel = crystal_energy_deposit(xhit,yhit,zhit,enedep,xcr_sel,ycr_sel,zcr_sel,rcr_sel,edep_sel);
       
       // for each selected crystal look at the total energy deposit
       for(int j=0;j<rcr_sel.size();j++){

	 // crystal coordinates
	 TVector3 sc(xcr_sel.at(j),ycr_sel.at(j),zcr_sel.at(j));
	 
	 //---- n-latent images ----//	 	 
	 float_t min_val = rnd.Gaus(ene_thr,sig_thr); // set threshold for energy deposit (empiric)
	 float nlatent_im=0;     // number of latent images
	 float sum_flen=0;       // sum of filament lengths per crystal
	 if(edep_sel.at(j)!=0){
	   if(edep_sel.at(j)>ene_thr || (edep_sel.at(j)<ene_thr && edep_sel.at(j)>min_val)){   // if energy deposit is over the cut
	     //cout << jentry << " " << j << " " << edep_sel.at(j) << endl;
	     nlatent_im = nlatent_image(ene_thr,sig_thr,edep_sel.at(j),rnd);  // evaluate number of latent images
	     he->Fill(edep_sel.at(j));
	     hn->Fill(nlatent_im);
	     
	   }
	   // cout << "hello " << edep_sel.at(j)*1000 << endl;
	   he0->Fill(edep_sel.at(j)); // keV
	 }
	 //-------------------------------//

	 //---- if latent images are formed ----//
	 if(nlatent_im>0){

	   // only cristals with latent images are stored
	   xcrs.push_back(xcr_sel.at(j));
	   ycrs.push_back(ycr_sel.at(j));
	   zcrs.push_back(zcr_sel.at(j));
	   rcrs.push_back(rcr_sel.at(j));
	   ecrs.push_back(edep_sel.at(j));
	   icrs.push_back(cr_index);
	   nlim.push_back(nlatent_im);
	   //cout << "lim " << nlatent_im << " " << xcrs.at(cr_index) << " " << ycrs.at(cr_index) << " " << zcrs.at(cr_index) << " " << rcrs.at(cr_index)  << endl;

	   //--- extract total filament 2D length per crystal ----//	   
	   float fTLen=0;          
	   if(nlatent_im==1)fTLen = hfil1->GetRandom();
	   if(nlatent_im>1)fTLen = hfil2->GetRandom();
	   //--------------------------------------------//

	   //---- filament lengths ----//
	   vector<float> vflen =  make_len_vector(nlatent_im,fTLen,hflen);
	   
	   //---- filament production ----//
	   TVector3 p0lim;   // vector with latent image coordinates

	   // loop on filaments
	   for(int k=0;k<nlatent_im;k++){
	     
	     //---- latent image start position ----//
	     rnd.SetSeed(0);
	     p0lim = start_latent_image(sc,rcrs.at(rcrs.size()-1),rnd);
	     
	     float fLen = vflen.at(k); // length of the filament 
	     float fmax = 0.5;         // max length of kinks in filament
	     sum_flen += fLen; 

	     //---- evalutate number of kinks in a single filament ----//
	     int nKinks = num_kinks(fLen,fmax,hflen->GetNbinsX(),wk1,wk2,wk3,wk4,wk5,rnd);
	     hrlen->Fill(fLen);
	     hrkink->Fill(nKinks);
	     //cout << flen << " " << nkinks << endl;
	     //------------------------------------------//
	     
	     //---- latent image info stored ----//
	     idlim.push_back(k+1);
	     crlim.push_back(j);
	     //cout << "lim " << event << " " << p0lim(0) << " " << p0lim(1) << p0lim(2) << endl;
	     xlim.push_back(p0lim(0));
	     ylim.push_back(p0lim(1));
	     zlim.push_back(p0lim(2));
	     flen.push_back(fLen);
	     nkinks.push_back(nKinks);
	     //-------------------------------//

	     //---- step lengths ----//
	     vector<float> klen =  make_len_vector(nKinks,fLen,hklen);
	     
	     vector<float> xf;
	     vector<float> yf;
	     vector<float> zf;
	     vector<float> kf;
	     vector<float> thef;
	     vector<float> phif;
	     
	     // make the filament structure
	     tie(xf,yf,zf,kf,thef,phif) =  make_filament(sc,p0lim, klen, rnd);
	     
	     std::random_shuffle ( kf.begin()+1, kf.end() ); // randomizzo i valori
	     
	     for(int l=0;l<xf.size();l++){
	       idfil.push_back(k+1);
	       crfil.push_back(j);
	       xfil.push_back(xf.at(l));
	       yfil.push_back(yf.at(l));
	       zfil.push_back(zf.at(l));
	       slen.push_back(kf.at(l));
	       if(kf.at(l)!=0){
		 //cout << "kf " << kf.at(l) << endl;
		 hslen->Fill(kf.at(l)/1000.);
		 hsphi->Fill(phif.at(l));
		 hsthe->Fill(cos(thef.at(l)));
	       }
	     }
	     // celle
	   }
	   
	   cr_index++;
	 }
	 
	 if(sum_flen!=0)htotlen->Fill(sum_flen);
       }
       /*
       if(!is_drawn){
	 vector<int> inx;
	 vector<int> iny;
	 vector<int> inz;
	 vector<int> outx;
	 vector<int> outy;
	 vector<int> outz;
	 //if(rcrs.size()>0){
	 //DrawEvent(event,xcrs,ycrs,zcrs,rcrs,xhit,yhit,zhit,xlim,ylim,zlim,nkinks,xfil,yfil,zfil);
	 //is_drawn=true;
	 //}
	 
	 //cout <<"srim ev " << ev << endl;
	 //tie(inx,iny,inz,outx,outy,outz) = MakeGrid(event,unit_cube,xcrs,ycrs,zcrs,rcrs,xfil,yfil,zfil,slen);
	 //DrawGrid(event,xcrs,ycrs,zcrs,rcrs,nkinks,xfil,yfil,zfil,inx,iny,inz,outx,outy,outz);
       }
       */

       ncrs = rcrs.size();
       
       Tree->Fill();
       
       xhit.clear();
       yhit.clear();
       zhit.clear();
       enedep.clear();
       xcrs.clear();
       ycrs.clear();
       zcrs.clear();
       rcrs.clear();
       icrs.clear();
       nlim.clear();
       ecrs.clear();       
       idlim.clear();
       crlim.clear();
       xlim.clear();
       ylim.clear();
       zlim.clear();
       flen.clear();
       nkinks.clear();
       idfil.clear();
       crfil.clear();
       xfil.clear();
       yfil.clear();
       zfil.clear();
       slen.clear();
       thef.clear();
       phif.clear();
       xcr_sel.clear();
       ycr_sel.clear();
       zcr_sel.clear();
       rcr_sel.clear();
       edep_sel.clear();
       cout << "event " << event << endl;
       event++;
       if(event>10000)break;
     }
   }

      
   Tree->Write();
   he0->Write();
   he->Write();
   hn->Scale(1./hn->Integral());
   hn->Write();
   hndata->Write();
   hrlen->Scale(1./hrlen->Integral());
   hrlen->Write();
   htotlen->Scale(1./htotlen->Integral());
   htotlen->Write();
   hkink->Scale(1./hkink->Integral());
   hkink->Write();
   hrkink->Scale(1./hrkink->Integral());
   hrkink->Write();
   hdiff->Write();
   hslen->Scale(1./hslen->Integral());
   hslen->Write();
   hsthe->Write();
   hsphi->Write();
   hklen->Write();
   hflen->Write();
   htlen->Write();   
   fout->Close();

   TCanvas *cc = new TCanvas("cc","cc",1500,1000);
   cc->Divide(3,2);
   cc->cd(1);
   hn->Draw("hist");
   hn->GetYaxis()->SetRangeUser(0,hn->GetMaximum()*1.2);
   hn->GetXaxis()->SetTitle("number of latent images per crystal");
   hndata->Draw("sames");
   hn->SetLineColor(kBlue+1);
   hndata->SetLineColor(kRed);
   hn->SetLineWidth(2);
   hndata->SetLineWidth(2);
   hn->SetTitle("MC");
   hndata->SetTitle("data");
   cc->cd(1)->BuildLegend();

   cc->cd(2);
   htotlen->Draw("hist");
   htotlen->GetYaxis()->SetRangeUser(0,htotlen->GetMaximum()*1.2);
   htotlen->GetXaxis()->SetTitle("sum of filament length per crystal [#mum]");
   htlen->Draw("sames");
   htotlen->SetLineColor(kBlue+1);
   htlen->SetLineColor(kRed);
   htotlen->SetLineWidth(2);
   htlen->SetLineWidth(2);
   htotlen->SetTitle("MC");
   htlen->SetTitle("data");
   cc->cd(2)->BuildLegend();

   cc->cd(3);
   hrlen->Draw("hist");
   hrlen->GetYaxis()->SetRangeUser(0,hrlen->GetMaximum()*1.2);
   hrlen->GetXaxis()->SetTitle("filament length [#mum]");
   hflen->Draw("sames");
   hrlen->SetLineColor(kBlue+1);
   hflen->SetLineColor(kRed);
   hrlen->SetLineWidth(2);
   hflen->SetLineWidth(2);
   hrlen->SetTitle("MC");
   hflen->SetTitle("data");
   cc->cd(3)->BuildLegend();
   
   cc->cd(4);
   hrkink->Draw("hist");
   hrkink->GetYaxis()->SetRangeUser(0,hrkink->GetMaximum()*1.2);
   hrkink->GetXaxis()->SetTitle("number of kinks");
   hkink->Draw("sames");
   hrkink->SetLineColor(kBlue+1);
   hkink->SetLineColor(kRed);
   hrkink->SetLineWidth(2);
   hkink->SetLineWidth(2);
   hrkink->SetTitle("MC");
   hkink->SetTitle("data");
   cc->cd(4)->BuildLegend();

   cc->cd(5);
   hslen->Draw("hist");
   hslen->GetYaxis()->SetRangeUser(0,hslen->GetMaximum()*1.2);
   hslen->GetXaxis()->SetTitle("step length [#mum]");
   hklen->Draw("sames");
   hslen->SetLineColor(kBlue+1);
   hklen->SetLineColor(kRed);
   hslen->SetLineWidth(2);
   hklen->SetLineWidth(2);
   hslen->SetTitle("MC");
   hklen->SetTitle("data");
   cc->cd(5)->BuildLegend();
}

int arun(){
  // gSystem->Load(".L filament_model.C");
  SRIM_data srim;
  srim.Loop();
  return 0;
}
