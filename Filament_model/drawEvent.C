// V.Gentile 05/2020
//
void drawEvent(){

   //   Connect file generated in $ROOTSYS/test
   TFile fileIn("filament.root");
   //TTree* theTree = nullptr;    
   TTreeReader theReader("data",&fileIn);
   //TTreeReaderValue<Event> eventRV(theReader, "event");
   TTreeReaderValue<Int_t> event(theReader, "event");
   TTreeReaderArray<Float_t> xcrs(theReader, "xcrs");
   TTreeReaderArray<Float_t> ycrs(theReader, "ycrs");
   TTreeReaderArray<Float_t> zcrs(theReader, "zcrs");
   TTreeReaderArray<Float_t> rcrs(theReader, "rcrs");
   TTreeReaderArray<Float_t> xhit(theReader, "xhit");
   TTreeReaderArray<Float_t> yhit(theReader, "yhit");
   TTreeReaderArray<Float_t> zhit(theReader, "zhit");
   TTreeReaderArray<Float_t> xlim(theReader, "xlim");
   TTreeReaderArray<Float_t> ylim(theReader, "ylim");
   TTreeReaderArray<Float_t> zlim(theReader, "zlim");
   TTreeReaderArray<Float_t> xfil(theReader, "xfil");
   TTreeReaderArray<Float_t> yfil(theReader, "yfil");
   TTreeReaderArray<Float_t> zfil(theReader, "zfil");
   TTreeReaderArray<Float_t> nkinks(theReader, "nkinks");

   int an_event;
   cout << "Select event ID ";
   cin >> an_event;
   cout << "Drawing event " << an_event << endl;

   while(theReader.Next()){

     /*if (*nTracksRV < 587) {
         continue; // Check if we don't have too many tracks
	 }
	 auto event = eventRV.Get();      //Read complete accepted event
	 //in memory.
	 hnseg->Fill(event->GetNseg());   //Fill histogram with number of
	 
     */                          //segments.
    
     if(*event==an_event){
       if(rcrs.GetSize()==0){
	 cout << "No crystal sensitised" << endl;
	 return 0;
       }
       if(rcrs.GetSize()>0){
	 cout << "ev " << *event << " ncrs " << rcrs.GetSize() << endl;
	 
	 TCanvas *c1 = new TCanvas("c1","c1",500,500);
	 
	 gSystem->Load("libGeom");
	 //delete previous geometry objects in case this script is re-executed
	 if (gGeoManager) {
	   gGeoManager->GetListOfNodes()->Delete();
	   gGeoManager->GetListOfShapes()->Delete();
	 }
	 
	 TRandom3 rnd;
       c1->SetFillColor(1);
       TView *view = TView::CreateView(1);
       view->SetRange(-200,-200,-200,200,200,200);
       TBRIK *brik  = new TBRIK("BRIK","BRIK","void",100,100,100);
       
       
       Double_t x0= xcrs.At(0);
       Double_t y0= ycrs.At(0);
       Double_t z0= zcrs.At(0);
       TNode *node  = new TNode("NODE1","NODE1","BRIK");
       node->cd();
       TNode *node1[rcrs.GetSize()];
       for(int j=0;j<rcrs.GetSize();j++){
	 TSPHE *sph = new TSPHE("sph","sphere","air",0,rcrs.At(j),0,180,0,360);
	 //TNode *node1  = new TNode("NODE1","NODE1","sph",xcrs.at(j)-x0, ycrs.at(j)-y0,zcrs.at(j)-z0);
	 node1[j]  = new TNode("NODE1","NODE1","sph",xcrs.At(j)-x0, ycrs.At(j)-y0,zcrs.At(j)-z0);
	 node1[j]->SetFillStyle(3001);
	 //sph->SetFillColor(0);
	 node1[j]->SetLineColorAlpha(kWhite, 0.15);
	 //node1[j]->SetLineColor(0);
	 node->cd();
	 node->Draw("gl");
	 //node1[j]->SavePrimitive(os);
	 c1->Update();
	 //cout << "node " << xcrs.at(j)-x0 << " " << ycrs.at(j)-y0 << " " << zcrs.at(j)-z0 << endl;
       }
       
       TPolyLine3D *line =new TPolyLine3D(xhit.GetSize());
       for(int i=0;i<xhit.GetSize();i++){
	 //cout <<"hit node "<<  i << " " <<  xhit.at(i) << " " << yhit.at(i) << " " << zhit.at(i) <<  endl;
	 line->SetPoint(i,xhit.At(i)-x0, yhit.At(i)-y0,zhit.At(i)-z0);
       }
       line->Draw("");
       line->SetLineWidth(2);
       line->SetLineColor(kYellow);
       line->SetLineStyle(2);
       c1->Update();
       
       TPolyMarker3D *pm3d[xlim.GetSize()]; 
       for(int k=0;k<xlim.GetSize();k++){
	 pm3d[k] = new TPolyMarker3D(1);
	 pm3d[k]->SetName(Form("pm3d%d",k+1));
	 pm3d[k]->SetPoint(0,xlim.At(k)-x0, ylim.At(k)-y0,zlim.At(k)-z0);
	 pm3d[k]->SetMarkerSize(1.5);
	 pm3d[k]->SetMarkerStyle(8);
	 pm3d[k]->SetMarkerColor(2);
	 pm3d[k]->Draw();
	 c1->Update();
       }
       
       int istep=0;
       int iline=0;
       
       TPolyLine3D *l[nkinks.GetSize()];
       for(int i=0;i<nkinks.GetSize();i++){
	 l[i] = new TPolyLine3D(nkinks.At(i)+1);
	 for(int j=0;j<nkinks.At(i)+1;j++){
	   l[i]->SetPoint(j,xfil.At(istep)-x0, yfil.At(istep)-y0,zfil.At(istep)-z0);
	   istep++;
	 }
	 l[i]->Draw("PL");
	 l[i]->SetLineWidth(2);
	 //l[i]->SetMarkerStyle(8);
	 int icolor = rnd.Uniform(2,10);
	 l[i]->SetLineColor(icolor);
	 c1->Update();
       }

       c1->GetViewer3D("OPENGL");
       //c1->SaveAs(Form("events/ev%d.C",*event));
       
     }
     }
   }
   // hnseg->Draw();
}
