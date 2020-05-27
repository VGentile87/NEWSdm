// crystal_coordinates function reads a file with the coordinates and radius of AgBr crystals
// The function returns a set of vectors with the above mentioned info

void make_datafilament() {
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("make_datafilament.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%sC30V_segmentedfilaments.txt",dir.Data()));
  
  Int_t fid, kid, xpix, ypix;
  Float_t slen;
  Int_t nlines = 0;

  Float_t flen=0;
  Int_t tmp_kid=0;
  
  TFile *fout = new TFile("C30keV_filament.root","RECREATE");
  TTree *tree = new TTree("tree","data");
  tree->Branch("fid",&fid);
  tree->Branch("kid",&kid);
  tree->Branch("xpix",&xpix);
  tree->Branch("ypix",&ypix);
  tree->Branch("slen",&slen);

  TH1F * hflen = new TH1F("flen","flen",25,0,0.5);
  TH1F * hflen1 = new TH1F("flen1","flen1",25,0,0.5);
  TH1F * hflen2 = new TH1F("flen2","flen2",25,0,0.5);
  TH1F * hflen3 = new TH1F("flen3","flen3",25,0,0.5);
  TH1F * hflen4 = new TH1F("flen4","flen4",25,0,0.5);
  TH1F * hflen5 = new TH1F("flen5","flen5",25,0,0.5);
  TH1F * hslen = new TH1F("slen","slen",20,0,0.2);
  TH1F * hslen1 = new TH1F("slen1","slen1",20,0,0.2);
  TH1F * hslen2 = new TH1F("slen2","slen2",20,0,0.2);
  TH1F * hslen3 = new TH1F("slen3","slen3",20,0,0.2);
  TH1F * hslen4 = new TH1F("slen4","slen4",20,0,0.2);
  TH1F * hslen5 = new TH1F("slen5","slen5",20,0,0.2);

  TH1F * hkink = new TH1F("kink","kink",6,0,6);

  vector<vector<int>> vfid;
  vfid.resize(200);
  vector<vector<float>> vlen;
  vlen.resize(200);

  do{
    in >> fid >> kid >> xpix >> ypix >> slen;
    if (!in.good()) break;
    if(kid!=0){
      vfid[fid].push_back(kid);
      vlen[fid].push_back(slen);
      hslen->Fill(slen);
      flen +=slen;
    }
    tree->Fill();
    if(fid!=0 && kid==0){
      hflen->Fill(flen);
      if(tmp_kid==1)hflen1->Fill(flen);
      if(tmp_kid==2)hflen2->Fill(flen);
      if(tmp_kid==3)hflen3->Fill(flen);
      if(tmp_kid==4)hflen4->Fill(flen);
      if(tmp_kid==5)hflen5->Fill(flen);
      flen=0;
    }
    tmp_kid=kid;
    if(in.eof())hflen->Fill(flen);
    nlines++;
  }while(!in.eof());
  printf(" found %d points\n",nlines);


  for(int i=0;i<vfid.size();i++){
    for(int j=0;j<vfid.at(i).size();j++){
      if(vfid.at(i).size()==1)hslen1->Fill(vlen.at(i).at(j));
      if(vfid.at(i).size()==2)hslen2->Fill(vlen.at(i).at(j));
      if(vfid.at(i).size()==3)hslen3->Fill(vlen.at(i).at(j));
      if(vfid.at(i).size()==4)hslen4->Fill(vlen.at(i).at(j));
      if(vfid.at(i).size()==5)hslen5->Fill(vlen.at(i).at(j));
    }
  }

  hkink->SetBinContent(2,hflen1->GetEntries());
  hkink->SetBinContent(3,hflen2->GetEntries());
  hkink->SetBinContent(4,hflen3->GetEntries());
  hkink->SetBinContent(5,hflen4->GetEntries());
  hkink->SetBinContent(6,hflen5->GetEntries());
  
  in.close();
  tree->Write();
  hflen->Write();
  hflen1->Write();
  hflen2->Write();
  hflen3->Write();
  hflen4->Write();
  hflen5->Write();
  hslen->Write();
  hslen1->Write();
  hslen2->Write();
  hslen3->Write();
  hslen4->Write();
  hslen5->Write();
  hkink->Write();

  hflen1->Divide(hflen);
  hflen2->Divide(hflen);
  hflen3->Divide(hflen);
  hflen4->Divide(hflen);
  hflen5->Divide(hflen);

  ofstream log_w("nkinks_weigth.dat");
  for(int i=0;i<hflen->GetNbinsX();i++){
    log_w << hflen1->GetBinContent(i+1) << " " << hflen2->GetBinContent(i+1) << " " << hflen3->GetBinContent(i+1) << " " << hflen4->GetBinContent(i+1) << " " << hflen5->GetBinContent(i+1) << endl;
  }
  log_w.close();
  
  fout->Close();
  
  //return make_tuple(icr,xcr,ycr,zcr,rcr);
}
