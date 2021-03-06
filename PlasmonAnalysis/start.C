#include "TSystem.h"
#include "TROOT.h"
#include "Riostream.h"
#include <utility>
#include <vector>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <stdio.h>
#include "TString.h"
#include "TUnixSystem.h"
#include "TTree.h"
#include <map>


using namespace std;

void start(){

  char* pPath;
  pPath = getenv ("DM_ROOT");
  gROOT->ProcessLine("gSystem->Load(\"libDMRoot\")");
  gROOT->ProcessLine("Vdmr->MakeClass(\"myData\")");
  gSystem->Exec("rm myData.C");
  gSystem->Exec("cp $DM_ROOT/src/macros/myData_v5.C .");
  gSystem->Exec("cp $DM_ROOT/src/macros/settings.mac .");
  //gSystem->Exec("cp $DM_ROOT/src/macros/run_gif.C .");
  gSystem->Exec("cp $DM_ROOT/src/macros/node.sh .");
  
  TString dir = gSystem->UnixPathName(__FILE__);
  dir.ReplaceAll("start.C","");
  dir.ReplaceAll("/./","/");
  ifstream in;
  in.open(Form("%smyData.h",dir.Data()));
  
  ofstream out("MyData.h");
  
  int iline=0;
  vector<vector<string> > runDescription = vector<vector<string> >();
      if (in.is_open()){
	string line;
	out << line << endl;
	while (getline (in,line)){
	  iline++;
	  istringstream iss(line);
	  string a;
	  vector<string> tmpS;
	  if((iline<16 || iline>25) && iline!=384){
	    while(!iss.eof()){     
	      iss>>a;
	      tmpS.push_back(a);
	      out << a << " ";
	      
	    }
	    runDescription.push_back(tmpS);
	    tmpS.clear();
	    out << endl;
	  }
	  if(iline==16){
	    out <<"#include <TObject.h>" << endl;
	    out <<"#include <TClonesArray.h>" << endl;
	    out <<"#include \""<< pPath << "/include/DMRViewHeader.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRAffine2D.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRCluster.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRGrain.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRMicrotrack.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRImage.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRImage.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRImage.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRRun.h\"" << endl;
	    out <<"#include \""<< pPath << "/include/DMRRunHeader.h\"" << endl;
	  }
	  if(iline==384) out << "myData::myData(TTree *tree) : fChain(0)" << endl;
	}
	out.close();
	in.close();
      }
      gSystem->Exec("rm myData.h");
      gSystem->Exec("mv MyData.h myData.h");
      gROOT->ProcessLine(".q");
}
