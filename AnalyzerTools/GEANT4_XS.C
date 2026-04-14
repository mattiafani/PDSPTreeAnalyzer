#include "GEANT4_XS.h"

GEANT4_XS::GEANT4_XS() :
  IgnoreNoHist(false)
{

  histDir = TDirectoryHelper::GetTempDirectory("GEANT4_XS");

}

void GEANT4_XS::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");

  TDirectory* origDir = gDirectory;

  //==== Xsec
  TString Xsec_path = datapath+"/GEANT4_XS/";
  TString Xsec_path_xrootd = "/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/PDSPTreeAnalyzer/data/v1/GEANT4_XS/";
  TString is_dunegpvm_str = getenv("PDSPAna_isdunegpvm");
  if(is_dunegpvm_str == "TRUE"){
    Xsec_path_xrootd = TString(getenv("PDSPAna_WD")) + "/data/v1/GEANT4_XS/";
  }

  cout << "[GEANT4_XS::ReadHistograms] Open : " << Xsec_path+"histmap.txt" << endl;
  //cout << "[GEANT4_XS::ReadHistograms] Xsec_path  : " << Xsec_path << endl;
  string elline;
  ifstream in(Xsec_path+"histmap.txt");
  while(getline(in,elline)){
    std::istringstream is( elline );

    TString tstring_elline = elline;
    if(tstring_elline.Contains("#")) continue;
    //cout << "[GEANT4_XS::ReadHistograms] tstring_elline : " << tstring_elline << endl;
    TString a,b,c,d,e,f,g;
    is >> a; // Xsec,
    is >> b; // KE, P
    is >> c; // pion, proton
    is >> d; // <process>
    is >> e; // <filename>
    is >> f; // <histname>
    is >> g; // Class
    TFile *file = TFile::Open(Xsec_path_xrootd + e);
    //cout << "Xsec_path + d : " << Xsec_path + d << endl;
    //cout << "a : " << a << ", b : " << b << ", c : << " << c << ", d : " << d << ", e : " << e << ", f : " << f << ", g : " << g << endl;

    if(g=="TGraph"){
      histDir->cd();
      map_graph_Xsec[b + "_" + c + "_" + d] = (TGraph*)file -> Get(f) -> Clone();
    }
  }

  cout << "[GEANT4_XS::ReadHistograms] map_graph_Xsec : " << endl;
  for(std::map< TString, TGraph* >::iterator it = map_graph_Xsec.begin(); it != map_graph_Xsec.end(); it++){
    cout << "[GEANT4_XS::ReadHistograms] key = " << it -> first << endl;
  }

  cout << "[GEANT4_XS::ReadHistograms] END" << endl;
}

GEANT4_XS::~GEANT4_XS(){

}

void GEANT4_XS::SetIsData(bool b){
  IsData = b;
}
