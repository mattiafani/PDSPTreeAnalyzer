#include "MCCorrection.h"

MCCorrection::MCCorrection() :
  IgnoreNoHist(false)
{

  histDir = TDirectoryHelper::GetTempDirectory("MCCorrection");

}

void MCCorrection::ReadHistograms(){

  TString datapath = getenv("DATA_DIR");

  TDirectory* origDir = gDirectory;

  //==== Momentum reweight
  TString MomentumReweight_path = datapath + "/Momentum_reweight/";
  TString MomentumReweight_path_xrootd = "/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/PDSPTreeAnalyzer/data/v1/Momentum_reweight/";
  TString is_dunegpvm_str = getenv("PDSPAna_isdunegpvm");
  if(is_dunegpvm_str == "TRUE"){
    MomentumReweight_path_xrootd = TString(getenv("PDSPAna_WD")) + "/data/v1/Momentum_reweight/";
  }
  
  cout << "[MCCorrection::ReadHistograms] MomentumReweight_path  : " << MomentumReweight_path << endl;
  string elline;
  cout << "[MCCorrection::ReadHistograms] Open : " << MomentumReweight_path + "histmap.txt" << endl;
  ifstream in(MomentumReweight_path + "histmap.txt");
  while(getline(in,elline)){
    std::istringstream is( elline );

    TString tstring_elline = elline;
    if(tstring_elline.Contains("#")) continue;
    cout << "[MCCorrection::ReadHistograms] tstring_elline : " << tstring_elline << endl;
    TString a,b,c,d,e,f;
    is >> a; // Beam,
    is >> b; // SF, Eff
    is >> c; // <Flag>
    is >> d; // <rootfilename>
    is >> e; // <histname>
    is >> f; // Class
    TFile *file = TFile::Open(MomentumReweight_path_xrootd + d);
    //cout << "MomentumReweight_path + d : " << MomentumReweight_path + d << endl;
    cout << "a : " << a << ", b : " << b << ", c : " << c << ", d : " << d << ", e : " << e << ", f : " << f << endl;

    if(f=="TH1D"){
      cout << "[MCCorrection::ReadHistograms] Saving to map_hist_MomentumReweight[" + a + "_" + b + "_" + c + "]" << endl; 
      histDir->cd();
      map_hist_MomentumReweight[a + "_" + b + "_" + c] = (TH1D *)file -> Get(e) -> Clone();
      cout << "[MCCorrection::ReadHistograms] Saved to map_hist_MomentumReweight map" << endl;
    }
  }

  cout << "[MCCorrection::ReadHistograms] map_hist_MomentumReweight :" << endl;
  for(std::map< TString, TH1D* >::iterator it = map_hist_MomentumReweight.begin(); it != map_hist_MomentumReweight.end(); it++){
    cout << "[MCCorrection::ReadHistograms] key = " << it -> first << endl;
  }

  TString SCE_map_path_xrootd = "/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/PDSPTreeAnalyzer/data/v1/SCE/SCE_DataDriven_180kV_v4.root";
  if(is_dunegpvm_str == "TRUE"){
    SCE_map_path_xrootd = TString(getenv("PDSPAna_WD")) + "/data/v1/SCE/SCE_DataDriven_180kV_v4.root";
  }

  cout << "[MCCorrection::ReadHistograms] SCE_map_path_xrootd : " << SCE_map_path_xrootd << endl;
  TFile *sce_file = TFile::Open(SCE_map_path_xrootd);
  xneg = (TH3F*)sce_file->Get("Reco_ElecField_X_Neg");
  yneg = (TH3F*)sce_file->Get("Reco_ElecField_Y_Neg");
  zneg = (TH3F*)sce_file->Get("Reco_ElecField_Z_Neg");
  xpos = (TH3F*)sce_file->Get("Reco_ElecField_X_Pos");
  ypos = (TH3F*)sce_file->Get("Reco_ElecField_Y_Pos");
  zpos = (TH3F*)sce_file->Get("Reco_ElecField_Z_Pos");

  pos_hists = {xpos, ypos, zpos};
  neg_hists = {xneg, yneg, zneg};

  cout << "[MCCorrection::ReadHistograms] END" << endl;
}

MCCorrection::~MCCorrection(){

}

void MCCorrection::SetIsData(bool b){
  IsData = b;
}

double MCCorrection::MomentumReweight_SF(TString flag, double P_beam_inst, int sys){

  double value = 1.;
  double error = 0.;

  TH1D *this_hist = map_hist_MomentumReweight["Beam_SF_" + flag];
  if(!this_hist){
    if(IgnoreNoHist) return 1.;
    else{
      cerr << "[MCCorrection::MuonReco_SF] No "<< "Beam_SF_" + flag << endl;
      exit(EXIT_FAILURE);
    }
  }

  int this_bin = this_hist->FindBin(P_beam_inst);
  value = this_hist->GetBinContent(this_bin);
  error = this_hist->GetBinError(this_bin);

  return value+double(sys)*error;

}

double MCCorrection::dEdx_scaled(double MC_dEdx, double slope_fraction){
  double f_const_p0 = 0.997;
  double f_pol1_p0 = 0.945;
  double f_pol1_p1 = 0.021; // maximum slope
  double this_slope = f_pol1_p1 * slope_fraction;

  double scale = 1.;
  if(MC_dEdx < 2.44) scale = f_const_p0;
  //else scale = f_pol1_p0 + MC_dEdx * f_pol1_p1;
  else scale = this_slope * (MC_dEdx - 2.44) + f_const_p0; 

  //if(scale > 1.2) scale = 1.2;

  //cout << "[dEdx_res::dEdx_scaled] MC_dEdx : " << MC_dEdx << ", scale : " << scale << ", slope_fraction : " << slope_fraction << endl; 
  //return scale * MC_dEdx; // This is correct one


  double this_dEdx = MC_dEdx;
  if(MC_dEdx > 1.8 && MC_dEdx < 2.3) this_dEdx = MC_dEdx * 1.04;
  return this_dEdx; // FIXME : just for test
}

double MCCorrection::Use_Other_Mod_Box_Params(double dEdx, double Efield, double alpha_new, double beta_new, double calib_const_ratio){
  double alpha_default = 0.93;
  double beta_default = 0.212;
  double rho = 1.40;

  double exp_term = exp( calib_const_ratio * (beta_new / beta_default) * log(alpha_default + beta_default * dEdx / (rho * Efield)) );
  double new_dEdx = (exp_term - alpha_new) * rho * Efield / beta_new;

  return new_dEdx;
}

double MCCorrection::Use_Abbey_Recom_Params(double dEdx, double Efield, double calib_const_ratio){

  double alpha_default = 0.93;
  double beta_default = 0.212;
  double rho = 1.40;

  // == Change these values
  double alpha_Abbey = 0.91;
  double beta_Abbey = 0.214;
  calib_const_ratio = 0.958;

  double exp_term = exp( calib_const_ratio * (beta_Abbey / beta_default) * log(alpha_default + beta_default * dEdx / (rho * Efield)) );
  double new_dEdx = (exp_term - alpha_Abbey) * rho * Efield / beta_Abbey;
  
  //if(fabs(dEdx - new_dEdx) > 0.00001) cout << "[MCCorrection::Use_Abbey_Recom_Params] dEdx : " << dEdx << ", new_dEdx : " << new_dEdx << endl;

  return new_dEdx;
}

double MCCorrection::SCE_Corrected_E(double x, double y, double z){
  double E0value=0.4867;
  if (sce=="on"){
    if(x>=0){
      for (auto h : pos_hists) {
        if (h->GetXaxis()->FindBin(x) < 1 ||
            h->GetXaxis()->FindBin(x) > h->GetNbinsX()) {
	  std::cout << "x oob: " << x << std::endl;
        }
        if (h->GetYaxis()->FindBin(y) < 1 ||
            h->GetYaxis()->FindBin(y) > h->GetNbinsY()) {
	  std::cout << "y oob: " << y << std::endl;
        }
        if (h->GetZaxis()->FindBin(z) < 1 ||
            h->GetZaxis()->FindBin(z) > h->GetNbinsZ()) {
	  std::cout << "z oob: " << z << std::endl;
        }
      }
      float ex=E0value+E0value*xpos->Interpolate(x,y,z);
      float ey=E0value*ypos->Interpolate(x,y,z);
      float ez=E0value*zpos->Interpolate(x,y,z);
      return sqrt(ex*ex+ey*ey+ez*ez);
    }
    if(x<0){
      for (auto h : neg_hists) {
        if (h->GetXaxis()->FindBin(x) < 1 ||
            h->GetXaxis()->FindBin(x) > h->GetNbinsX()) {
	  std::cout << "x oob: " << x << std::endl;
        }
        if (h->GetYaxis()->FindBin(y) < 1 ||
            h->GetYaxis()->FindBin(y) > h->GetNbinsY()) {
	  std::cout << "y oob: " << y << std::endl;
	}
        if (h->GetZaxis()->FindBin(z) < 1 ||
            h->GetZaxis()->FindBin(z) > h->GetNbinsZ()) {
	  std::cout << "z oob: " << z << std::endl;
	}
      }
      float ex=E0value+E0value*xneg->Interpolate(x,y,z);
      float ey=E0value*yneg->Interpolate(x,y,z);
      float ez=E0value*zneg->Interpolate(x,y,z);
      return sqrt(ex*ex+ey*ey+ez*ez);
    }
  }
  return E0value;
}
