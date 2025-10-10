#ifndef AnalyzerCore_h
#define AnalyzerCore_h

#include "Math/Vector4D.h"
#include "TString.h"
#include "TMath.h"
#include "TH3.h"
#include "TH2D.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include <sstream>
#include "TRandom.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TFitResult.h"

#include "PDSPTree.h"
#include "../DataFormats/Event.h"
#include "../DataFormats/Daughter.h"
#include "../DataFormats/TrueDaughter.h"
#include "../DataFormats/MCSSegment.h"
#include "../AnalyzerTools/BetheBloch.h"
#include "../AnalyzerTools/GEANT4_XS.h"
#include "../AnalyzerTools/MCCorrection.h"
#include "../AnalyzerTools/Fittings.h"

//#define M_Z 91.1876
//#define M_W 80.379
#define M_mu 105.65837
#define M_neutron 939.565
#define M_proton 938.272
#define M_pion 139.570
#define M_Kaon 493.677
#define M_e 0.510998
#define M_pizero 134.976
#define L_beamline 28.575 // meter
#define v_light 299792458.0 // meter per sec

class AnalyzerCore {

public:

  AnalyzerCore();
  virtual ~AnalyzerCore();

  PDSPTree evt;
  Long64_t MaxEvent, NSkipEvent;
  int LogEvery;
  TString MCSample;
  double Beam_Momentum;
  vector<TString> Userflags;

  bool debug_mode = false;

  virtual void initializeAnalyzer(){

  };

  virtual void executeEvent(){

  };

  //==================
  // Read Tree
  //==================
  virtual void SetTreeName(){
    TString tname = "/pduneana/beamana";
    fChain = new TChain(tname);
  }

  virtual void AddFile(TString filename){
    fChain->Add(filename);
  }

  Int_t GetEntry(Long64_t entry);

  TChain *fChain;

  void Init();
  void Init_evt();

  std::string printcurrunttime(){

    std::stringstream out;
    TDatime datime;
    out << datime.GetYear()<<"-"<<AddZeroToTime(datime.GetMonth())<<"-"<<AddZeroToTime(datime.GetDay())<<" "<<AddZeroToTime(datime.GetHour())<<":"<<AddZeroToTime(datime.GetMinute())<<":"<<AddZeroToTime(datime.GetSecond());
    return out.str();

  }

  std::string AddZeroToTime(int twodigit){
    if(twodigit<10){
      return "0"+std::to_string(twodigit);
    }
    else{
      return std::to_string(twodigit);
    }
  }

  void Loop();

  //==================
  // Get Particles
  //==================
  Event GetEvent();
  std::vector<TrueDaughter> GetAllTrueDaughters();
  std::vector<TrueDaughter> GetPionTrueDaughters(const vector<TrueDaughter>& in);
  std::vector<TrueDaughter> GetProtonTrueDaughters(const vector<TrueDaughter>& in);
  std::vector<TrueDaughter> GetNeutronTrueDaughters(const vector<TrueDaughter>& in);
  std::vector<Daughter> GetAllDaughters();
  std::vector<Daughter> GetDaughters(const vector<Daughter>& in, int cut_Nhit = 10,double cut_beam_dist = 10., double cut_trackScore = 0.5);
  std::vector<Daughter> GetPions(const vector<Daughter>& in);
  std::vector<Daughter> GetStoppingPions(const vector<Daughter>& in, double chi2_pion_cut = 6.);
  std::vector<Daughter> GetProtons(const vector<Daughter>& in);
  std::vector<Daughter> GetTruePions(const vector<Daughter>& in);
  std::vector<Daughter> GetTrueProtons(const vector<Daughter>& in);
  std::vector<Daughter> GetTrueNeutrons(const vector<Daughter>& in);
  
  //==================
  // AnalyzerTools
  //==================
  std::map< int, BetheBloch* > map_BB;
  MCCorrection *MCCorr;
  GEANT4_XS *G4Xsec;
  Fittings *Fitter;
  void initializeAnalyzerTools();

  //==================
  // Set Beam Variables
  //==================
  double beam_TOF = -999.;
  double P_beam_TOF = -999.;
  double P_beam_inst = -999.;
  double KE_beam_inst = -999.;
  double exp_trk_len_beam_inst = -999.;
  double trk_len_ratio = -999.;
  double mass_beam = 139.57;
  double P_ff_reco = -999.;
  double KE_ff_reco = -999.;
  double KE_end_reco = -999.;
  double E_end_reco = -999.;
  double Get_true_ffKE();
  double Get_true_beamlen();
  double P_ff_true = -999.;
  double KE_ff_true = -999.;
  //double KE_ff_true = -999.;

  int pandora_slice_pdg;
  void SetPandoraSlicePDG(int pdg);
  double daughter_michel_score;
  double beam_costh;
  double beam_TPC_theta;
  double beam_TPC_phi;
  double delta_X_spec_TPC;
  double delta_Y_spec_TPC;
  double cos_delta_spec_TPC;
  double chi2_proton;
  double chi2_pion;
  double chi2_muon;
  std::map< int, TProfile* > map_profile;
  int GetPiParType();
  int GetPi2ParType();
  int pi_type = 0;
  TString pi_type_str = "";
  int GetPParType();
  int p_type = 0;
  TString p_type_str = "";
  int pi_truetype = 0;
  int GetPiTrueType();
  
  //==================
  // Event Selections
  //==================
  bool IsData = false;
  bool Pass_Beam_PID(int PID);
  bool PassBeamScraperCut() const;
  double P_beam_inst_scale = 1.;
  double beam_momentum_low = 0.;
  double beam_momentum_high = 10000.;
  bool PassBeamMomentumWindowCut() const;
  bool PassPandoraSliceCut() const;
  bool PassCaloSizeCut() const;  
  bool PassAPA3Cut(double cut = 215.);
  bool PassMichelScoreCut(double cut = 0.55);
  const double beam_angleX_data = 100.464;
  const double beam_angleY_data = 103.442;
  const double beam_angleZ_data = 17.6633;
  const double beam_angleX_mc = 101.579;
  const double beam_angleY_mc = 101.212;
  const double beam_angleZ_mc = 16.5822;
  bool PassBeamCosCut(double cut = 0.95);
  bool PassBeamStartZCut(double cut = -5.);
  bool PassProtonVetoCut(double cut = 80.);
  bool PassMuonVetoCut(double cut = 8.);
  bool PassStoppedPionVetoCut(double cut = 13.);

  //==================
  // Additional Functions
  //==================
  double TOF_to_P(double TOF, int PID);
  double Convert_P_Spectrometer_to_P_ff(double P_beam_inst, TString particle, TString key, int syst);
  double Particle_chi2(const vector<double> & dEdx, const vector<double> & ResRange, int PID, double dEdx_res_frac = 1.);
  double Particle_chi2_skip(const vector<double> & dEdx, const vector<double> & ResRange, int PID, double dEdx_res_frac);
  double max_additional_res_length_pion = 450.; // == [cm] 
  double max_additional_res_length_proton = 120.; // == [cm]
  double res_length_step_pion = 1.0; // == [cm]
  double res_length_step_proton = 0.2; // == [cm]
  double dEdx_truncate_upper_pion = 5.;
  double dEdx_truncate_bellow_pion = 0.5;
  double dEdx_truncate_upper_proton = 20.;
  double dEdx_truncate_bellow_proton = 0.2;
  double Fit_HypTrkLength_Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam);
  double Fit_HypTrkLength_Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam);
  double KE_Hypfit_Gaussian(const Daughter in, int PID);
  double KE_Hypfit_Likelihood(const Daughter in, int PID);
  double Get_EQE_NC_Pion(double P_pion, double cos_theta, double E_binding, int which_sol);
  double Get_EQE_NC_Proton(double P_proton, double cos_theta, double E_binding, int which_sol);
  bool Is_EQE(double window);
  
  //==================
  // MCS
  //==================
  vector<MCSSegment> SplitIntoSegments(const vector<TVector3> & hits, double segment_size);
  TVector3 RotateToZaxis(TVector3 reference_vec, TVector3 input);
  double HL_kappa_a = 0.022;
  double HL_kappa_c = 9.078;
  double HL_sigma_res = 0.001908;
  double HL_epsilon = 0.038;
  double MCS_Get_HL_Sigma(double segment_size, double P, double mass);
  double MCS_Get_Likelihood(double HL_sigma, double delta_angle);
  double MCS_Likelihood_Fitting(const vector<MCSSegment> segments, double segment_size, int PID);

  //==================
  //===Fitting
  //==================
  static Double_t langaufun(Double_t *x, Double_t *par) {
    Double_t invsq2pi = 0.398942280401;

    Double_t np = 500.0;
    Double_t sc = 5.0;
    Double_t xx;
    Double_t mpc;
    Double_t fland;
    Double_t sum = 0.0;
    Double_t xlow,xupp;
    Double_t step;
    Double_t i;

    mpc=par[1];
    xlow = x[0] - sc * par[3];
    xupp = x[0] + sc * par[3];
    step = (xupp-xlow)/np;

    for(i=1.0; i<=np/2; i++) {
      xx = xlow + (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
      xx = xupp - (i-.5) * step;
      fland = TMath::Landau(xx,mpc,par[0]) / par[0];
      sum += fland * TMath::Gaus(x[0],xx,par[3]);
    }

    return (par[2] * step * sum * invsq2pi / par[3]);
  }

  TF1 *langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status, TString FunName);

  //==================
  //===Plotting
  //==================
  std::map< TString, TH1D* > maphist_TH1D;
  std::map< TString, TH2D* > maphist_TH2D;
  std::map< TString, TH3D* > maphist_TH3D;

  TH1D* GetHist1D(TString histname);
  TH2D* GetHist2D(TString histname);
  TH3D* GetHist3D(TString histname);

  void FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  void FillHist(TString histname, double value, double weight, int n_bin, double *xbins);
  void FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double x_min, double x_max,
                int n_biny, double y_min, double y_max);
  void FillHist(TString histname,
                double value_x, double value_y,
                double weight,
                int n_binx, double *xbins,
                int n_biny, double *ybins);
  void FillHist(TString histname,
                double value_x, double value_y, double value_z,
                double weight,
                int n_binx, double x_min, double x_max,
                int n_biny, double y_min, double y_max,
                int n_binz, double z_min, double z_max);
  void FillHist(TString histname,
                double value_x, double value_y, double value_z,
                double weight,
                int n_binx, double *xbins,
                int n_biny, double *ybins,
                int n_binz, double *zbins);

  //==== JSFillHist : 1D
  std::map< TString, std::map<TString, TH1D*> > JSmaphist_TH1D;
  TH1D* JSGetHist1D(TString suffix, TString histname);
  void JSFillHist(TString suffix, TString histname, double value, double weight, int n_bin, double x_min, double x_max);
  //==== JSFillHist : 2D
  std::map< TString, std::map<TString, TH2D*> > JSmaphist_TH2D;
  TH2D* JSGetHist2D(TString suffix, TString histname);
  void JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double x_min, double x_max,
                  int n_biny, double y_min, double y_max);
  void JSFillHist(TString suffix, TString histname,
                  double value_x, double value_y,
                  double weight,
                  int n_binx, double *xbins,
                  int n_biny, double *ybins);

  virtual void WriteHist();
  
  //==== Output rootfile 
  void SwitchToTempDir();
  TFile *outfile;
  TTree *fEventTree = nullptr;
  void SetOutfilePath(TString outname);

};

namespace pi{  const unsigned int nIntTypes = 9;
  const char intTypeName[nIntTypes+1][100] = {"Data",
					      "PiInel",
					      "PiElas",
					      "Muon",
					      "misID:cosmic",
					      "misID:p",
					      "misID:pi",
					      "misID:mu",
					      "misID:e/#gamma",
					      "misID:other"};
  enum intType{
    kData,
    kPiInel,
    kPiElas,
    kMuon,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}

namespace pi2{  const unsigned int nIntTypes = 14;
  const char intTypeName[nIntTypes+1][100] = {"Data",
					      "PiElas",
					      "PiRes",
					      "PiQE",
					      "PiDCEX",
                                              "PiABS",
                                              "PiCEX",
					      "Muon",
					      "misID:cosmic",
					      "misID:p",
					      "misID:pi",
					      "misID:mu",
					      "misID:e/#gamma",
					      "misID:other"};
  enum intType{
    kData,
    kPiElas,
    kPiRes,
    kPiQE,
    kPiDCEX,
    kPiABS,
    kPiCEX,
    kMuon,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}

namespace pitrue{  const unsigned int nIntTypes = 5;
  const char intTypeName[nIntTypes+1][100] = {"Elas",
					      "PiQE",
                                              "PiABS",
                                              "PiCEX",
					      "Other"};
  enum intType{
    kElas,
    kQE,
    kABS,
    kCEX,
    kOther
  };
}

namespace p{
  const unsigned int nIntTypes = 8;
  const char intTypeName[nIntTypes+1][100] = {"Data",
					      "PInel",
					      "PElas",
					      "misID:cosmic",
					      "misID:p",
					      "misID:pi",
					      "misID:mu",
					      "misID:e/#gamma",
					      "misID:other"};
  enum intType{
    kData,
    kPInel,
    kPElas,
    kMIDcosmic,
    kMIDp,
    kMIDpi,
    kMIDmu,
    kMIDeg,
    kMIDother
  };
}
#endif
