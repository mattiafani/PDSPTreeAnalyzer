#ifndef pionqe0p5_h
#define pionqe0p5_h

#include "AnalyzerCore.h"

class pionqe0p5 : public AnalyzerCore {
   public:
    void initializeAnalyzer();
    void executeEvent();

    //==================
    // QE variables
    //==================
    double QE_Q2 = -999;
    double QE_KEPi0 = -999;
    double QE_KEPi1 = -999;
    double QE_AngPi = -999;
    double QE_nu = -999;
    double QE_EQE = -999;
    double recoQE_Q2 = -999;
    double recoQE_KEPi0 = -999;
    double recoQE_KEPi1 = -999;
    double recoQE_AngPi = -999;
    double recoQE_nu = -999;
    double recoQE_EQE = -999;

    //==================
    // beam sel
    //==================
    double cos_spec_tpc_at_z = -999.;
    double delta_x_tpc_spec_at_z = -999;
    double delta_y_tpc_spec_at_z = -999;
    void Set_delta_XY_spec_TPC_at_Z(const vector<double>& calo_X, const vector<double>& calo_Y, const vector<double>& calo_Z, double Z_min);

    bool Pass_BeamStartZ(double N_sigma = 2.);
    double beam_start_z_mu_mc = 0.143;
    double beam_start_z_sigma_mc = 0.219;
    double beam_start_z_mu_data = 3.252;
    double beam_start_z_sigma_data = 1.076;

    bool Pass_beam_delta_X_cut(double N_sigma = 2.);
    bool Pass_beam_delta_Y_cut(double N_sigma = 2.);
    double Beam_delta_X_mu_mc = 1.573;
    double Beam_delta_X_sigma_mc = 1.439;
    double Beam_delta_X_mu_data = -2.647;
    double Beam_delta_X_sigma_data = 1.802;
    double Beam_delta_Y_mu_mc = -0.644;
    double Beam_delta_Y_sigma_mc = 1.473;
    double Beam_delta_Y_mu_data = -2.269;
    double Beam_delta_Y_sigma_data = 2.314;

    double Beam_delta_X_at_z10_mu_mc = -1.762;
    double Beam_delta_X_at_z10_sigma_mc = 1.737;
    double Beam_delta_X_at_z10_mu_data = 3.300;
    double Beam_delta_X_at_z10_sigma_data = 2.107;
    double Beam_delta_Y_at_z10_mu_mc = 0.768;
    double Beam_delta_Y_at_z10_sigma_mc = 1.833;
    double Beam_delta_Y_at_z10_mu_data = 2.498;
    double Beam_delta_Y_at_z10_sigma_data = 2.449;

    //==================
    // beam reco
    //==================
    double KE_ff_subt = 30.;
    void MuonKELoss(TString beam_selec_str, double weight);
    double GetBeamRRatZ10cm(const vector<double>& ResRange, const vector<double>& calo_Z);

    //==================
    // daughter selection
    //==================
    std::vector<Daughter> SelectLooseChargedPions(const vector<Daughter>& in);
    std::vector<Daughter> SelectLooseNeutralPions(const vector<Daughter>& in);

    void StudySecondaryNeutralPions(const vector<Daughter> daughters_all);
    void FillRecoPionPlots(TString daughter_sec_str, const vector<Daughter> pions, double weight);

    //==================
    // true beam study
    //==================
    void FillTrueBeamPlots(TString sel_str, const vector<Daughter> pions, double weight);
    void FillBeamTrueTrajPlots(TString sel_str, TString n_pi_str, double trklen_reco_over_truth, double weight);

    //==================
    // thin-slice cross-section
    //==================
    // void FillXsecHistograms(double weight);
    // double Get_true_tpc_len();
    // bool IsInelasticSignal_reco();
    // bool IsInelasticSignal_true();

    // True track length in the TPC, computed from true_beam_traj_Z/endZ
    double Get_true_tpc_len();

    // Returns true if pi_type is a signal inelastic category (reco side)
    bool IsInelasticSignal_reco();

    // Returns true if pi_truetype is a signal inelastic category (true side)
    bool IsInelasticSignal_true();

    // Main histogram-filling method — call at end of executeEvent() after all cuts
    void FillXsecHistograms(double weight);

    pionqe0p5();
    ~pionqe0p5();
};

#endif
