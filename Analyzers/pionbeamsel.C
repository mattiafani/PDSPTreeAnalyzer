#include "TLorentzVector.h"
#include "pionbeamsel.h"

void pionbeamsel::initializeAnalyzer() {
    cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
    debug_mode = true;
    debug_mode = false;
}

void pionbeamsel::executeEvent() {
    pi_type = GetPi2ParType();
    pi_type_str = Form("%d", pi_type);
    FillHist("beam_cut_flow", 0.5, 1., 20, 0., 20.);

    pi_truetype = GetPiTrueType();
    // if (pi_truetype == pitrue::kQE) FillQEPlots("QE_All");

    P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
    KE_beam_inst = map_BB[211]->MomentumtoKE(P_beam_inst);
    exp_trk_len_beam_inst = map_BB[211]->RangeFromKESpline(KE_beam_inst);
    trk_len_ratio = evt.reco_beam_alt_len / exp_trk_len_beam_inst;
    mass_beam = 139.57;
    KE_ff_reco = KE_beam_inst - 30.;  // -- checked upstream KE loss using stopping muons, central value is 30 MeV but we need to assign syst. uncert. on it
    KE_end_reco = -999.;
    E_end_reco = KE_end_reco + mass_beam;

    if (P_beam_inst > 522. && P_beam_inst < 523. && pi_type == 1) {
        cout << Form("P 522 MeV/c!!! (run,subrun,evt) = (%d, %d, %d), true_beam_P(X,Y,Z) = (%f, %f, %f)", evt.run, evt.subrun, evt.event, evt.true_beam_startPx * 1000., evt.true_beam_startPy * 1000., evt.true_beam_startPz * 1000.) << endl;
    }

    if (P_beam_inst > 483. && P_beam_inst < 484. && pi_type == 3) {
        cout << Form("P 483 MeV/c!!! (run,subrun,evt) = (%d, %d, %d), true_beam_P(X,Y,Z) = (%f, %f, %f)", evt.run, evt.subrun, evt.event, evt.true_beam_startPx * 1000., evt.true_beam_startPy * 1000., evt.true_beam_startPz * 1000.) << endl;
    }

    double P_reweight = 1.;
    if (!IsData) P_reweight = MCCorr->MomentumReweight_SF("Pion_PID_P_0p5", P_beam_inst, 0.);
    FillHist("beam_inst_P_x2", evt.beam_inst_P * 2. * 1000., 1., 100., 0., 1000.);

    // -- 1. Beam instruments

    // CUT: Momentum window
    if (P_beam_inst < 420. || P_beam_inst > 600.) return;
    // if(!PassBeamMomentumWindowCut()) return;

    // CUT: Beam_PID
    if (!Pass_Beam_PID(211)) return;
    FillBeamPlots("Beam_PID", P_reweight);

    // CUT: Beam_scraper
    if (!PassBeamScraperCut()) return;
    FillHist("beam_cut_flow", 1.5, 1., 20, 0., 20.);
    FillBeamPlots("Beam_scraper", P_reweight);
    if (trk_len_ratio > 1.) {
        MuonKELoss("Beam_scraper", 1.);
    }

    // CUT: Beam_calo || Beam_coll_hits || Beam_calo_wire
    if (evt.reco_beam_calo_wire->empty()) return;
    Set_delta_XY_spec_TPC_at_Z((*evt.reco_beam_calo_X), (*evt.reco_beam_calo_Y), (*evt.reco_beam_calo_Z), Z_beam_end_cut);
    FillBeamPlots("Beam_collhits", P_reweight);

    // CUT: Beam_reco_trk
    if (evt.reco_beam_type != pandora_slice_pdg) return;
    // Uncomment this line to also includes a trk_len_ratio cut -> TODO: make this a real cut, properly integrated in the code
    // if (trk_len_ratio >= 0.9) return;
    FillBeamPlots("Beam_recotrk", P_reweight);

    // double true_ff_X, true_ff_Y, true_ff_Z;
    // if(!IsData) get_true_XYZ_at_FF(true_ff_X, true_ff_Y, true_ff_Z);
    /*
    if(evt.reco_beam_calo_wire->empty()){
      cout << Form("(run,subrun,evt) = (%d, %d, %d): ", evt.run, evt.subrun, evt.event) << ", PDG: " << evt.true_beam_PDG << ", true_beam_startP: " << evt.true_beam_startP * 1000. <<
        Form(", true_beam_FF(X,Y,Z): (%f, %f, %f)", true_ff_X, true_ff_Y, true_ff_Z) << Form(", true_beam_end(X,Y,Z): (%f, %f, %f)", evt.true_beam_endX, evt.true_beam_endY, evt.true_beam_endZ) << ", evt.reco_beam_type: " << evt.reco_beam_type << ", evt.reco_beam_calo_wire->empty(): " <<
        evt.reco_beam_calo_wire->empty() << ", reco_beam_dQdX_SCE.size: " << (*evt.reco_beam_dQdX_SCE).size() << ", reco_beam_PFP_nHits: " << evt.reco_beam_PFP_nHits << ", " << (*evt.true_beam_endProcess) << endl;
    }
    */
    // if(!Pass_BeamStartZ(2.)) return;
    // FillBeamPlots("Beam_startZ", P_reweight);

    // CUT: Beam_endZ ???
    if (evt.reco_beam_calo_endZ < Z_beam_end_cut) return;
    double rr_at_z_cut = GetBeamRRatZ((*evt.reco_beam_resRange_SCE), (*evt.reco_beam_calo_Z), Z_beam_end_cut);
    KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, rr_at_z_cut);
    E_end_reco = KE_end_reco + mass_beam;
    FillBeamPlots("Beam_endZ", P_reweight);

    // CUT: Beam_deltaXY
    if (!Pass_beam_delta_X_cut(2.)) return;
    if (!Pass_beam_delta_Y_cut(2.)) return;
    FillBeamPlots("Beam_deltaXY", P_reweight);

    // CUT: Beam_chi2proton (modified 20260508)
    // if (chi2_proton > 280. || chi2_proton < 160.) return;
    if (chi2_proton > 300. || chi2_proton < 140.) return;
    FillBeamPlots("Beam_chi2proton", P_reweight);
}

void pionbeamsel::FillBeamPlots(TString beam_selec_str, double weight) {
    double beam_reco_as_trk = -2.;
    if (evt.reco_beam_type == pandora_slice_pdg)
        beam_reco_as_trk = 0.5;
    else
        beam_reco_as_trk = -0.5;

    double beam_calo_size = -2.;
    if (!(evt.reco_beam_calo_wire->empty()))
        beam_calo_size = 0.5;
    else
        beam_calo_size = -0.5;

    // == Fit results after calo-size cut
    double Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - beam_start_z_mu_data) / beam_start_z_sigma_data;
    double Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_data) / Beam_delta_X_sigma_data;
    double Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_data) / Beam_delta_Y_sigma_data;
    double Beam_delta_X_at_z10_over_sigma = (delta_x_tpc_spec_at_z - Beam_delta_X_at_z10_mu_data) / Beam_delta_X_at_z10_sigma_data;
    double Beam_delta_Y_at_z10_over_sigma = (delta_y_tpc_spec_at_z - Beam_delta_Y_at_z10_mu_data) / Beam_delta_Y_at_z10_sigma_data;
    if (!IsData) {
        Beam_startZ_over_sigma = (evt.reco_beam_calo_startZ - beam_start_z_mu_mc) / beam_start_z_sigma_mc;
        Beam_delta_X_over_sigma = (delta_X_spec_TPC - Beam_delta_X_mu_mc) / Beam_delta_X_sigma_mc;
        Beam_delta_Y_over_sigma = (delta_Y_spec_TPC - Beam_delta_Y_mu_mc) / Beam_delta_Y_sigma_mc;
        Beam_delta_X_at_z10_over_sigma = (delta_x_tpc_spec_at_z - Beam_delta_X_at_z10_mu_mc) / Beam_delta_X_at_z10_sigma_mc;
        Beam_delta_Y_at_z10_over_sigma = (delta_y_tpc_spec_at_z - Beam_delta_Y_at_z10_mu_mc) / Beam_delta_Y_at_z10_sigma_mc;
        JSFillHist(beam_selec_str, beam_selec_str + "_true_beam_startP", evt.true_beam_startP * 1000., weight, 2000., 0., 2000.);
        JSFillHist(beam_selec_str, beam_selec_str + "_true_beam_startP_" + pi_type_str, evt.true_beam_startP * 1000., weight, 2000., 0., 2000.);
    }

    double Z_dir_sign = evt.reco_beam_calo_endZ - evt.reco_beam_calo_startZ;
    if (Z_dir_sign < 0.)
        Z_dir_sign = -0.5;
    else
        Z_dir_sign = 0.5;

    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX", evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY", evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ", evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ", evt.reco_beam_calo_endZ, weight, 1000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Z_dir_sign", Z_dir_sign, weight, 2., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_beam_inst_X", evt.beam_inst_X, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_beam_inst_Y", evt.beam_inst_Y, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst", P_beam_inst, weight, 2000., 0., 2000.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh", beam_costh, weight, 2000., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton", chi2_proton, weight, 10000., 0., 1000.);

    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_reco_as_trk_" + pi_type_str, beam_reco_as_trk, weight, 2., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_calo_size_" + pi_type_str, beam_calo_size, weight, 2., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startX_" + pi_type_str, evt.reco_beam_calo_startX, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startY_" + pi_type_str, evt.reco_beam_calo_startY, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_" + pi_type_str, evt.reco_beam_calo_startZ, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_Z_dir_sign_" + pi_type_str, Z_dir_sign, weight, 2., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_inst_X_" + pi_type_str, evt.beam_inst_X, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_inst_Y_" + pi_type_str, evt.beam_inst_Y, weight, 10000., -100., 900.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_" + pi_type_str, delta_X_spec_TPC, weight, 2000., -100., 100.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_" + pi_type_str, delta_Y_spec_TPC, weight, 2000., -100., 100.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_x_tpc_spec_at_z_" + pi_type_str, delta_x_tpc_spec_at_z, weight, 2000., -100., 100.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_y_tpc_spec_at_z_" + pi_type_str, delta_y_tpc_spec_at_z, weight, 2000., -100., 100.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_cos_delta_spec_TPC_" + pi_type_str, cos_delta_spec_TPC, weight, 2000., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_cos_spec_tpc_at_z_" + pi_type_str, cos_spec_tpc_at_z, weight, 2000., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_costh_" + pi_type_str, beam_costh, weight, 2000., -1., 1.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_theta_" + pi_type_str, beam_TPC_theta, weight, 5000., -1., 4.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_TPC_phi_" + pi_type_str, beam_TPC_phi, weight, 8000., -4., 4.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_endZ_" + pi_type_str, evt.reco_beam_calo_endZ, weight, 1100., -100., 1000.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_alt_len_" + pi_type_str, evt.reco_beam_alt_len, weight, 1000., 0., 1000.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_P_beam_inst_" + pi_type_str, P_beam_inst, weight, 2000., 0., 2000.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_chi2_proton_" + pi_type_str, chi2_proton, weight, 10000., 0., 1000.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_trk_len_ratio_" + pi_type_str, trk_len_ratio, weight, 1000., 0, 10.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_ff_" + pi_type_str, KE_ff_reco, weight, 2000., 0., 2000.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_end_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);

    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_startZ_over_sigma_" + pi_type_str, Beam_startZ_over_sigma, weight, 2000., -10., 10.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_X_over_sigma, weight, 2000., -10., 10.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_spec_TPC_over_sigma_" + pi_type_str, Beam_delta_Y_over_sigma, weight, 2000., -10., 10.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_X_at_z10_over_sigma_" + pi_type_str, Beam_delta_X_at_z10_over_sigma, weight, 2000., -10., 10.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_delta_Y_at_z10_over_sigma_" + pi_type_str, Beam_delta_Y_at_z10_over_sigma, weight, 2000., -10., 10.);
}

void pionbeamsel::MuonKELoss(TString beam_selec_str, double weight) {
    double this_muon_beam_inst_KE = map_BB[13]->MomentumtoKE(P_beam_inst);

    double this_rr_at_z_cut = GetBeamRRatZ((*evt.reco_beam_resRange_SCE), (*evt.reco_beam_calo_Z), Z_beam_end_cut);
    double this_muon_KE_z_cut = map_BB[13]->KEFromRangeSpline(this_rr_at_z_cut);

    double this_KELoss = this_muon_beam_inst_KE - this_muon_KE_z_cut;
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KELoss_" + pi_type_str, this_KELoss, weight, 2000., -100., 100.);
    JSFillHist(beam_selec_str, beam_selec_str + "_Beam_KE_beam_inst_vs_Beam_KELoss_" + pi_type_str, this_muon_beam_inst_KE, this_KELoss, weight, 100., 0., 1000., 200., -100., 100.);
}

double pionbeamsel::GetBeamRRatZ(const vector<double>& ResRange, const vector<double>& calo_Z, double Z_fv = 30.) {
    double out = -1.;
    if (ResRange.size() < 1 || calo_Z.size() < 1) return out;

    int this_N_calo = calo_Z.size();
    // cout << "ResRange.size(): " << ResRange.size() << ", calo_Z.size(): " << calo_Z.size() << endl; // -- confirmed that the two vectors have exactly the same sizes

    for (int i = 1; i < this_N_calo; i++) {
        // cout << "ResRange " << i << ": " << ResRange.at(i) << ", calo_Z: " << calo_Z.at(i) << endl; // -- confirmed that rr is in descending order
        double this_calo_Z = calo_Z.at(i);
        if (this_calo_Z > Z_fv) {
            double this_rr = ResRange.at(i);
            double prev_rr = ResRange.at(i - 1);
            double prev_calo_Z = calo_Z.at(i - 1);

            double this_ratio = (Z_fv - prev_calo_Z) / (this_calo_Z - prev_calo_Z);
            double this_delta_rr = prev_rr - this_rr;
            out = prev_rr - this_delta_rr * this_ratio;
            break;
        }
    }

    return out;
}

bool pionbeamsel::Pass_BeamStartZ(double N_sigma) {
    double this_beam_start_Z = evt.reco_beam_calo_startZ;
    bool out = false;
    if (IsData) {
        out = fabs(beam_start_z_mu_data - this_beam_start_Z) < (N_sigma * beam_start_z_sigma_data);
    } else {
        out = fabs(beam_start_z_mu_mc - this_beam_start_Z) < (N_sigma * beam_start_z_sigma_mc);
    }

    return out;
}

bool pionbeamsel::Pass_beam_delta_X_cut(double N_sigma) {
    bool out = false;
    double Beam_delta_X_over_sigma = (delta_x_tpc_spec_at_z - Beam_delta_X_at_z10_mu_data) / Beam_delta_X_at_z10_sigma_data;
    if (!IsData) Beam_delta_X_over_sigma = (delta_x_tpc_spec_at_z - Beam_delta_X_at_z10_mu_mc) / Beam_delta_X_at_z10_sigma_mc;

    if (fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

    return out;
}

bool pionbeamsel::Pass_beam_delta_Y_cut(double N_sigma) {
    bool out = false;
    double Beam_delta_Y_over_sigma = (delta_y_tpc_spec_at_z - Beam_delta_Y_at_z10_mu_data) / Beam_delta_Y_at_z10_sigma_data;
    if (!IsData) Beam_delta_Y_over_sigma = (delta_y_tpc_spec_at_z - Beam_delta_Y_at_z10_mu_mc) / Beam_delta_Y_at_z10_sigma_mc;

    if (fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

    return out;
}

void pionbeamsel::Set_delta_XY_spec_TPC_at_Z(const vector<double>& calo_X, const vector<double>& calo_Y, const vector<double>& calo_Z, double Z_min) {
    double X_spec_at_Z = evt.beam_inst_X + Z_min * evt.beam_inst_dirX / evt.beam_inst_dirZ;
    double Y_sepc_at_Z = evt.beam_inst_Y + Z_min * evt.beam_inst_dirY / evt.beam_inst_dirZ;

    size_t fv_index = 0;
    for (; fv_index < calo_Z.size(); ++fv_index) {
        if (calo_Z[fv_index] > Z_min) break;
    }
    if (fv_index == 0) ++fv_index;

    double x1 = calo_X[fv_index];
    double y1 = calo_Y[fv_index];
    double z1 = calo_Z[fv_index];

    double x0 = calo_X[fv_index - 1];
    double y0 = calo_Y[fv_index - 1];
    double z0 = calo_Z[fv_index - 1];

    double xl = calo_X.back();
    double yl = calo_Y.back();
    double zl = calo_Z.back();

    // project to the FV face
    double x = (Z_min - z0) * (x1 - x0) / (z1 - z0) + x0;
    double y = (Z_min - z0) * (y1 - y0) / (z1 - z0) + y0;

    delta_x_tpc_spec_at_z = x - X_spec_at_Z;
    delta_y_tpc_spec_at_z = y - Y_sepc_at_Z;

    double r_end = sqrt((xl - x) * (xl - x) + (yl - y) * (yl - y) + (zl - Z_min) * (zl - Z_min));
    cos_spec_tpc_at_z = ((xl - x) * evt.beam_inst_dirX + (yl - y) * evt.beam_inst_dirY + (zl - Z_min) * evt.beam_inst_dirZ) / r_end;
}

void pionbeamsel::get_true_XYZ_at_FF(double& true_ff_X, double& true_ff_Y, double& true_ff_Z) {
    double this_X = -999.;
    double this_Y = -999.;
    double this_Z = -999.;
    int start_idx = -1;

    for (int i_true_hit = 0; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++) {
        if ((*evt.true_beam_traj_Z).at(i_true_hit) >= 0) {
            start_idx = i_true_hit - 1;
            if (start_idx < 0) start_idx = -1;
            break;
        }
    }
    if (start_idx >= 0) {
        this_X = (*evt.true_beam_traj_X).at(start_idx + 1);
        this_Y = (*evt.true_beam_traj_Y).at(start_idx + 1);
        this_Z = (*evt.true_beam_traj_Z).at(start_idx + 1);
    }

    true_ff_X = this_X;
    true_ff_Y = this_Y;
    true_ff_Z = this_Z;
}

pionbeamsel::pionbeamsel() {
}

pionbeamsel::~pionbeamsel() {
}
