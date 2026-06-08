#include "TLorentzVector.h"
#include "pionqe0p5.h"

static const double kSliceThickness = 2.5; //5.;  // 10.;  // cm
static const int kNSlices = 100; // 50;

void pionqe0p5::initializeAnalyzer() {
    cout << "[[PionAnalyzer::initializeAnalyzer]] Beam Momentum : " << Beam_Momentum << endl;
    debug_mode = true;
    debug_mode = false;
}

void pionqe0p5::executeEvent() {
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
    KE_ff_reco = KE_beam_inst - KE_ff_subt;  // -- checked upstream KE loss using stopping muons, central value is 30 MeV but we need to assign syst. uncert. on it
    KE_end_reco = -999.;
    E_end_reco = KE_end_reco + mass_beam;

    double P_reweight = 1.;
    if (!IsData) P_reweight = MCCorr->MomentumReweight_SF("Pion_PID_P_0p5", P_beam_inst, 0.);

    // -- 1. Beam instruments
    if (P_beam_inst < 420. || P_beam_inst > 600.) return;
    // if(!PassBeamMomentumWindowCut()) return;

    if (!Pass_Beam_PID(211)) return;
    FillHist("beam_cut_flow_" + pi_type_str, 0.5, 1., 20, 0., 20.);

    if (!PassBeamScraperCut()) return;
    FillHist("beam_cut_flow_" + pi_type_str, 1.5, 1., 20, 0., 20.);

    if (evt.reco_beam_calo_wire->empty()) return;
    Set_delta_XY_spec_TPC_at_Z((*evt.reco_beam_calo_X), (*evt.reco_beam_calo_Y), (*evt.reco_beam_calo_Z), 10.);
    FillHist("beam_cut_flow_" + pi_type_str, 2.5, 1., 20, 0., 20.);

    if (evt.reco_beam_type != pandora_slice_pdg) return;
    FillHist("beam_cut_flow_" + pi_type_str, 3.5, 1., 20, 0., 20.);

    if (evt.reco_beam_calo_endZ < 10.) return;
    double rr_at_z10cm = GetBeamRRatZ10cm((*evt.reco_beam_resRange_SCE), (*evt.reco_beam_calo_Z));
    KE_end_reco = map_BB[211]->KEAtLength(KE_ff_reco, rr_at_z10cm);
    E_end_reco = KE_end_reco + mass_beam;
    FillHist("beam_cut_flow_" + pi_type_str, 4.5, 1., 20, 0., 20.);

    if (!Pass_beam_delta_X_cut(2.)) return;
    if (!Pass_beam_delta_Y_cut(2.)) return;
    FillHist("beam_cut_flow_" + pi_type_str, 5.5, 1., 20, 0., 20.);

    if (chi2_proton > 300. || chi2_proton < 140.) return;
    FillHist("beam_cut_flow_" + pi_type_str, 6.5, 1., 20, 0., 20.);

    vector<Daughter> daughters_all = GetAllDaughters();
    vector<Daughter> loose_charged_pions = SelectLooseChargedPions(daughters_all);
    vector<Daughter> loose_neutral_pions = SelectLooseNeutralPions(daughters_all);

    // == Check allshower and alltrack multiplicities
    if ((*evt.reco_daughter_allTrack_ID).size() != (*evt.reco_daughter_allShower_energy).size()) {
        cout << Form("N(daugh. all track): %zu, N(daugh. all shower): %zu", (*evt.reco_daughter_allTrack_ID).size(), (*evt.reco_daughter_allShower_energy).size()) << endl;
    }

    // == Study secondary particle selection
    StudySecondaryNeutralPions(daughters_all);

    // == Study secondary particle multiplicities
    FillRecoPionPlots("loose_charged_pion", loose_charged_pions, P_reweight);

    int N_loose_charged = (int)loose_charged_pions.size();
    int N_loose_neutral = (int)loose_neutral_pions.size();
    FillXsecHistograms(P_reweight, N_loose_charged, N_loose_neutral);
}

void pionqe0p5::StudySecondaryNeutralPions(const vector<Daughter> daughters_all) {
    vector<TrueDaughter> true_daughters_all = GetAllTrueDaughters();
    bool has_nupi = false;

    for (unsigned int i = 0; i < true_daughters_all.size(); i++) {
        TrueDaughter this_true_daughter = true_daughters_all.at(i);
        int this_PDG = this_true_daughter.PDG();
        double this_startP = this_true_daughter.startP() * 1000.;
        // cout << Form("%d, PDG: %d, startP: %f", i, this_true_daughter.PDG(), this_true_daughter.startP() * 1000.) << endl;
        if (abs(this_PDG) == 211) {
            JSFillHist("true_sec", "true_sec_charged_pion_start_P", this_startP, 1., 1000., 0., 1000.);
        } else if (this_PDG == 111) {
            JSFillHist("true_sec", "true_sec_neutral_pion_start_P", this_startP, 1., 1000., 0., 1000.);
            has_nupi = true;
        } else if (this_PDG == 2212) {
            JSFillHist("true_sec", "true_sec_proton_start_P", this_startP, 1., 1000., 0., 1000.);
        } else
            continue;
    }

    JSFillHist("Daughter_multiplicity", "all", daughters_all.size(), 1., 10., -0.5, 9.5);
    if (has_nupi) JSFillHist("Daughter_multiplicity", "all_has_nupi", daughters_all.size(), 1., 10., -0.5, 9.5);

    for (unsigned int i = 0; i < daughters_all.size(); i++) {
        Daughter this_daughter = daughters_all.at(i);
        int this_PdgID = this_daughter.PFP_true_byHits_PDG();
        TString particle_str = "";
        if (this_PdgID == 2212)
            particle_str = "proton";
        else if (abs(this_PdgID) == 211)
            particle_str = "chpion";
        else if (abs(this_PdgID) == 13)
            particle_str = "muon";
        // else if(this_PdgID == 111)  particle_str = "nupion";
        else
            particle_str = "other";

        if (this_PdgID == 111) cout << "nupi daughter" << endl;

        JSFillHist("Daughter_" + particle_str, particle_str + "_Beam_Dist_allShower", this_daughter.Beam_Dist_allShower(), 1., 100., 0., 100.);
        JSFillHist("Daughter_" + particle_str, particle_str + "_allShower_energy", this_daughter.allShower_energy(), 1., 1000., 0., 1000.);
        JSFillHist("Daughter_" + particle_str, particle_str + "_PFP_trackScore", this_daughter.PFP_trackScore(), 1., 100., 0., 1.);

        if (has_nupi) {
            JSFillHist("Daughter_" + particle_str, particle_str + "_PFP_trackScore_has_nupi", this_daughter.PFP_trackScore(), 1., 100., 0., 1.);
        }
    }
}

void pionqe0p5::FillRecoPionPlots(TString daughter_sec_str, const vector<Daughter> pions, double weight) {
    JSFillHist(daughter_sec_str, daughter_sec_str + "_N_reco_pions_" + pi_type_str, pions.size(), weight, 10., -0.5, 9.5);
    if (pions.size() == 0) {
        JSFillHist(daughter_sec_str, daughter_sec_str + "_Beam_KE_end_0pi_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);
    } else {
        JSFillHist(daughter_sec_str, daughter_sec_str + "_Beam_KE_end_least1pi_" + pi_type_str, KE_end_reco, weight, 2000., 0., 2000.);
    }
}

void pionqe0p5::FillTrueBeamPlots(TString sel_str, const vector<Daughter> pions, double weight) {
    // == study track breaking
    if (!IsData) {
        double true_beam_len = Get_true_beamlen();
        double trklen_reco_over_truth = -1.;
        if (true_beam_len > 0.) trklen_reco_over_truth = evt.reco_beam_alt_len / true_beam_len;

        double KE_ff_true = Get_true_ffKE();
        double trklen_KE_ff_true = -9999.;
        if (KE_ff_true > 0.) trklen_KE_ff_true = map_BB[211]->RangeFromKE(KE_ff_true);

        double trklen_reco_over_KE_ff = evt.reco_beam_alt_len / trklen_KE_ff_true;

        JSFillHist(sel_str, sel_str + "_trklen_reco_over_truth_" + pi_type_str, trklen_reco_over_truth, weight, 1000., 0., 10.);
        JSFillHist(sel_str, sel_str + "_trklen_reco_over_KE_ff_" + pi_type_str, trklen_reco_over_KE_ff, weight, 1000., 0., 10.);

        TString n_pi_str = "";
        if (pions.size() == 0)
            n_pi_str = "0pi";
        else
            n_pi_str = "least1pi";

        JSFillHist(sel_str, sel_str + "_trklen_reco_over_truth_" + n_pi_str + "_" + pi_type_str, trklen_reco_over_truth, weight, 1000., 0., 10.);
        JSFillHist(sel_str, sel_str + "_KE_end_reco_vs_trklen_reco_over_truth_" + n_pi_str + "_" + pi_type_str, KE_end_reco, trklen_reco_over_truth, weight, 10., 0., 500., 150., 0., 1.5);
        JSFillHist(sel_str, sel_str + "_KE_end_reco_vs_trklen_reco_over_KE_ff_" + n_pi_str + "_" + pi_type_str, KE_end_reco, trklen_reco_over_KE_ff, weight, 10., 0., 500., 150., 0., 1.5);

        FillBeamTrueTrajPlots(sel_str, n_pi_str, trklen_reco_over_truth, 1.);
    }
}

void pionqe0p5::FillBeamTrueTrajPlots(TString sel_str, TString n_pi_str, double trklen_reco_over_truth, double weight) {
    if (abs(evt.true_beam_PDG) != 13 && abs(evt.true_beam_PDG) != 211 & abs(evt.true_beam_PDG) != 2212) return;
    int start_idx = -1;
    for (int i_true_hit = 0; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++) {
        if ((*evt.true_beam_traj_Z).at(i_true_hit) >= 0) {
            start_idx = i_true_hit - 1;
            if (start_idx < 0) start_idx = -1;
            break;
        }
    }

    double min_cos_seg = 1.1;
    if (start_idx >= 0) {
        for (int i_true_hit = start_idx + 2; i_true_hit < (*evt.true_beam_traj_Z).size() - 1; i_true_hit++) {
            TVector3 this_vec((*evt.true_beam_traj_X).at(i_true_hit + 1) - (*evt.true_beam_traj_X).at(i_true_hit),
                              (*evt.true_beam_traj_Y).at(i_true_hit + 1) - (*evt.true_beam_traj_Y).at(i_true_hit),
                              (*evt.true_beam_traj_Z).at(i_true_hit + 1) - (*evt.true_beam_traj_Z).at(i_true_hit));

            TVector3 prev_vec((*evt.true_beam_traj_X).at(i_true_hit) - (*evt.true_beam_traj_X).at(i_true_hit - 1),
                              (*evt.true_beam_traj_Y).at(i_true_hit) - (*evt.true_beam_traj_Y).at(i_true_hit - 1),
                              (*evt.true_beam_traj_Z).at(i_true_hit) - (*evt.true_beam_traj_Z).at(i_true_hit - 1));

            double this_cos_theta_seg = this_vec.Dot(prev_vec) / (this_vec.Mag() * prev_vec.Mag());
            if (this_cos_theta_seg < min_cos_seg) min_cos_seg = this_cos_theta_seg;
        }
    }

    JSFillHist(sel_str, sel_str + "_trklen_reco_over_truth_vs_min_cos_seg_" + pi_type_str, trklen_reco_over_truth, min_cos_seg, weight, 20., 0.5, 1.5, 100., -1., 1.);
    JSFillHist(sel_str, sel_str + "_trklen_reco_over_truth_vs_min_cos_seg_" + n_pi_str + "_" + pi_type_str, trklen_reco_over_truth, min_cos_seg, weight, 20., 0.5, 1.5, 100., -1., 1.);
}

double pionqe0p5::GetBeamRRatZ10cm(const vector<double>& ResRange, const vector<double>& calo_Z) {
    double out = -1.;
    if (ResRange.size() < 1 || calo_Z.size() < 1) return out;

    int this_N_calo = calo_Z.size();
    // cout << "ResRange.size(): " << ResRange.size() << ", calo_Z.size(): " << calo_Z.size() << endl; // -- confirmed that the two vectors have exactly the same sizes

    for (int i = 1; i < this_N_calo; i++) {
        // cout << "ResRange " << i << ": " << ResRange.at(i) << ", calo_Z: " << calo_Z.at(i) << endl; // -- confirmed that rr is in descending order
        double this_calo_Z = calo_Z.at(i);
        if (this_calo_Z > 10.) {
            double this_rr = ResRange.at(i);
            double prev_rr = ResRange.at(i - 1);
            double prev_calo_Z = calo_Z.at(i - 1);

            double this_ratio = (10. - prev_calo_Z) / (this_calo_Z - prev_calo_Z);
            double this_delta_rr = prev_rr - this_rr;
            out = prev_rr - this_delta_rr * this_ratio;
            // cout << "this_rr_at_z10cm: " << this_rr_at_z10cm << endl; // -- confirmed that the output value is reasonable
            break;
        }
    }

    return out;
}

bool pionqe0p5::Pass_BeamStartZ(double N_sigma) {
    double this_beam_start_Z = evt.reco_beam_calo_startZ;
    bool out = false;
    if (IsData) {
        out = fabs(beam_start_z_mu_data - this_beam_start_Z) < (N_sigma * beam_start_z_sigma_data);
    } else {
        out = fabs(beam_start_z_mu_mc - this_beam_start_Z) < (N_sigma * beam_start_z_sigma_mc);
    }

    return out;
}

bool pionqe0p5::Pass_beam_delta_X_cut(double N_sigma) {
    bool out = false;
    double Beam_delta_X_over_sigma = (delta_x_tpc_spec_at_z - Beam_delta_X_at_z10_mu_data) / Beam_delta_X_at_z10_sigma_data;
    if (!IsData) Beam_delta_X_over_sigma = (delta_x_tpc_spec_at_z - Beam_delta_X_at_z10_mu_mc) / Beam_delta_X_at_z10_sigma_mc;

    if (fabs(Beam_delta_X_over_sigma) < N_sigma) out = true;

    return out;
}

bool pionqe0p5::Pass_beam_delta_Y_cut(double N_sigma) {
    bool out = false;
    double Beam_delta_Y_over_sigma = (delta_y_tpc_spec_at_z - Beam_delta_Y_at_z10_mu_data) / Beam_delta_Y_at_z10_sigma_data;
    if (!IsData) Beam_delta_Y_over_sigma = (delta_y_tpc_spec_at_z - Beam_delta_Y_at_z10_mu_mc) / Beam_delta_Y_at_z10_sigma_mc;

    if (fabs(Beam_delta_Y_over_sigma) < N_sigma) out = true;

    return out;
}

void pionqe0p5::Set_delta_XY_spec_TPC_at_Z(const vector<double>& calo_X, const vector<double>& calo_Y, const vector<double>& calo_Z, double Z_min) {
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

std::vector<Daughter> pionqe0p5::SelectLooseChargedPions(const vector<Daughter>& in) {
    vector<Daughter> out;
    double cut_trackScore = 0.5;
    double cut_chi2_proton = 60.;
    double cut_trk_len_upper = 180.;
    double cut_trk_len_lower = 10.;
    for (unsigned int i = 0; i < in.size(); i++) {
        Daughter this_in = in.at(i);
        double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
        if (this_in.PFP_trackScore() > cut_trackScore && this_chi2 > cut_chi2_proton && this_in.allTrack_alt_len() < cut_trk_len_upper && this_in.allTrack_alt_len() > cut_trk_len_lower) out.push_back(this_in);
    }
    return out;
}

std::vector<Daughter> pionqe0p5::SelectLooseNeutralPions(const vector<Daughter>& in) {
    vector<Daughter> out;
    double cut_trackScore = 0.3;
    for (unsigned int i = 0; i < in.size(); i++) {
        Daughter this_in = in.at(i);
        if (this_in.PFP_trackScore() < cut_trackScore) out.push_back(this_in);
    }
    return out;
}

double pionqe0p5::Get_true_tpc_len() {
    if (IsData) return -1.;
    if (!evt.true_beam_traj_Z || evt.true_beam_traj_Z->empty()) return -1.;

    // Use the 3D trajectory PATH length from the TPC front face (Z = 0) to the
    // end point, consistent with the reco side (reco_beam_alt_len is a 3D length)
    // and with the thin-slice formula's assumption that one slice = 5 cm of argon
    // traversed. The previous version used the Z-projection (true_beam_endZ - z_ff),
    // which underestimates the argon path by ~1/cos(theta) and biased the truth
    // thin-slice cross section high, growing with depth as the track scatters.
    double true_len = Get_true_beamlen();  // accumulated 3D length from Z = 0
    if (true_len < 0.) return -1.;
    return true_len;
}

// Reco-side signal definitions for data return true so every reco interactor
// fills N_int; the MC purity correction handles channel decomposition.
// MC-side definitions gate on pi_truetype.
bool pionqe0p5::IsAbsSignal_reco() { return IsData ? true : (pi_truetype == pitrue::kABS); }
bool pionqe0p5::IsAbsSignal_true() { return (pi_truetype == pitrue::kABS); }

bool pionqe0p5::IsCexSignal_reco() { return IsData ? true : (pi_truetype == pitrue::kCEX); }
bool pionqe0p5::IsCexSignal_true() { return (pi_truetype == pitrue::kCEX); }

bool pionqe0p5::IsOtherSignal_reco() {
    return IsData ? true : (pi_truetype == pitrue::kQE || pi_truetype == pitrue::kOther);
}
bool pionqe0p5::IsOtherSignal_true() {
    return (pi_truetype == pitrue::kQE || pi_truetype == pitrue::kOther);
}

bool pionqe0p5::IsInelSignal_reco() {
    return IsData ? true : (pi_truetype == pitrue::kQE || pi_truetype == pitrue::kABS || pi_truetype == pitrue::kCEX || pi_truetype == pitrue::kOther);
}
bool pionqe0p5::IsInelSignal_true() {
    return (pi_truetype == pitrue::kQE ||
            pi_truetype == pitrue::kABS ||
            pi_truetype == pitrue::kCEX ||
            pi_truetype == pitrue::kOther);
}

// Four parallel xsec channels:
//   inel       : always fills (total inelastic)
//   abslike    : reco N(loose charged) == 0 && N(loose neutral) == 0
//   cexlike    : reco N(loose charged) == 0 && N(loose neutral) >= 1
//   otherlike  : reco N(loose charged) >= 1

void pionqe0p5::FillXsecHistograms(double weight, int N_loose_charged, int N_loose_neutral) {
    double reco_len = evt.reco_beam_alt_len;
    if (reco_len < 0.) return;

    int recoSliceID = (int)(reco_len / kSliceThickness);
    if (recoSliceID >= kNSlices) recoSliceID = kNSlices - 1;

    double range_ff = map_BB[211]->RangeFromKESpline(KE_ff_reco);

    // N_inc fills: every slice from 0 to recoSliceID
    for (int s = 0; s <= recoSliceID; s++) {
        JSFillHist("Xsec", Form("Xsec_N_inc_reco_%d", pi_type),
                   (double)s, weight, kNSlices, -0.5, kNSlices - 0.5);
        double range_at_mid = range_ff - (s + 0.5) * kSliceThickness;
        if (range_at_mid > 0.) {
            double KE_at_mid = map_BB[211]->KEFromRangeSpline(range_at_mid);
            JSFillHist("Xsec", Form("Xsec_KE_reco_slice%d", s),
                       KE_at_mid, weight, 500, 0., 1000.);
        }
    }

    // Reco pool assignment from loose daughter multiplicities
    bool reco_abslike = (N_loose_charged == 0) && (N_loose_neutral == 0);
    bool reco_cexlike = (N_loose_charged == 0) && (N_loose_neutral >= 1);
    bool reco_otherlike = (N_loose_charged >= 1);

    // N_int fills: total inelastic always; exclusive channels gated on reco pool
    if (IsInelSignal_reco())
        JSFillHist("Xsec", Form("Xsec_N_int_reco_inel_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
    if (reco_abslike && IsAbsSignal_reco())
        JSFillHist("Xsec", Form("Xsec_N_int_reco_abslike_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
    if (reco_cexlike && IsCexSignal_reco())
        JSFillHist("Xsec", Form("Xsec_N_int_reco_cexlike_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
    if (reco_otherlike && IsOtherSignal_reco())
        JSFillHist("Xsec", Form("Xsec_N_int_reco_otherlike_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);

    // Purity decomposition (MC only): denominator is all selected events in the
    // reco pool; numerator is events whose true channel matches the pool name.
    if (!IsData) {
        JSFillHist("Xsec", Form("Xsec_purity_all_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        if (IsInelSignal_true())
            JSFillHist("Xsec", Form("Xsec_purity_sig_inel_%d", pi_type),
                       (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);

        if (reco_abslike) {
            JSFillHist("Xsec", Form("Xsec_purity_all_abslike_%d", pi_type),
                       (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
            if (IsAbsSignal_true())
                JSFillHist("Xsec", Form("Xsec_purity_sig_abslike_%d", pi_type),
                           (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        }
        if (reco_cexlike) {
            JSFillHist("Xsec", Form("Xsec_purity_all_cexlike_%d", pi_type),
                       (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
            if (IsCexSignal_true())
                JSFillHist("Xsec", Form("Xsec_purity_sig_cexlike_%d", pi_type),
                           (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        }
        if (reco_otherlike) {
            JSFillHist("Xsec", Form("Xsec_purity_all_otherlike_%d", pi_type),
                       (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
            if (IsOtherSignal_true())
                JSFillHist("Xsec", Form("Xsec_purity_sig_otherlike_%d", pi_type),
                           (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        }
    }

    if (!IsData) {
        double true_len = Get_true_tpc_len();
        if (true_len < 0.) return;

        int trueSliceID = (int)(true_len / kSliceThickness);
        if (trueSliceID >= kNSlices) trueSliceID = kNSlices - 1;

        double KE_ff_true = Get_true_ffKE();
        double range_ff_true = (KE_ff_true > 0.) ? map_BB[211]->RangeFromKESpline(KE_ff_true) : -1.;

        for (int s = 0; s <= trueSliceID; s++) {
            JSFillHist("Xsec", Form("Xsec_N_inc_true_%d", pi_type),
                       (double)s, weight, kNSlices, -0.5, kNSlices - 0.5);
            if (range_ff_true > 0.) {
                double range_at_mid_true = range_ff_true - (s + 0.5) * kSliceThickness;
                if (range_at_mid_true > 0.) {
                    double KE_true_at_mid = map_BB[211]->KEFromRangeSpline(range_at_mid_true);
                    JSFillHist("Xsec", Form("Xsec_KE_true_slice%d_%d", s, pi_type),
                               KE_true_at_mid, weight, 500, 0., 1000.);
                }
            }
        }

        if (IsInelSignal_true())
            JSFillHist("Xsec", Form("Xsec_N_int_true_inel_%d", pi_type),
                       (double)trueSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        if (IsAbsSignal_true())
            JSFillHist("Xsec", Form("Xsec_N_int_true_abslike_%d", pi_type),
                       (double)trueSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        if (IsCexSignal_true())
            JSFillHist("Xsec", Form("Xsec_N_int_true_cexlike_%d", pi_type),
                       (double)trueSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        if (IsOtherSignal_true())
            JSFillHist("Xsec", Form("Xsec_N_int_true_otherlike_%d", pi_type),
                       (double)trueSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);

        // Migration: reco slice vs true slice, reco information.
        // _inel uses any reco signal (i.e. always for data, IsInelSignal_reco for MC);
        // exclusive channels also require true type to match for a clean efficiency matrix.
        if (IsInelSignal_reco())
            JSFillHist("Xsec", Form("Xsec_migration_inel_%d", pi_type),
                       (double)recoSliceID, (double)trueSliceID, weight,
                       kNSlices, -0.5, kNSlices - 0.5,
                       kNSlices, -0.5, kNSlices - 0.5);
        if (reco_abslike && IsAbsSignal_true())
            JSFillHist("Xsec", Form("Xsec_migration_abslike_%d", pi_type),
                       (double)recoSliceID, (double)trueSliceID, weight,
                       kNSlices, -0.5, kNSlices - 0.5,
                       kNSlices, -0.5, kNSlices - 0.5);
        if (reco_cexlike && IsCexSignal_true())
            JSFillHist("Xsec", Form("Xsec_migration_cexlike_%d", pi_type),
                       (double)recoSliceID, (double)trueSliceID, weight,
                       kNSlices, -0.5, kNSlices - 0.5,
                       kNSlices, -0.5, kNSlices - 0.5);
        if (reco_otherlike && IsOtherSignal_true())
            JSFillHist("Xsec", Form("Xsec_migration_otherlike_%d", pi_type),
                       (double)recoSliceID, (double)trueSliceID, weight,
                       kNSlices, -0.5, kNSlices - 0.5,
                       kNSlices, -0.5, kNSlices - 0.5);
    }
}

/////////////////////////////////

pionqe0p5::pionqe0p5() {
}

pionqe0p5::~pionqe0p5() {
}
