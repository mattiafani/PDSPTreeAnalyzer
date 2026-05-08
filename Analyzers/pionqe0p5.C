#include "TLorentzVector.h"
#include "pionqe0p5.h"

static const double kSliceThickness = 5.;  // 10.;  // cm
static const int kNSlices = 50;

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

    // == Study broken tracks
    // FillTrueBeamPlots("Beam_chi2proton", loose_pions, P_reweight);

    // /////////////////////////////
    // // ---- TEMPORARY: verify true_beam_endZ vs last traj point ----
    // static int debug_count = 0;
    // if (!IsData && debug_count < 20 && abs(evt.true_beam_PDG) == 211) {
    //     double last_traj_Z = evt.true_beam_traj_Z->empty() ? -999. : evt.true_beam_traj_Z->back();
    //     cout << Form("Evt %d | true_beam_endZ = %.2f | traj_Z.back() = %.2f | diff = %.3f | endProcess = %s",
    //                  evt.event,
    //                  evt.true_beam_endZ,
    //                  last_traj_Z,
    //                  evt.true_beam_endZ - last_traj_Z,
    //                  evt.true_beam_endProcess->c_str())
    //          << endl;
    //     debug_count++;
    // }
    // // ---- END TEMPORARY ----
    // /////////////////////////////

    // // ---- TEMPORARY 2:
    // static int dbg2 = 0;
    // if (!IsData && dbg2 < 30 && IsInelasticSignal_reco()) {
    //     cout << Form("Evt %d: pi_type=%d pi_truetype=%d endProcess=%s recoLen=%.1f",
    //                  evt.event, pi_type, pi_truetype,
    //                  evt.true_beam_endProcess->c_str(),
    //                  evt.reco_beam_alt_len)
    //          << endl;
    //     dbg2++;
    // }
    // // ---- END TEMPORARY 2 ----
    // /////////////////////////////

    FillXsecHistograms(P_reweight);
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

// ------------------------------------------------------------
// Get_true_tpc_len
//
// Returns the true track length inside the TPC, defined as:
//   from the first trajectory point with Z >= 0 (front face)
//   to true_beam_endZ.
// Returns -1 if the front face crossing cannot be found or if
// the pion ends before entering the TPC.
// ------------------------------------------------------------
double pionqe0p5::Get_true_tpc_len() {
    if (IsData) return -1.;
    if (!evt.true_beam_traj_Z || evt.true_beam_traj_Z->empty()) return -1.;

    // Find the first traj point at or past Z = 0 (TPC front face)
    double z_ff = -1.;
    for (int i = 0; i < (int)evt.true_beam_traj_Z->size(); i++) {
        double z = evt.true_beam_traj_Z->at(i);
        if (z >= 0.) {
            // Interpolate to Z = 0 if needed
            if (i > 0) {
                double z_prev = evt.true_beam_traj_Z->at(i - 1);
                // z_prev < 0 < z: interpolate linearly
                z_ff = 0.;  // we define the entry point as z = 0
            } else {
                z_ff = z;
            }
            break;
        }
    }
    if (z_ff < 0.) return -1.;  // never entered TPC

    double z_end = evt.true_beam_endZ;
    if (z_end < 0.) return -1.;  // interacted before entering TPC

    double true_len = z_end - z_ff;
    if (true_len < 0.) return -1.;

    return true_len;
}

// ------------------------------------------------------------
// IsInelasticSignal_reco
//
// For the reco side: a pion is counted as an inelastic interactor
// if its true type is one of the beam-matched pion categories.
// kPiElas includes both true elastic AND pi→mu (endProcess is
// "pi+Inelastic" but daughter is a muon), so we include it.
// All six pi2 pion categories (1–6) are signal; background
// categories (misID, muon) are not.
// ------------------------------------------------------------

bool pionqe0p5::IsInelasticSignal_reco() {
    if (IsData) return true;  // for data, count all track endings
    return (pi_truetype == pitrue::kQE ||
            pi_truetype == pitrue::kABS ||
            pi_truetype == pitrue::kCEX);
}

// bool pionqe0p5::IsInelasticSignal_reco() {
//     return (pi_type == pi2::kPiElas ||
//             pi_type == pi2::kPiRes ||
//             pi_type == pi2::kPiQE ||
//             pi_type == pi2::kPiDCEX ||
//             pi_type == pi2::kPiABS ||
//             pi_type == pi2::kPiCEX);
// }

// ------------------------------------------------------------
// IsInelasticSignal_true
//
// For the true side: inelastic = QE, ABS, or CEX.
// kElas covers true elastic scattering and range-out → NOT signal.
// kOther covers non-pion beam particles → NOT signal.
// ------------------------------------------------------------
bool pionqe0p5::IsInelasticSignal_true() {
    return (pi_truetype == pitrue::kQE ||
            pi_truetype == pitrue::kABS ||
            pi_truetype == pitrue::kCEX);
}

// ------------------------------------------------------------
// FillXsecHistograms
//
// Called after all beam selection cuts pass.
// Fills the following histograms (all in TDirectory "Xsec"):
//
//  Reco side (per pi2 category, suffix = "_<pi_type>"):
//   Xsec_N_inc_reco_<N>   : slice ID for incident pions
//   Xsec_N_int_reco_<N>   : slice ID for interacting pions (reco signal)
//   Xsec_KE_reco_slice<S> : KE distribution of incident pions in slice S
//
//  True side (MC only, same category suffix):
//   Xsec_N_inc_true_<N>   : true slice ID for incident pions
//   Xsec_N_int_true_<N>   : true slice ID for inelastic pions
//   Xsec_KE_true_slice<S>_<N> : true KE distribution in true slice S
//
//  Purity (MC only):
//   Xsec_purity_all_<N>   : all events after cuts, in reco slice ID
//   Xsec_purity_signal_<N>: signal events (beam-matched pion) in reco slice ID
//
//  Migration (MC only):
//   Xsec_migration_<N>    : 2D(recoSliceID, trueSliceID) for signal events
// ------------------------------------------------------------
void pionqe0p5::FillXsecHistograms(double weight) {
    // ----------------------------------------------------------
    // 1. Compute reco slice ID
    //    The reco track length is reco_beam_alt_len.
    //    SliceID = floor(track_length / kSliceThickness)
    //    Clamp to [0, kNSlices-1] to avoid out-of-range fills.
    // ----------------------------------------------------------
    double reco_len = evt.reco_beam_alt_len;
    if (reco_len < 0.) return;  // sanity check

    int recoSliceID = (int)(reco_len / kSliceThickness);
    if (recoSliceID >= kNSlices) recoSliceID = kNSlices - 1;

    // KE at the front face (already computed in executeEvent)
    double KE_reco = KE_ff_reco;

    // KE at each slice midpoint (for the KE-axis calibration histograms)
    // We approximate it as KE_ff_reco minus the energy lost traversing
    // to the midpoint of the slice, using Bethe-Bloch range lookup.
    // KE_at_slice_mid(s) = KE from range: range(KE_ff) - (s + 0.5) * t
    double range_ff = map_BB[211]->RangeFromKESpline(KE_ff_reco);

    // ----------------------------------------------------------
    // 2. Fill N_inc (reco): pion is incident on every slice from
    //    0 up to and including recoSliceID.
    // ----------------------------------------------------------
    for (int s = 0; s <= recoSliceID; s++) {
        JSFillHist("Xsec", Form("Xsec_N_inc_reco_%d", pi_type),
                   (double)s, weight, kNSlices, -0.5, kNSlices - 0.5);

        // KE at midpoint of slice s
        double range_at_mid = range_ff - (s + 0.5) * kSliceThickness;
        if (range_at_mid > 0.) {
            double KE_at_mid = map_BB[211]->KEFromRangeSpline(range_at_mid);
            JSFillHist("Xsec", Form("Xsec_KE_reco_slice%d", s),
                       KE_at_mid, weight, 500, 0., 1000.);
        }
    }

    // ----------------------------------------------------------
    // 3. Fill N_int (reco): only the last slice, and only if
    //    the pion is classified as an inelastic interactor.
    //    "Inelastic" here = any of the 6 beam-matched pion types
    //    (kPiElas through kPiCEX). Pions that range out before
    //    the last slice boundary will not be classified as
    //    inelastic here — they have pi_type = kPiElas with
    //    endProcess != "pi+Inelastic", but that still maps to
    //    kPiElas. We include them because the cross-section
    //    formula accounts for them via the denominator structure;
    //    the unfolding and purity correction will handle it.
    //    Background categories (misID, muon) are NOT filled here.
    // ----------------------------------------------------------
    // if (IsInelasticSignal_reco()) {
    //     JSFillHist("Xsec", Form("Xsec_N_int_reco_%d", pi_type),
    //                (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
    // }
    bool fill_int = IsData ? true : IsInelasticSignal_reco();
    if (fill_int) {
        JSFillHist("Xsec", Form("Xsec_N_int_reco_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
    }

    // ----------------------------------------------------------
    // 4. Purity histograms (MC only)
    //    purity_all: all events (any category) in this reco slice
    //    purity_signal: only beam-matched pion events
    // ----------------------------------------------------------
    if (!IsData) {
        JSFillHist("Xsec", Form("Xsec_purity_all_%d", pi_type),
                   (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);

        if (IsInelasticSignal_reco()) {
            JSFillHist("Xsec", Form("Xsec_purity_signal_%d", pi_type),
                       (double)recoSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        }
    }

    // ----------------------------------------------------------
    // 5. True side (MC only)
    // ----------------------------------------------------------
    if (!IsData) {
        double true_len = Get_true_tpc_len();
        if (true_len < 0.) return;  // pion didn't enter TPC

        int trueSliceID = (int)(true_len / kSliceThickness);
        if (trueSliceID >= kNSlices) trueSliceID = kNSlices - 1;

        double KE_ff_true = Get_true_ffKE();
        double range_ff_true = (KE_ff_true > 0.) ? map_BB[211]->RangeFromKESpline(KE_ff_true) : -1.;

        // -- N_inc (true)
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

        // -- N_int (true): only fill if this is a true inelastic interaction
        if (IsInelasticSignal_true()) {
            JSFillHist("Xsec", Form("Xsec_N_int_true_%d", pi_type),
                       (double)trueSliceID, weight, kNSlices, -0.5, kNSlices - 0.5);
        }

        // -- Migration matrix: reco vs true slice ID
        //    Only fill for beam-matched pion events (signal categories)
        if (!IsData && IsInelasticSignal_reco()) {
            JSFillHist("Xsec", Form("Xsec_migration_%d", pi_type),
                       (double)recoSliceID, (double)trueSliceID, weight,
                       kNSlices, -0.5, kNSlices - 0.5,
                       kNSlices, -0.5, kNSlices - 0.5);
        }

    }  // end MC-only block
}

/////////////////////////////////

pionqe0p5::pionqe0p5() {
}

pionqe0p5::~pionqe0p5() {
}
