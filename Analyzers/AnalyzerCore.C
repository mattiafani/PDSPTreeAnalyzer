#include "AnalyzerCore.h"
#include "PDSPTree.h"

AnalyzerCore::AnalyzerCore(){
  MaxEvent = -1;
  NSkipEvent = 0;
  LogEvery = 1000;
  MCSample = "";
  Beam_Momentum = -1.;
  Userflags.clear();
  outfile = NULL;
  MCCorr = new MCCorrection();
  G4Xsec = new GEANT4_XS();
  map_BB[13] = new BetheBloch(13);
  map_BB[211] = new BetheBloch(211);
  map_BB[321] = new BetheBloch(321);
  map_BB[2212] = new BetheBloch(2212);
  Fitter = new Fittings();
}

AnalyzerCore::~AnalyzerCore(){

  //=== hist maps
  
  for(std::map< TString, TH1D* >::iterator mapit = maphist_TH1D.begin(); mapit!=maphist_TH1D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH1D.clear();

  for(std::map< TString, TH2D* >::iterator mapit = maphist_TH2D.begin(); mapit!=maphist_TH2D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH2D.clear();

  for(std::map< TString, TH3D* >::iterator mapit = maphist_TH3D.begin(); mapit!=maphist_TH3D.end(); mapit++){
    delete mapit->second;
  }
  maphist_TH3D.clear();

  //==== output rootfile
  
  if(outfile) outfile->Close();
  delete outfile;

  //==== Tools
  delete MCCorr;
  delete G4Xsec;
  //delete map_BB;

  if (!fChain) return;
  delete fChain->GetCurrentFile();
  cout << "[AnalyzerCore::~AnalyzerCore] JOB FINISHED " << printcurrunttime() << endl;

}
//==================
//==== Read Trees
//==================
Int_t AnalyzerCore::GetEntry(Long64_t entry)
{
  // Read contents of entry 
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

void AnalyzerCore::Loop(){

  Long64_t nentries = fChain->GetEntries();
  if(MaxEvent>0){
    nentries = std::min(nentries,MaxEvent);
  }

  cout << "[AnalyzerCore::Loop] MaxEvent = " << MaxEvent << endl;
  cout << "[AnalyzerCore::Loop] NSkipEvent = " << NSkipEvent << endl;
  cout << "[AnalyzerCore::Loop] LogEvery = " << LogEvery << endl;
  cout << "[AnalyzerCore::Loop] MCSample = " << MCSample << endl;
  cout << "[AnalyzerCore::Loop] Beam_Momentum = " << Beam_Momentum << endl;
  cout << "[AnalyzerCore::Loop] Userflags = {" << endl;
  for(unsigned int i=0; i<Userflags.size(); i++){
    cout << "[AnalyzerCore::Loop]   \"" << Userflags.at(i) << "\"," << endl;
  }
  cout << "[AnalyzerCore::Loop] }" << endl;

  cout << "[AnalyzerCore::Loop] Event Loop Started " << printcurrunttime() << endl;

  for(Long64_t jentry=0; jentry<nentries;jentry++){

    if(jentry<NSkipEvent){
      continue;
    }

    if(jentry%LogEvery==0){
      cout << "[AnalyzerCore::Loop RUNNING] " << jentry << "/" << nentries << " ("<<100.*jentry/nentries<<" %) @ " << printcurrunttime() << endl;
    }

    fChain->GetEntry(jentry);
    Init_evt();
    executeEvent();
  }

  cout << "[AnalyzerCore::Loop] LOOP END " << printcurrunttime() << endl;

}


//==== Attach the historams to ai different direcotry, not outfile
//==== We will write these histograms in WriteHist() to outfile
void AnalyzerCore::SwitchToTempDir(){

  gROOT->cd();
  TDirectory *tempDir = NULL;
  int counter = 0;
  while (!tempDir) {
    //==== First, let's find a directory name that doesn't exist yet
    std::stringstream dirname;
    dirname << "AnalyzerCore" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    //==== Let's try to make this directory
    tempDir = gROOT->mkdir((dirname.str()).c_str());
  }
  tempDir->cd();

}

void AnalyzerCore::SetOutfilePath(TString outname){
  outfile = new TFile(outname,"RECREATE");
};


Event AnalyzerCore::GetEvent(){

  Event ev;
  ev.SetBeam_Momentum(Beam_Momentum);

  return ev;

}

//==================
// Get Daughters
//==================
std::vector<TrueDaughter> AnalyzerCore::GetAllTrueDaughters(){

  vector<TrueDaughter> out;
  for (size_t i = 0; i < evt.true_beam_daughter_PDG->size(); i++){
    TrueDaughter this_Daughter;
    this_Daughter.Set_IsEmpty(false);
    this_Daughter.Set_PDG((*evt.true_beam_daughter_PDG).at(i));
    this_Daughter.Set_ID((*evt.true_beam_daughter_ID).at(i));
    this_Daughter.Set_startX((*evt.true_beam_daughter_startX).at(i));
    this_Daughter.Set_startY((*evt.true_beam_daughter_startY).at(i));
    this_Daughter.Set_startZ((*evt.true_beam_daughter_startZ).at(i));
    this_Daughter.Set_startPx((*evt.true_beam_daughter_startPx).at(i));
    this_Daughter.Set_startPy((*evt.true_beam_daughter_startPy).at(i));
    this_Daughter.Set_startPz((*evt.true_beam_daughter_startPz).at(i));
    this_Daughter.Set_startP((*evt.true_beam_daughter_startP).at(i));
    this_Daughter.Set_endX((*evt.true_beam_daughter_endX).at(i));
    this_Daughter.Set_endY((*evt.true_beam_daughter_endY).at(i));
    this_Daughter.Set_endZ((*evt.true_beam_daughter_endZ).at(i));
    this_Daughter.Set_Process((*evt.true_beam_daughter_Process).at(i));
    this_Daughter.Set_endProcess((*evt.true_beam_daughter_endProcess).at(i));
    out.push_back(this_Daughter);
  }

  return out;
}

std::vector<TrueDaughter> AnalyzerCore::GetPionTrueDaughters(const vector<TrueDaughter>& in){
  vector<TrueDaughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    TrueDaughter this_in = in.at(i);
    if(abs(this_in.PDG()) == 211) out.push_back(this_in);
  }

  return out;
}

std::vector<TrueDaughter> AnalyzerCore::GetProtonTrueDaughters(const vector<TrueDaughter>& in){
  vector<TrueDaughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    TrueDaughter this_in = in.at(i);
    if(this_in.PDG() == 2212) out.push_back(this_in);
  }

  return out;
}

std::vector<TrueDaughter> AnalyzerCore::GetNeutronTrueDaughters(const vector<TrueDaughter>& in){
  vector<TrueDaughter> out;
  for(unsigned int i = 0; i < in.size(); i++){
    TrueDaughter this_in = in.at(i);
    if(this_in.PDG() == 2112) out.push_back(this_in);
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetAllDaughters(){

  // ==== Beam direction using allTrack info : no beam_allTrack trees in AltSCE data...
  TVector3 reco_beam_end(evt.reco_beam_endX, evt.reco_beam_endY, evt.reco_beam_endZ);
  /*
  TVector3 reco_unit_beam(evt.reco_beam_allTrack_trackDirX, evt.reco_beam_allTrack_trackDirY, evt.reco_beam_allTrack_trackDirZ);
  reco_unit_beam = (1. / reco_unit_beam.Mag() ) * reco_unit_beam;
  */

  // ==== Beam direction using calo info
  TVector3 pt0(evt.reco_beam_calo_startX,
	       evt.reco_beam_calo_startY,
	       evt.reco_beam_calo_startZ);
  TVector3 pt1(evt.reco_beam_calo_endX,
	       evt.reco_beam_calo_endY,
	       evt.reco_beam_calo_endZ);
  TVector3 dir = pt1 - pt0;
  TVector3 reco_unit_beam = (1. / dir.Mag() ) * dir;

  vector<Daughter> out;
  for (size_t i = 0; i < evt.reco_daughter_allTrack_ID->size(); i++){
    if((*evt.reco_daughter_allTrack_dQdX_SCE)[i].empty()) continue;
    //if((*evt.reco_daughter_dQdX_SCE)[i].empty()) continue; // == for new MC from Jake

    Daughter this_Daughter;
    //cout << "[PionXsec::GetAllDaughters] i : " << i << endl;
    this_Daughter.SetIsEmpty(false);
    if(evt.MC){
      this_Daughter.Set_PFP_true_byHits_PDG((*evt.reco_daughter_PFP_true_byHits_PDG).at(i));
      this_Daughter.Set_PFP_true_byHits_ID((*evt.reco_daughter_PFP_true_byHits_ID).at(i));
      this_Daughter.Set_PFP_true_byHits_origin((*evt.reco_daughter_PFP_true_byHits_origin).at(i));
      this_Daughter.Set_PFP_true_byHits_parID((*evt.reco_daughter_PFP_true_byHits_parID).at(i));
      this_Daughter.Set_PFP_true_byHits_parPDG((*evt.reco_daughter_PFP_true_byHits_parPDG).at(i));
      this_Daughter.Set_PFP_true_byHits_process((*evt.reco_daughter_PFP_true_byHits_process).at(i));
      this_Daughter.Set_PFP_true_byHits_sharedHits((*evt.reco_daughter_PFP_true_byHits_sharedHits).at(i));
      this_Daughter.Set_PFP_true_byHits_emHits((*evt.reco_daughter_PFP_true_byHits_emHits).at(i));
      this_Daughter.Set_PFP_true_byHits_len((*evt.reco_daughter_PFP_true_byHits_len).at(i));
      this_Daughter.Set_PFP_true_byHits_startX((*evt.reco_daughter_PFP_true_byHits_startX).at(i));
      this_Daughter.Set_PFP_true_byHits_startY((*evt.reco_daughter_PFP_true_byHits_startY).at(i));
      this_Daughter.Set_PFP_true_byHits_startZ((*evt.reco_daughter_PFP_true_byHits_startZ).at(i));
      this_Daughter.Set_PFP_true_byHits_endX((*evt.reco_daughter_PFP_true_byHits_endX).at(i));
      this_Daughter.Set_PFP_true_byHits_endY((*evt.reco_daughter_PFP_true_byHits_endY).at(i));
      this_Daughter.Set_PFP_true_byHits_endZ((*evt.reco_daughter_PFP_true_byHits_endZ).at(i));
      this_Daughter.Set_PFP_true_byHits_startPx((*evt.reco_daughter_PFP_true_byHits_startPx).at(i));
      this_Daughter.Set_PFP_true_byHits_startPy((*evt.reco_daughter_PFP_true_byHits_startPy).at(i));
      this_Daughter.Set_PFP_true_byHits_startPz((*evt.reco_daughter_PFP_true_byHits_startPz).at(i));
      this_Daughter.Set_PFP_true_byHits_startP((*evt.reco_daughter_PFP_true_byHits_startP).at(i));
      this_Daughter.Set_PFP_true_byHits_startE((*evt.reco_daughter_PFP_true_byHits_startE).at(i));
      this_Daughter.Set_PFP_true_byHits_endProcess((*evt.reco_daughter_PFP_true_byHits_endProcess).at(i));
      this_Daughter.Set_PFP_true_byHits_purity((*evt.reco_daughter_PFP_true_byHits_purity).at(i));
      this_Daughter.Set_PFP_true_byHits_completeness((*evt.reco_daughter_PFP_true_byHits_completeness).at(i));
      this_Daughter.Set_PFP_true_byE_PDG((*evt.reco_daughter_PFP_true_byE_PDG).at(i));
      this_Daughter.Set_PFP_true_byE_len((*evt.reco_daughter_PFP_true_byE_len).at(i));
    }

    this_Daughter.Set_PFP_ID((*evt.reco_daughter_PFP_ID).at(i));
    this_Daughter.Set_PFP_nHits((*evt.reco_daughter_PFP_nHits).at(i));
    this_Daughter.Set_PFP_nHits_collection((*evt.reco_daughter_PFP_nHits_collection).at(i));
    this_Daughter.Set_PFP_trackScore((*evt.reco_daughter_PFP_trackScore).at(i));
    this_Daughter.Set_PFP_emScore((*evt.reco_daughter_PFP_emScore).at(i));
    this_Daughter.Set_PFP_michelScore((*evt.reco_daughter_PFP_michelScore).at(i));
    this_Daughter.Set_PFP_trackScore_collection((*evt.reco_daughter_PFP_trackScore_collection).at(i));
    this_Daughter.Set_PFP_emScore_collection((*evt.reco_daughter_PFP_emScore_collection).at(i));
    this_Daughter.Set_PFP_michelScore_collection((*evt.reco_daughter_PFP_michelScore_collection).at(i));
    //this_Daughter.Set_spacePts_X((*evt.reco_daughter_spacePts_X).at(i)); // == Not exist in data
    //this_Daughter.Set_spacePts_Y((*evt.reco_daughter_spacePts_Y).at(i));
    //this_Daughter.Set_spacePts_Z((*evt.reco_daughter_spacePts_Z).at(i));
    this_Daughter.Set_allTrack_ID((*evt.reco_daughter_allTrack_ID).at(i));
    this_Daughter.Set_allTrack_EField_SCE((*evt.reco_daughter_allTrack_EField_SCE).at(i));
    this_Daughter.Set_allTrack_calo_X((*evt.reco_daughter_allTrack_calo_X).at(i));
    this_Daughter.Set_allTrack_calo_Y((*evt.reco_daughter_allTrack_calo_Y).at(i));
    this_Daughter.Set_allTrack_calo_Z((*evt.reco_daughter_allTrack_calo_Z).at(i));
    this_Daughter.Set_allTrack_resRange_SCE((*evt.reco_daughter_allTrack_resRange_SCE).at(i));
    this_Daughter.Set_allTrack_resRange_SCE_plane0((*evt.reco_daughter_allTrack_resRange_plane0).at(i)); // == FIXME, to SCE
    this_Daughter.Set_allTrack_resRange_SCE_plane1((*evt.reco_daughter_allTrack_resRange_plane1).at(i)); // == FIXME, to SCE

    // == commented out for Jake's new 0.5 GeV/c MC
    /*
    this_Daughter.Set_allTrack_calibrated_dEdX_SCE((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE).at(i));
    this_Daughter.Set_allTrack_calibrated_dEdX_SCE_plane0((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE_plane0).at(i));
    this_Daughter.Set_allTrack_calibrated_dEdX_SCE_plane1((*evt.reco_daughter_allTrack_calibrated_dEdX_SCE_plane1).at(i));
    */

    this_Daughter.Set_allTrack_Chi2_proton((*evt.reco_daughter_allTrack_Chi2_proton).at(i));
    this_Daughter.Set_allTrack_Chi2_pion((*evt.reco_daughter_allTrack_Chi2_pion).at(i));
    this_Daughter.Set_allTrack_Chi2_muon((*evt.reco_daughter_allTrack_Chi2_muon).at(i));
    this_Daughter.Set_allTrack_Chi2_ndof((*evt.reco_daughter_allTrack_Chi2_ndof).at(i));
    this_Daughter.Set_allTrack_Chi2_ndof_pion((*evt.reco_daughter_allTrack_Chi2_ndof_pion).at(i));
    this_Daughter.Set_allTrack_Chi2_ndof_muon((*evt.reco_daughter_allTrack_Chi2_ndof_muon).at(i));
    this_Daughter.Set_allTrack_Theta((*evt.reco_daughter_allTrack_Theta).at(i));
    this_Daughter.Set_allTrack_Phi((*evt.reco_daughter_allTrack_Phi).at(i));
    // == FIXME
    //this_Daughter.Set_allTrack_startDirX((*evt.reco_daughter_allTrack_startDirX).at(i));
    //this_Daughter.Set_allTrack_startDirY((*evt.reco_daughter_allTrack_startDirY).at(i));
    //this_Daughter.Set_allTrack_startDirZ((*evt.reco_daughter_allTrack_startDirZ).at(i));
    this_Daughter.Set_allTrack_alt_len((*evt.reco_daughter_allTrack_alt_len).at(i));
    this_Daughter.Set_allTrack_startX((*evt.reco_daughter_allTrack_startX).at(i));
    this_Daughter.Set_allTrack_startY((*evt.reco_daughter_allTrack_startY).at(i));
    this_Daughter.Set_allTrack_startZ((*evt.reco_daughter_allTrack_startZ).at(i));
    this_Daughter.Set_allTrack_endX((*evt.reco_daughter_allTrack_endX).at(i));
    this_Daughter.Set_allTrack_endY((*evt.reco_daughter_allTrack_endY).at(i));
    this_Daughter.Set_allTrack_endZ((*evt.reco_daughter_allTrack_endZ).at(i));
    this_Daughter.Set_allTrack_vertex_michel_score((*evt.reco_daughter_allTrack_vertex_michel_score).at(i));
    this_Daughter.Set_allTrack_vertex_nHits((*evt.reco_daughter_allTrack_vertex_nHits).at(i));

    this_Daughter.Set_allShower_ID((*evt.reco_daughter_allShower_ID).at(i));
    this_Daughter.Set_allShower_len((*evt.reco_daughter_allShower_len).at(i));
    this_Daughter.Set_allShower_startX((*evt.reco_daughter_allShower_startX).at(i));
    this_Daughter.Set_allShower_startY((*evt.reco_daughter_allShower_startY).at(i));
    this_Daughter.Set_allShower_startZ((*evt.reco_daughter_allShower_startZ).at(i));
    this_Daughter.Set_allShower_dirX((*evt.reco_daughter_allShower_dirX).at(i));
    this_Daughter.Set_allShower_dirY((*evt.reco_daughter_allShower_dirY).at(i));
    this_Daughter.Set_allShower_dirZ((*evt.reco_daughter_allShower_dirZ).at(i));
    this_Daughter.Set_allShower_energy((*evt.reco_daughter_allShower_energy).at(i));

    this_Daughter.Set_pandora_type((*evt.reco_daughter_pandora_type).at(i));

    TVector3 unit_daughter((*evt.reco_daughter_allTrack_endX).at(i) - (*evt.reco_daughter_allTrack_startX).at(i),
                           (*evt.reco_daughter_allTrack_endY).at(i) - (*evt.reco_daughter_allTrack_startY).at(i),
                           (*evt.reco_daughter_allTrack_endZ).at(i) - (*evt.reco_daughter_allTrack_startZ).at(i) );
    unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;
    double cos_theta = cos(unit_daughter.Angle(reco_unit_beam));
    this_Daughter.Set_Beam_Cos(cos_theta);

    TVector3 reco_daughter_start((*evt.reco_daughter_allTrack_startX).at(i), (*evt.reco_daughter_allTrack_startY).at(i), (*evt.reco_daughter_allTrack_startZ).at(i) );
    double dist_beam_end = (reco_daughter_start - reco_beam_end).Mag();
    this_Daughter.Set_Beam_Dist(dist_beam_end);

    TVector3 reco_daughter_start_allShower((*evt.reco_daughter_allShower_startX).at(i), (*evt.reco_daughter_allShower_startY).at(i), (*evt.reco_daughter_allShower_startZ).at(i) );
    double dist_beam_end_allShower = (reco_daughter_start_allShower - reco_beam_end).Mag();
    this_Daughter.Set_Beam_Dist_allShower(dist_beam_end_allShower);
    
    out.push_back(this_Daughter);
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetDaughters(const vector<Daughter>& in, int cut_Nhit, double cut_beam_dist, double cut_trackScore){

  vector<Daughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Dist() < cut_beam_dist){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetPions(const vector<Daughter>& in){

  vector<Daughter> out;

  double cut_cos_beam = 0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 60.;
  double cut_startZ = 215.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 > cut_chi2_proton
       && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist && this_in.allTrack_startZ() < cut_startZ){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetStoppingPions(const vector<Daughter>& in, double chi2_pion_cut){

  vector<Daughter> out;

  double cut_cos_beam = 0.95;
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_emScore = 0.5;
  double cut_chi2_proton = 60.;
  double cut_startZ = 215.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() / this_in.allTrack_Chi2_ndof();
    double this_chi2_pion = this_in.allTrack_Chi2_pion() / this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_in.PFP_emScore() < cut_emScore && this_chi2 > cut_chi2_proton
       && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist && this_in.allTrack_startZ() < cut_startZ && this_chi2_pion < chi2_pion_cut){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetProtons(const vector<Daughter>& in){

  vector<Daughter> out;

  double cut_cos_beam = 0.9962; // For protons with scattering angle > 5 deg, TMath::ACos(0.9962) *180. / TMath::Pi() : (double) 4.9965133
  double cut_beam_dist = 10.;
  double cut_trackScore = 0.5;
  double cut_chi2_proton = 10.;
  int cut_Nhit = 20;
  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    double this_chi2 = this_in.allTrack_Chi2_proton() /this_in.allTrack_Chi2_ndof();
    if(this_in.PFP_trackScore() > cut_trackScore && this_chi2 < cut_chi2_proton && this_in.PFP_nHits() > cut_Nhit && this_in.Beam_Cos() < cut_cos_beam && this_in.Beam_Dist() < cut_beam_dist){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetTruePions(const vector<Daughter>& in){

  vector<Daughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(abs(this_true_PID) == 211){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetTrueProtons(const vector<Daughter>& in){

  vector<Daughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(this_true_PID == 2212){
      out.push_back(this_in);
    }
  }

  return out;
}

std::vector<Daughter> AnalyzerCore::GetTrueNeutrons(const vector<Daughter>& in){

  vector<Daughter> out;

  for(unsigned int i = 0; i < in.size(); i++){
    Daughter this_in = in.at(i);
    int this_true_PID = this_in.PFP_true_byHits_PDG();
    if(this_true_PID == 2112){
      out.push_back(this_in);
    }
  }

  return out;
}

//==================
// Set Beam Variables
//==================
double AnalyzerCore::Get_true_ffKE(){
 
  double KE_ff_true = -9999.;

  if(abs(evt.true_beam_PDG) != 13 && abs(evt.true_beam_PDG) != 211 & abs(evt.true_beam_PDG) != 2212) return KE_ff_true;
  //cout << "===========[AnalyzerCore::Get_true_ffKE]" << endl;

  vector<double> true_trklen_accum;
  true_trklen_accum.reserve((*evt.true_beam_traj_Z).size());
  //cout << "true_trklen_accum.size() : " << true_trklen_accum.size() << endl;
  int start_idx = -1;
  for(int i_true_hit = 0; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++){
    //cout << "for i_true_hit : " << i_true_hit << "/" << (*evt.true_beam_traj_Z).size() << endl;
    if((*evt.true_beam_traj_Z).at(i_true_hit) >= 0){
      start_idx = i_true_hit - 1;
      if (start_idx < 0) start_idx = -1;
      break;
    }
    true_trklen_accum[i_true_hit] = 0.;
  }
  if (start_idx >= 0){
    //cout << "if (start_idx >= 0)" << endl;
    double true_trklen = -9999.;
    for (int i_true_hit = start_idx + 1; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++){
      if (i_true_hit == start_idx+1) {
	true_trklen = sqrt( pow( (*evt.true_beam_traj_X).at(i_true_hit)-(*evt.true_beam_traj_X).at(i_true_hit-1), 2)
			    + pow( (*evt.true_beam_traj_Y).at(i_true_hit)-(*evt.true_beam_traj_Y).at(i_true_hit-1), 2)
			    + pow( (*evt.true_beam_traj_Z).at(i_true_hit)-(*evt.true_beam_traj_Z).at(i_true_hit-1), 2)
			    ) * (*evt.true_beam_traj_Z).at(i_true_hit)/((*evt.true_beam_traj_Z).at(i_true_hit)-(*evt.true_beam_traj_Z).at(i_true_hit-1));
      }
      else{
	true_trklen += sqrt( pow( (*evt.true_beam_traj_X).at(i_true_hit)-(*evt.true_beam_traj_X).at(i_true_hit-1), 2)
			     + pow( (*evt.true_beam_traj_Y).at(i_true_hit)-(*evt.true_beam_traj_Y).at(i_true_hit-1), 2)
			     + pow( (*evt.true_beam_traj_Z).at(i_true_hit)-(*evt.true_beam_traj_Z).at(i_true_hit-1), 2)
			     );
      }
      true_trklen_accum[i_true_hit] = true_trklen;
    }
    //cout << "for loop ended" << endl;
    double this_dEdx = map_BB[abs(evt.true_beam_PDG)] -> meandEdx((*evt.true_beam_traj_KE)[start_idx+1]);
    KE_ff_true = (*evt.true_beam_traj_KE).at(start_idx+1) + this_dEdx * true_trklen_accum[start_idx+1];
  }

  return KE_ff_true;
}

double AnalyzerCore::Get_true_beamlen(){

  double beamlen_true = -9999.;
  if(abs(evt.true_beam_PDG) != 13 && abs(evt.true_beam_PDG) != 211 & abs(evt.true_beam_PDG) != 2212) return beamlen_true;
  int start_idx = -1;
  for(int i_true_hit = 0; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++){
    if((*evt.true_beam_traj_Z).at(i_true_hit) >= 0){
      start_idx = i_true_hit - 1;
      if (start_idx < 0) start_idx = -1;
      break;
    }
  }
  if (start_idx >= 0){
    beamlen_true = 0.;
    for (int i_true_hit = start_idx + 1; i_true_hit < (*evt.true_beam_traj_Z).size(); i_true_hit++){
      if (i_true_hit == start_idx+1) {
        beamlen_true = sqrt( pow( (*evt.true_beam_traj_X).at(i_true_hit)-(*evt.true_beam_traj_X).at(i_true_hit-1), 2)
                            + pow( (*evt.true_beam_traj_Y).at(i_true_hit)-(*evt.true_beam_traj_Y).at(i_true_hit-1), 2)
                            + pow( (*evt.true_beam_traj_Z).at(i_true_hit)-(*evt.true_beam_traj_Z).at(i_true_hit-1), 2)
                            ) * (*evt.true_beam_traj_Z).at(i_true_hit)/((*evt.true_beam_traj_Z).at(i_true_hit)-(*evt.true_beam_traj_Z).at(i_true_hit-1));
      }
      else{
        beamlen_true += sqrt( pow( (*evt.true_beam_traj_X).at(i_true_hit)-(*evt.true_beam_traj_X).at(i_true_hit-1), 2)
                             + pow( (*evt.true_beam_traj_Y).at(i_true_hit)-(*evt.true_beam_traj_Y).at(i_true_hit-1), 2)
                             + pow( (*evt.true_beam_traj_Z).at(i_true_hit)-(*evt.true_beam_traj_Z).at(i_true_hit-1), 2)
                             );
      }
    }
  }

  return beamlen_true;
}

void AnalyzerCore::SetPandoraSlicePDG(int pdg){
  pandora_slice_pdg = pdg;
};

int AnalyzerCore::GetPiParType(){

  if (!evt.MC){
    return pi::kData;
  }
  else if (evt.event%2){ // divide half of MC as fake data
    return pi::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){ // the true beam track is not selected
    if (evt.reco_beam_true_byE_origin == 2) {
      return pi::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){ // the selected track is a pion (but not true beam pion, so it is a secondary pion)
      return pi::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return pi::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return pi::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return pi::kMIDeg;
    }
    else {
      return pi::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == -13){
    return pi::kMuon;
  }
  else if (evt.true_beam_PDG == 211){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      return pi::kPiInel;
    }
    else return pi::kPiElas;
  }

  return pi::kMIDother;
}

int AnalyzerCore::GetPi2ParType(){

  if (!evt.MC){
    return pi2::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){ // the true beam track is not selected
    if (evt.reco_beam_true_byE_origin == 2) {
      return pi2::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){ // the selected track is a pion (but not true beam pion, so it is a secondary pion)
      return pi2::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return pi2::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return pi2::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return pi2::kMIDeg;
    }
    else {
      return pi2::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == -13){
    return pi2::kMuon;
  }
  else if (evt.true_beam_PDG == 211){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      int npip = 0;
      int npim = 0;
      int npi0 = 0;
      int np = 0;
      int nn = 0;
      int nmup = 0;
      for (size_t i = 0; i<evt.true_beam_daughter_PDG->size(); ++i){
        int pdg = evt.true_beam_daughter_PDG->at(i);
        if (pdg == 211){
          ++npip;
        }
        else if (pdg == -211){
          ++npim;
        }
        else if (pdg == 111){
          ++npi0;
        }
        else if (pdg == 2212){
          ++np;
        }
        else if (pdg == 2112){
          ++nn;
        }
        else if (pdg == -13 || pdg == 13){
          ++nmup;
        }
      }
      int npis = npip + npim + npi0;
      if (npis > 1){
	return pi2::kPiRes;
      }
      else if (npip){
        return pi2::kPiQE;
      }
      else if (npim){
	return pi2::kPiDCEX;
      }
      else if (npi0){
        return pi2::kPiCEX;
      }
      else if (nmup){
        return pi2::kPiElas;
      }
      else{
        return pi2::kPiABS;
      }
    }
    else{      
      return pi2::kPiElas;
    }
  }

  return pi2::kMIDother;
}

int AnalyzerCore::GetPiTrueType(){

  if (!evt.MC){
    return pitrue::kOther;
  }
  else if (evt.true_beam_PDG == 211 && evt.true_beam_endZ > 0){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      int npip = 0;
      int npim = 0;
      int npi0 = 0;
      int np = 0;
      int nn = 0;
      int nmup = 0;
      for (size_t i = 0; i<evt.true_beam_daughter_PDG->size(); ++i){
        int pdg = evt.true_beam_daughter_PDG->at(i);
        if (pdg == 211){
          ++npip;
        }
        else if (pdg == -211){
          ++npim;
        }
        else if (pdg == 111){
          ++npi0;
        }
        else if (pdg == 2212){
          ++np;
        }
        else if (pdg == 2112){
          ++nn;
        }
        else if (pdg == -13 || pdg == 13){
          ++nmup;
        }
      }
      if (npip || npim){
        return pitrue::kQE;
      }
      else if (npi0){
        return pitrue::kCEX;
      }
      else if (nmup){
        return pitrue::kElas;
      }
      else{
        return pitrue::kABS;
      }
    }
    else{      
      return pitrue::kElas;
    }
  }

  return pitrue::kOther;
}

int AnalyzerCore::GetPParType(){

  if (!evt.MC){
    return p::kData;
  }
  else if (evt.event%2){
    return p::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){
    if (evt.reco_beam_true_byE_origin == 2) {
      return p::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){
      return p::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return p::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return p::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return p::kMIDeg;
    }
    else {
      return p::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == 2212){
    if ((*evt.true_beam_endProcess) == "protonInelastic"){
      return p::kPInel;
    }
    else return p::kPElas;
  }

  return p::kMIDother;
}

//==================
// Event Selections
//==================
bool AnalyzerCore::Pass_Beam_PID(int PID){

  // == Study only proton and pion/muon beams
  if(PID != 2212 && PID != 211 && abs(PID) != 13) return false;

  vector<int> PDGID_candiate_vec;
  PDGID_candiate_vec.clear();
  if(PID == 2212){
    PDGID_candiate_vec.push_back(2212);
  }

  if(PID == 211 || abs(PID) == 13){
    PDGID_candiate_vec.push_back(211);
    PDGID_candiate_vec.push_back(13);
    PDGID_candiate_vec.push_back(-13);
  }

  if (evt.reco_reconstructable_beam_event == 0) return false; // First, remove empty events
  if (evt.MC){
    for (size_t i = 0; i < PDGID_candiate_vec.size(); ++i){
      if (evt.true_beam_PDG == PDGID_candiate_vec[i]) return true; // Truth matched
    }
    return false;
  }
  else{ // data
    if (evt.beam_inst_trigger == 8) return false; // Is cosmics
    if (evt.beam_inst_nMomenta != 1 || evt.beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<PDGID_candiate_vec.size(); ++i){
      for (size_t j = 0; j<evt.beam_inst_PDG_candidates->size(); ++j){
        if ((*evt.beam_inst_PDG_candidates)[j] == PDGID_candiate_vec[i]) return true; // Matching data PID result and PDGID candidates
      }
    }
    return false;
  }
}

bool AnalyzerCore::PassBeamScraperCut() const{

  // == Define beam plug circle in [cm] unit
  double center_x = 0.;
  double center_y = 0.;
  double radius = 0.;
  double N_sigma = 0.;
  if(evt.MC){
    center_x = -29.6;
    center_y = 422.;
    radius = 4.8;
    N_sigma = 1.4;
  }
  else{
    center_x = -32.16;
    center_y = 422.7;
    radius = 4.8;
    N_sigma = 1.2;
  }

  double this_distance = sqrt( pow(evt.beam_inst_X - center_x , 2.) + pow(evt.beam_inst_Y - center_y, 2.) );
  bool out = this_distance < radius * N_sigma;

  return out;
}

bool AnalyzerCore::PassBeamMomentumWindowCut() const{

  bool out = false;
  double P_beam_inst = evt.beam_inst_P * 1000. * P_beam_inst_scale;
  
  //cout << "[AnalyzerCore::PassBeamMomentumWindowCut] P_beam_inst : " << P_beam_inst << ", P_beam_inst_scale : " << P_beam_inst_scale << endl;
  out = ((P_beam_inst > beam_momentum_low) && (P_beam_inst < beam_momentum_high));
  return out;
}

bool AnalyzerCore::PassPandoraSliceCut() const{
  bool out = false;
  out = (evt.reco_beam_type == pandora_slice_pdg);
  return out;
}

bool AnalyzerCore::PassCaloSizeCut() const{ // Require hits information in the collection plane
  bool out = false;
  out = !(evt.reco_beam_calo_wire->empty());
  return out;
}

bool AnalyzerCore::PassAPA3Cut(double cut){ // only use track in the first TPC
  bool out = false;
  out = evt.reco_beam_calo_endZ < cut;
  return out;
}

bool AnalyzerCore::PassMichelScoreCut(double cut){ // further veto muon tracks according to Michel score
  bool out = false;
  out = daughter_michel_score < cut;
  return out;
}

bool AnalyzerCore::PassBeamCosCut(double cut){
  bool out = false;
  out = beam_costh > cut;
  return out;
}

bool AnalyzerCore::PassBeamStartZCut(double cut){
  bool out = false;
  out = evt.reco_beam_calo_startZ > cut;
  return out;
}

bool AnalyzerCore::PassProtonVetoCut(double cut){ // to remove proton background
  bool out = false;
  out = chi2_proton > cut;
  return out;
}

bool AnalyzerCore::PassMuonVetoCut(double cut){ // to remove muon background
  bool out = false;
  out = chi2_muon > cut;
  return out;
}

bool AnalyzerCore::PassStoppedPionVetoCut(double cut){ // to remove stopped pion background
  bool out = false;
  out = chi2_pion > cut;
  return out;
}

//==================
// Initialize
//==================
void AnalyzerCore::initializeAnalyzerTools(){

  
  MCCorr->SetIsData(IsData);
  MCCorr->ReadHistograms();
  
  G4Xsec->ReadHistograms();
}

void AnalyzerCore::Init(){

  cout << "Let initiallize!" << endl;
  evt.Init_PDSPTree(fChain);

  IsData = MCSample.Contains("Data");

  // == Additional Root files
  TString datapath = getenv("DATA_DIR");
  TString datapath_xrootd = "/Users/sungbino/OneDrive/OneDrive/ProtoDUNE-SP/PionKI/PDSPTreeAnalyzer/data/v1/";
  TString is_dunegpvm_str = getenv("PDSPAna_isdunegpvm");
  if(is_dunegpvm_str == "TRUE"){
    datapath_xrootd = "xroot://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune/persistent/users/sungbino/PDSP_data/";
  }
  cout << "[AnalyzerCore::Init] Open : " << datapath_xrootd << endl;
  TFile *file_profile = TFile::Open(datapath_xrootd + "/dEdx_profiles/dEdxrestemplates.root");
  map_profile[13] = (TProfile *)file_profile -> Get("dedx_range_mu");
  map_profile[211] = (TProfile *)file_profile -> Get("dedx_range_pi");
  map_profile[321] = (TProfile *)file_profile -> Get("dedx_range_ka");
  map_profile[2212] = (TProfile *)file_profile -> Get("dedx_range_pro");
  cout << "[[AnalyzerCore::Init]] Called Profiles" << endl;

  // == Beam Window cut
  if(!IsData) P_beam_inst_scale = 1.; // FIXME: Jake's additional MC sample has various weights
  beam_momentum_low = Beam_Momentum * 1000. * 0.8;
  beam_momentum_high = Beam_Momentum * 1000. * 1.2;
  cout << "[[AnalyzerCore::Init]] Called beam window cuts ["  << beam_momentum_low << ", " << beam_momentum_high << "]" << endl;

  // == Pandora
  SetPandoraSlicePDG(13);
  cout << "[[AnalyzerCore::Init]] Set Pandora sllice PDG" << endl;

}

void AnalyzerCore::Init_evt(){
  daughter_michel_score = -999.;
  beam_costh = -999;
  beam_TPC_theta = -999.;
  beam_TPC_phi = -999.;
  delta_X_spec_TPC = -999.;
  delta_Y_spec_TPC = -999.;
  cos_delta_spec_TPC = -999.;
  chi2_proton = -1.;
  if (!evt.reco_beam_calo_wire->empty()){
    if (evt.reco_beam_vertex_nHits) daughter_michel_score = evt.reco_beam_vertex_michel_score_weight_by_charge;

    TVector3 pt0(evt.reco_beam_calo_startX,
                 evt.reco_beam_calo_startY,
                 evt.reco_beam_calo_startZ);
    TVector3 pt1(evt.reco_beam_calo_endX,
                 evt.reco_beam_calo_endY,
                 evt.reco_beam_calo_endZ);
    TVector3 dir = pt1 - pt0;
    dir = dir.Unit();
    if (!IsData){
      TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180),
                       cos(beam_angleY_mc*TMath::Pi()/180),
                       cos(beam_angleZ_mc*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }

    else{
      TVector3 beamdir(cos(beam_angleX_data*TMath::Pi()/180),
                       cos(beam_angleY_data*TMath::Pi()/180),
                       cos(beam_angleZ_data*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }

    beam_TPC_theta = dir.Theta();
    beam_TPC_phi = dir.Phi();
    delta_X_spec_TPC = evt.beam_inst_X - evt.reco_beam_calo_startX;
    delta_Y_spec_TPC = evt.beam_inst_Y - evt.reco_beam_calo_startY;

    TVector3 spec_dir(evt.beam_inst_dirX, evt.beam_inst_dirY, evt.beam_inst_dirZ);
    spec_dir = spec_dir.Unit();
    cos_delta_spec_TPC = dir.Dot(spec_dir);

    chi2_proton = evt.reco_beam_Chi2_proton/evt.reco_beam_Chi2_ndof;
    //chi2_pion = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 211, 1.);
    //chi2_muon = Particle_chi2( (*evt.reco_beam_calibrated_dEdX_SCE), (*evt.reco_beam_resRange_SCE), 13, 1.);
  }
}

//==================
// Additional Functions
//==================
double AnalyzerCore::TOF_to_P(double TOF, int PID){

  double P_measured = -999.;
  double mass = -999.;
  if(PID == 2212) mass = M_proton;
  else if(abs(PID) == 13) mass = M_mu;
  else if(abs(PID) == 211) mass = M_pion;
  else return P_measured;

  //cout << "[AnalyzerCore::TOF_to_P] " << endl;
  double L_over_ct = 1.e09 * L_beamline / (v_light * TOF);
  P_measured = L_over_ct * mass / (sqrt(1. - pow(L_over_ct, 2.)));

  return P_measured; // in MeV/c unit
}


double AnalyzerCore::Convert_P_Spectrometer_to_P_ff(double P_beam_inst, TString particle, TString key, int syst){

  double out = P_beam_inst;
  double delta_P = 0.;
  double p0 = 0., p1 = 0., p2 = 0.;

  if(particle.Contains("pion")){
    if(key == "AllTrue"){
      if(syst == 0){
        p0 = 182.7;
        p1 = -0.4893;
        p2 = 0.0003186;
      }
      else if(syst == -1){
        p0 = 140.2;
        p1 = -0.3979;
        p2 = 0.0002678;
      }
      else if(syst == 1){
        p0 = 235.2;
        p1 = -0.6022;
        p2 = 0.0003807;
      }
      else{
        return out;
      }
    }
  }
  else{
    return out;
  }

  delta_P = p0 + p1 * out + p2 * out * out;

  return out - delta_P;
}

double AnalyzerCore::Particle_chi2(const vector<double> & dEdx, const vector<double> & ResRange, int PID, double dEdx_res_frac){

  //cout << "[AnalyzerCore::Particle_chi2] Start" << endl; 
  if(PID != 2212 && PID != 13 && PID != 211){
    return 99999.;
  }
  if( dEdx.size() < 1 || ResRange.size() < 1 ) return 88888.;

  int N_max_hits = 1000;
  int this_N_calo = dEdx.size();
  int this_N_hits = min(N_max_hits, this_N_calo);
  int N_skip = 1;
  double dEdx_truncate_upper = 1000.;
  double dEdx_truncate_bellow = 0.;
  double this_chi2 = 0.;
  int npt = 0;
  for(int j = N_skip; j < this_N_hits - N_skip; j++){

    double dEdx_measured = dEdx.at(j);
    if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue;

    double this_res_length = ResRange.at(j);
    int bin = map_profile[PID] -> FindBin( this_res_length );
    if( bin >= 1 && bin <= map_profile[PID]->GetNbinsX() ){
      double template_dedx = map_profile[PID]->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( map_profile[PID]->GetBinContent( bin - 1 ) + map_profile[PID]->GetBinContent( bin + 1 ) ) / 2.;
      }

      double template_dedx_err = map_profile[PID]->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( map_profile[PID]->GetBinError( bin - 1 ) + map_profile[PID]->GetBinError( bin + 1 ) ) / 2.;
      }

      double dedx_res = 0.04231 + 0.0001783 * dEdx_measured * dEdx_measured;
      dedx_res *= dEdx_measured * dEdx_res_frac;

      this_chi2 += ( pow( (dEdx_measured - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }
  if(npt == 0) return 77777.;

  return this_chi2 / npt;
}

double AnalyzerCore::Particle_chi2_skip(const vector<double> & dEdx, const vector<double> & ResRange, int PID, double dEdx_res_frac){

  //cout << "[AnalyzerCore::Particle_chi2] Start" << endl;                                                                                                                                                                                     
  if(PID != 2212 && PID != 13 && PID != 211){
    return 99999.;
  }
  int N_skip = 4;

  if( dEdx.size() < N_skip + 1 || ResRange.size() < N_skip + 1 ) return 88888.;

  int N_max_hits = 1000;
  int this_N_calo = dEdx.size();
  int this_N_hits = min(N_max_hits, this_N_calo);
  double dEdx_truncate_upper = 1000.;
  double dEdx_truncate_bellow = 0.;
  double this_chi2 = 0.;
  int npt = 0;
  for(int j = N_skip; j < this_N_hits - N_skip; j++){

    double dEdx_measured = dEdx.at(j);
    if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue;

    double this_res_length = ResRange.at(j);
    int bin = map_profile[PID] -> FindBin( this_res_length );
    if( bin >= 1 && bin <= map_profile[PID]->GetNbinsX() ){
      double template_dedx = map_profile[PID]->GetBinContent( bin );
      if( template_dedx < 1.e-6 ){
        template_dedx = ( map_profile[PID]->GetBinContent( bin - 1 ) + map_profile[PID]->GetBinContent( bin + 1 ) ) / 2.;
      }

      double template_dedx_err = map_profile[PID]->GetBinError( bin );
      if( template_dedx_err < 1.e-6 ){
        template_dedx_err = ( map_profile[PID]->GetBinError( bin - 1 ) + map_profile[PID]->GetBinError( bin + 1 ) ) / 2.;
      }

      double dedx_res = 0.04231 + 0.0001783 * dEdx_measured * dEdx_measured;
      dedx_res *= dEdx_measured * dEdx_res_frac;

      this_chi2 += ( pow( (dEdx_measured - template_dedx), 2 ) / ( pow(template_dedx_err, 2) + pow(dedx_res, 2) ) );

      ++npt;
    }
  }
  if(npt == 0) return 77777.;

  return this_chi2 / npt;
}

double AnalyzerCore::Fit_HypTrkLength_Gaussian(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam){

  //cout << "[HadAna::Fit_dEdx_Residual_Length] Start" << endl;
  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }

  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  int N_skip = 3;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 10){
    //cout << "[HadAna::Fit_dEdx_Residual_Length] Too small number of hits!" << endl;
    return -9999.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  vector<double> chi2_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits
      double this_res_length = ResRange.at(j) + this_additional_res_length;

      double this_KE = map_BB[PID]->KEFromRangeSpline(this_res_length);
      //double dEdx_theory = dEdx_Bethe_Bloch(this_KE, this_mass);
      double dEdx_theory = map_BB[PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(save_graph){
      // == Save vectors for graphes
      chi2_vector.push_back(this_chi2);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  // == Save graph
  if(save_graph){
    // == Vectors for graphes
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = N_skip; i < this_N_hits - N_skip; i++){
      double this_range_original = ResRange.at(i);
      range_original.push_back(this_range_original);

      double this_range_bestfit = ResRange.at(i) + best_additional_res_length;
      range_bestfit.push_back(this_range_bestfit);
      range_reco.push_back(ResRange.at(i));
      dEdx_ordered.push_back(dEdx.at(i));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits - 1));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &chi2_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    chi2_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
  }

  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  if(this_is_beam) original_res_length = ResRange.at(0); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    best_chi2 = 99999.;
    return -9999.;
  }
  else if(best_chi2 > 99990.){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    best_chi2 = 99999.;
    return -9999.;
  }
  else if(best_chi2 < 1.0e-11){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    best_chi2 = 99999.;
    return -9999.;
  }

  return best_total_res_length;
}

double AnalyzerCore::Fit_HypTrkLength_Likelihood(const vector<double> & dEdx, const vector<double> & ResRange, int PID, bool save_graph, bool this_is_beam){
  // == only for charged pions and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }

  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  int N_skip = 3;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow =dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length =max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialization
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;
  int this_N_calo = dEdx.size();
  if(this_N_calo <= 10){
    return -9999.; // == Too small number of hits
  }
  int this_N_hits = this_N_calo; // == Use how many hits
  int i_bestfit = -1;

  // == Fit
  vector<double> m2lnL_vector;
  vector<double> additional_res_legnth_vector;
  for(int i = 0; i < res_length_trial; i++){

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_BB[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      
      // == Likelihood
      double this_pitch = fabs(ResRange.at(j - 1) - ResRange.at(j + 1)) / 2.0;
      double this_likelihood = map_BB[PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      this_m2lnL += (-2.0) * log(this_likelihood);
    }
    //this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(save_graph){
      // == Save vectors for graphes
      m2lnL_vector.push_back(this_m2lnL);
      additional_res_legnth_vector.push_back(this_additional_res_length);
    }
  }

  // == Save graph
  if(save_graph){
    vector<double> range_original;
    vector<double> range_bestfit;
    vector<double> range_reco;
    vector<double> dEdx_ordered;
    for(int i = N_skip; i < this_N_hits - N_skip; i++){
      double this_range_original = ResRange.at(i);
      range_original.push_back(this_range_original);

      double this_range_bestfit = ResRange.at(i) + best_additional_res_length;
      range_bestfit.push_back(this_range_bestfit);
      range_reco.push_back(ResRange.at(i));
      dEdx_ordered.push_back(dEdx.at(i));
    }
    TGraph *dEdx_gr = new TGraph(this_N_hits - 10, &range_original[0], &dEdx_ordered[0]);
    dEdx_gr -> SetName(Form("dEdx_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits - 1));
    dEdx_gr -> Write();
    delete dEdx_gr;

    TGraph *dEdx_bestfit_gr = new TGraph(this_N_hits - 10,&range_bestfit[0], &dEdx_ordered[0]);
    dEdx_bestfit_gr -> SetName(Form("dEdx_bestfit_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_bestfit_gr -> Write();
    delete dEdx_bestfit_gr;

    TGraph *dEdx_reco_gr = new TGraph(this_N_hits - 10,&range_reco[0], &dEdx_ordered[0]);
    dEdx_reco_gr -> SetName(Form("dEdx_reco_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    dEdx_reco_gr -> Write();
    delete dEdx_reco_gr;

    TGraph *chi2_gr = new TGraph(additional_res_legnth_vector.size(), &additional_res_legnth_vector[0], &m2lnL_vector[0]);
    chi2_gr -> SetName(Form("Chi2_Run%d_Evt%d_Nhit%d", evt.run, evt.event, this_N_hits));
    chi2_gr -> Write();
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    delete chi2_gr;
  }

  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  if(this_is_beam) original_res_length = ResRange.at(0); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]

  //cout << "this_N_hits : " << this_N_hits << ", original_res_length : " << original_res_length << ", best_additional_res_length : " << best_additional_res_length << ", best_total_res_length : " << best_total_res_length << endl;

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : no mimumum" << endl;
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    return -9999.;
  }
  else if(best_m2lnL > 99990.){
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 > 99990." << endl;
    return -9999.;
  }
  else if(best_m2lnL < 1.0e-11){
    m2lnL_vector.clear();
    additional_res_legnth_vector.clear();
    //cout << "[HadAna::Fit_Beam_Hit_dEdx_Residual_Length] Fit failed : best_chi2 < 1.0e-11" << endl;
    return -9999.;
  }

  // == Return
  m2lnL_vector.clear();
  additional_res_legnth_vector.clear();
  return best_total_res_length;
}

double AnalyzerCore::KE_Hypfit_Gaussian(const Daughter in, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  int N_skip = 3;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_chi2 = 99999.;
  vector<double> dEdx = in.allTrack_calibrated_dEdX_SCE();
  vector<double> ResRange = in.allTrack_resRange_SCE();

  int this_N_calo = dEdx.size();
  if(this_N_calo <= 10){
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  // == Fit
  for(int i = 0; i < res_length_trial; i++){
    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_chi2 = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip hits
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_BB[PID]->KEFromRangeSpline(this_res_length);
      //double dEdx_theory = dEdx_Bethe_Bloch(this_KE, this_mass);
      double dEdx_theory = map_BB[PID]->meandEdx(this_KE);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate
      // == Gaussian approx.
      //double dEdx_theory_err = dEdx_theory * 0.02;
      this_chi2 += pow(dEdx_measured - dEdx_theory, 2);
    }
    this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_chi2 < best_chi2){
      best_chi2 = this_chi2;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
  }

  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length;

  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_chi2 > 99990.){ // == Fit failed : best_chi2 > 99990."
    return -6666.;
  }
  else if(best_chi2 < 1.0e-11){ // == Fit failed : best_chi2 < 1.0e-11
    return -5555.;
  }

  double best_KE = map_BB[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}

double AnalyzerCore::KE_Hypfit_Likelihood(const Daughter in, int PID){

  // == PID input : mass hypothesis, valid only for muons, charged pions, and protons
  if(!(PID == 13 || PID == 2212 || PID == 211)){
    return -9999.;
  }
  // == Tunable parameters
  double min_additional_res_length = 0.;
  double max_additional_res_length = max_additional_res_length_pion;
  double res_length_step = res_length_step_pion;
  int N_skip = 3;
  double dEdx_truncate_upper = dEdx_truncate_upper_pion;
  double dEdx_truncate_bellow = dEdx_truncate_bellow_pion;
  if(PID == 2212){
    max_additional_res_length = max_additional_res_length_proton;
    dEdx_truncate_upper = dEdx_truncate_upper_proton;
    dEdx_truncate_bellow = dEdx_truncate_bellow_proton;
    res_length_step = res_length_step_proton;
  }
  int res_length_trial = (max_additional_res_length - min_additional_res_length) / res_length_step;

  // == Initialize
  double best_additional_res_length = -0.1;
  double best_m2lnL = 99999.;
  vector<double> dEdx = in.allTrack_calibrated_dEdX_SCE();
  vector<double> ResRange = in.allTrack_resRange_SCE();

  int this_N_calo = dEdx.size();
  if(this_N_calo <= 10){
    return -8888.; // == Too small number of hits
  }
  int i_bestfit = -1;
  int this_N_hits = this_N_calo;

  double default_m2lnL = -1.;
  // == Fit
  for(int i = 0; i < res_length_trial; i++){

    double this_additional_res_length = min_additional_res_length + (i + 0.) * res_length_step;
    double this_m2lnL = 0.;
    for(int j = N_skip; j < this_N_hits - N_skip; j++){ // == Do not use first and last N_skip this
      double this_res_length = ResRange.at(j) + this_additional_res_length;
      double this_KE = map_BB[PID]->KEFromRangeSpline(this_res_length);
      double dEdx_measured = dEdx.at(j);
      if(dEdx_measured < dEdx_truncate_bellow || dEdx_measured > dEdx_truncate_upper) continue; // == Truncate

      // == Likelihood
      double this_pitch = fabs(ResRange.at(j - 1) - ResRange.at(j + 1)) / 2.0;
      double this_likelihood = map_BB[PID] -> dEdx_PDF(this_KE, this_pitch, dEdx_measured);
      if(this_likelihood > 1e-6)this_m2lnL += (-2.0) * log(this_likelihood);
    }
    //this_chi2 = this_chi2 / (this_N_hits + 0.); // == chi2 / n.d.f
    if(this_m2lnL < best_m2lnL){
      best_m2lnL = this_m2lnL;
      best_additional_res_length = this_additional_res_length;
      i_bestfit = i;
    }
    if(i == 0){
      default_m2lnL = this_m2lnL;
    }
  }

  //cout << "default_m2lnL : " << default_m2lnL << ", best_m2lnL : " << best_m2lnL << endl;
  // == Result
  double original_res_length = ResRange.at(this_N_calo - 1); // == [cm]
  double best_total_res_length = best_additional_res_length + original_res_length; // == [cm]
  // == Define fitting failed cases
  if(i_bestfit == res_length_trial - 1){ // == Fit failed : no mimumum
    return -7777.;
  }
  else if(best_m2lnL > 99990.){ // == Fit failed : best_chi2 > 99990.
    return -6666.;
  }
  else if(best_m2lnL < 1.0e-11){ // == Fit failed : best_chi2 < 1.0e-11
    return -5555.;
  }

  // == Return
  double best_KE = map_BB[PID] -> KEFromRangeSpline(best_total_res_length);
  return best_KE;
}


double AnalyzerCore::Get_EQE_NC_Pion(double P_pion, double cos_theta, double E_binding, int which_sol){

  double E_pion = sqrt( pow(P_pion, 2.0) + pow(M_pion, 2.0) );

  double A = M_proton - E_binding - E_pion;
  double B = pow(M_pion, 2.) - pow(P_pion, 2.) - pow(M_proton, 2.);

  // == ax^2 + bx + c = 0
  double a = 4. * (A*A - P_pion * P_pion * cos_theta * cos_theta);
  double b = 4. * A * (A*A + B);
  double c = pow(A*A + B, 2.) + 4. * M_pion * M_pion * P_pion * P_pion * cos_theta * cos_theta;

  double numer1 = (-1.) * b;
  double numer_sqrt = sqrt(b*b - 4. * a * c);
  double denom = 2. * a;

  double EQE = (numer1 + (which_sol + 0.) * numer_sqrt ) / denom;

  return EQE;
}

double AnalyzerCore::Get_EQE_NC_Proton(double P_proton, double cos_theta, double E_binding, int which_sol){

  double E_proton = sqrt( pow(P_proton, 2.0) + pow(M_proton, 2.0) );
  
  double A = M_proton - E_binding - E_proton;

  // == ax^2 + bx + c = 0 
  double a = 4. * (A*A - P_proton * P_proton * cos_theta * cos_theta);
  double b = 4. * A * (A*A - P_proton * P_proton);
  double c = pow(A*A - P_proton * P_proton, 2.) + 4. * M_pion * M_pion * P_proton * P_proton * cos_theta * cos_theta;

  double numer1 = (-1.) * b;
  double numer_sqrt = sqrt(b*b - 4. * a * c);
  double denom = 2. * a;

  double EQE = (numer1 + (which_sol + 0.) * numer_sqrt ) / denom;

  return EQE;
}

bool AnalyzerCore::Is_EQE(double window){
  bool is_QE = false;

  TVector3 unit_beam(evt.true_beam_endPx, evt.true_beam_endPy, evt.true_beam_endPz);
  double P_beam = evt.true_beam_endP * 1000.;
  double m_beam = evt.true_beam_mass;
  double E_beam = sqrt(P_beam * P_beam + m_beam * m_beam);

  int N_daugther_piplus = 0;
  for(unsigned int i = 0; i < (*evt.true_beam_daughter_ID).size(); i++){
    int this_daughter_PID = (*evt.true_beam_daughter_PDG).at(i);
    if(this_daughter_PID == 211){
      N_daugther_piplus++;

      TVector3 unit_daughter((*evt.true_beam_daughter_startPx).at(i), (*evt.true_beam_daughter_startPy).at(i), (*evt.true_beam_daughter_startPz).at(i) );
      unit_daughter = (1./ unit_daughter.Mag() ) * unit_daughter;
      double cos_theta = cos(unit_beam.Angle(unit_daughter));
      double P_daughter = (*evt.true_beam_daughter_startP).at(i) * 1000.;
      double EQE_NC_4 = Get_EQE_NC_Pion(P_daughter, cos_theta, 4., -1.);
      double this_EQEmE = EQE_NC_4 - E_beam;
      if(this_EQEmE > (-1.) * window && this_EQEmE < window){
        is_QE = true;
      }
    }
  }

  return is_QE;
}

//==================
// MCS
//==================
vector<MCSSegment> AnalyzerCore::SplitIntoSegments(const vector<TVector3> & hits, double segment_size){

  vector<MCSSegment> out;
  if(hits.size() < 6) return out;

  TVector3 hit_start = hits.at(0);
  int hit_start_index = 0;
  for(unsigned int i_hit = 2; i_hit < hits.size(); i_hit++){
    vector<TVector3> this_hit_collection;
    double this_segment_size = 0.;
    for(unsigned int j_hit = hit_start_index; j_hit < i_hit; j_hit++){
      this_hit_collection.push_back(hits.at(j_hit));
      if(j_hit < i_hit - 1){
	double this_pitch = (hits.at(j_hit + 1) - hits.at(j_hit)).Mag();
	this_segment_size = this_segment_size + this_pitch;
      }
    }
    if(this_hit_collection.size() < 3) continue;

    TVector3 reco_start = hits.at(hit_start_index);
    TVector3 reco_end = hits.at(i_hit);

    vector<TVector3> this_fitted = Fitter -> Linear3Dfit(this_hit_collection);
    if(this_fitted.size() == 0) continue;

    if(i_hit == hits.size() - 1){
      //cout << "[AnalyzerCore::SplitIntoSegments] Last hit" << endl;
      MCSSegment this_segment;
      TVector3 fitted_start = this_fitted.at(0);
      TVector3 fitted_end = this_fitted.at(0) +this_fitted.at(1);
      TVector3 fitted_vec = this_fitted.at(1);
      this_segment.SetMCSSegment(reco_start, reco_end, fitted_start, fitted_end, fitted_vec, this_segment_size);

      if(fitted_vec.Mag() > 0.) out.push_back(this_segment);
    }
    else if(this_segment_size > segment_size){
      //cout << "[AnalyzerCore::SplitIntoSegments] New segment with i_hit : " << i_hit << endl;
      MCSSegment this_segment;
      TVector3 fitted_start = this_fitted.at(0);
      TVector3 fitted_end = this_fitted.at(0) + this_fitted.at(1);
      TVector3 fitted_vec = this_fitted.at(1);
      this_segment.SetMCSSegment(reco_start, reco_end, fitted_start, fitted_end, fitted_vec, this_segment_size);

      out.push_back(this_segment);
      hit_start_index = i_hit;
      //cout << "[AnalyzerCore::SplitIntoSegments] hit_start_index : " << hit_start_index << endl;
    }

    this_hit_collection.clear();
    this_fitted.clear();
  }

  return out;
}

TVector3 AnalyzerCore::RotateToZaxis(TVector3 reference_vec, TVector3 input_vec){
  double theta = reference_vec.Theta();
  TVector3 ortho_line_on_xy(1., -1. * reference_vec.X() / reference_vec.Y(), 0.);
  TVector3 cross = reference_vec.Cross(ortho_line_on_xy);
  if(cross.Z() < 0.) theta = -1. * theta;
  // double dot = reference_vec.Dot(ortho_line_on_xy);
  //cout << "[AnalyzerCore::RotateToZaxis] dot : " << dot << ", cross.Z() : " << cross.Z() << endl;
  TVector3 out(input_vec);
  out.Rotate(-1 * theta, ortho_line_on_xy);
  return out;
}

double AnalyzerCore::MCS_Get_HL_Sigma(double segment_size, double P, double mass){
  double kappa = ( (HL_kappa_a / (P*P)) + HL_kappa_c );
  double one_over_pbeta = pow(P*P + mass * mass, 0.5) / (P*P);
  double root_term = pow(segment_size / 14., 0.5);

  double sigma_HL = kappa * one_over_pbeta * root_term * (1 + HL_epsilon * log(segment_size / 14.));
  double out = pow(sigma_HL * sigma_HL + HL_sigma_res * HL_sigma_res, 0.5);
  return out;
}

double AnalyzerCore::MCS_Get_Likelihood(double HL_sigma, double delta_angle){
  double out = TMath::Log(HL_sigma) + 0.5 * pow(delta_angle / HL_sigma, 2);
  return out;
}

double AnalyzerCore::MCS_Likelihood_Fitting(const vector<MCSSegment> segments, double segment_size, int PID){

  cout << "[AnalyzerCore::MCS_Likelihood_Fitting] Start" << endl;

  double out = -9999.;
  if(segments.size() < 3) return out;

  double P_step = 2.; // -- [MeV]
  int N_step = 1000;
  double total_range = 0.;
  for(unsigned int i = 0; i < segments.size(); i++){
    double this_segment_range = segments.at(i).Range();
    total_range = total_range + this_segment_range;
  }

  double KE_from_range = map_BB[abs(PID)] -> KEFromRangeSpline(total_range);
  double P_from_range = map_BB[abs(PID)] -> KEtoMomentum(KE_from_range);
  double best_P = -999.;
  double best_mlogL = 99999999.;
  int i_best_mlogL = -1.;
  for(int i = 0; i < N_step; i++){
    vector<double> this_P_vec;
    vector<double> this_theta_xz_vec;
    vector<double> this_theta_yz_vec;
    double this_P = P_from_range + P_step * (i + 0.);
    this_P_vec.push_back(this_P);
    for(int j = 0; j < segments.size() - 1; j++){
      TVector3 this_vec = segments.at(j).FittedVec();
      TVector3 next_vec = segments.at(j + 1).FittedVec();
      TVector3 rotated_this_vec = RotateToZaxis(this_vec, this_vec);
      TVector3 rotated_next_vec = RotateToZaxis(this_vec, next_vec);
      double this_theta_xz = TMath::ATan(rotated_next_vec.X() / rotated_next_vec.Z());
      double this_theta_yz = TMath::ATan(rotated_next_vec.Y() / rotated_next_vec.Z());
      this_theta_xz_vec.push_back(this_theta_xz);
      this_theta_yz_vec.push_back(this_theta_yz);

      if(j < segments.size() - 2){
	double this_partial_range = segments.at(j).Range();
	double this_KE = map_BB[abs(PID)] -> MomentumtoKE(this_P);
	double next_KE = map_BB[abs(PID)] -> KEAtLength(this_KE, this_partial_range);
	double next_P = map_BB[abs(PID)] -> KEtoMomentum(next_KE);
	this_P_vec.push_back(next_P);
	this_P = next_P;
      }
    }

    for(unsigned int j = 0; j < this_theta_xz_vec.size(); j++){
      cout << Form("%d P : %.2f, theta_xz : %.2e, theta_yz :  %.2e", j, this_P_vec.at(j), this_theta_xz_vec.at(j), this_theta_yz_vec.at(j)) << endl;
    }

    this_P_vec.clear();
    this_theta_xz_vec.clear();
    this_theta_yz_vec.clear();
  }

  return out;
}


//=================
//==== Fittings
//=================
TF1 *AnalyzerCore::langaufit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF, Int_t *Status, TString FunName){

  Int_t i;
  TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
  if (ffitold) delete ffitold;

  TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
  ffit->SetParameters(startvalues);
  ffit->SetParNames("Width","MPV","Area","GSigma");

  for (i=0; i<4; i++) {
    ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
  }

  TFitResultPtr fitres = his->Fit(FunName,"RBOSQN");
  ffit->GetParameters(fitparams);
  for (i=0; i<4; i++) {
    fiterrors[i] = ffit->GetParError(i);
  }

  ChiSqr[0] = ffit->GetChisquare();
  NDF[0] = ffit->GetNDF();
  Status[0] = fitres->CovMatrixStatus();

  return (ffit);
}

//==================
//==== Plotting
//==================
TH1D* AnalyzerCore::GetHist1D(TString histname){

  TH1D *h = NULL;
  std::map<TString, TH1D*>::iterator mapit = maphist_TH1D.find(histname);
  if(mapit != maphist_TH1D.end()) return mapit->second;

  return h;

}

TH2D* AnalyzerCore::GetHist2D(TString histname){

  TH2D *h = NULL;
  std::map<TString, TH2D*>::iterator mapit = maphist_TH2D.find(histname);
  if(mapit != maphist_TH2D.end()) return mapit->second;

  return h;

}

TH3D* AnalyzerCore::GetHist3D(TString histname){

  TH3D *h = NULL;
  std::map<TString, TH3D*>::iterator mapit = maphist_TH3D.find(histname);
  if(mapit != maphist_TH3D.end()) return mapit->second;

  return h;

}

void AnalyzerCore::FillHist(TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

void AnalyzerCore::FillHist(TString histname, double value, double weight, int n_bin, double *xbins){

  TH1D *this_hist = GetHist1D(histname);
  if( !this_hist ){
    this_hist = new TH1D(histname, "", n_bin, xbins);
    this_hist->SetDirectory(NULL);
    maphist_TH1D[histname] = this_hist;
  }

  this_hist->Fill(value, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y,
			    double weight,
			    int n_binx, double x_min, double x_max,
			    int n_biny, double y_min, double y_max){

  TH2D *this_hist = GetHist2D(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    this_hist->SetDirectory(NULL);
    maphist_TH2D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y,
			    double weight,
			    int n_binx, double *xbins,
			    int n_biny, double *ybins){

  TH2D *this_hist = GetHist2D(histname);
  if( !this_hist ){
    this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
    this_hist->SetDirectory(NULL);
    maphist_TH2D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y, double value_z,
			    double weight,
			    int n_binx, double x_min, double x_max,
			    int n_biny, double y_min, double y_max,
			    int n_binz, double z_min, double z_max){

  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    this_hist = new TH3D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max, n_binz, z_min, z_max);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, weight);

}

void AnalyzerCore::FillHist(TString histname,
			    double value_x, double value_y, double value_z,
			    double weight,
			    int n_binx, double *xbins,
			    int n_biny, double *ybins,
			    int n_binz, double *zbins){

  TH3D *this_hist = GetHist3D(histname);
  if( !this_hist ){
    this_hist = new TH3D(histname, "", n_binx, xbins, n_biny, ybins, n_binz, zbins);
    this_hist->SetDirectory(NULL);
    maphist_TH3D[histname] = this_hist;
  }

  this_hist->Fill(value_x, value_y, value_z, weight);

}

TH1D* AnalyzerCore::JSGetHist1D(TString suffix, TString histname){

  TH1D *h = NULL;

  std::map< TString, std::map<TString, TH1D*> >::iterator mapit = JSmaphist_TH1D.find(suffix);
  if(mapit==JSmaphist_TH1D.end()){
    return h;
  }
  else{

    std::map<TString, TH1D*> this_maphist = mapit->second;
    std::map<TString, TH1D*>::iterator mapitit = this_maphist.find(histname);
    if(mapitit != this_maphist.end()) return mapitit->second;

  }

  return h;

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname, double value, double weight, int n_bin, double x_min, double x_max){

  TH1D *this_hist = JSGetHist1D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH1D(histname, "", n_bin, x_min, x_max);
    (JSmaphist_TH1D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value, weight);

}

TH2D* AnalyzerCore::JSGetHist2D(TString suffix, TString histname){

  TH2D *h = NULL;

  std::map< TString, std::map<TString, TH2D*> >::iterator mapit = JSmaphist_TH2D.find(suffix);
  if(mapit==JSmaphist_TH2D.end()){
    return h;
  }
  else{

    std::map<TString, TH2D*> this_maphist = mapit->second;
    std::map<TString, TH2D*>::iterator mapitit = this_maphist.find(histname);
    if(mapitit != this_maphist.end()) return mapitit->second;

  }

  return h;

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname,
			      double value_x, double value_y,
			      double weight,
			      int n_binx, double x_min, double x_max,
			      int n_biny, double y_min, double y_max){

  TH2D *this_hist = JSGetHist2D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH2D(histname, "", n_binx, x_min, x_max, n_biny, y_min, y_max);
    (JSmaphist_TH2D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::JSFillHist(TString suffix, TString histname,
			      double value_x, double value_y,
			      double weight,
			      int n_binx, double *xbins,
			      int n_biny, double *ybins){

  TH2D *this_hist = JSGetHist2D(suffix, histname);
  if( !this_hist ){

    this_hist = new TH2D(histname, "", n_binx, xbins, n_biny, ybins);
    (JSmaphist_TH2D[suffix])[histname] = this_hist;

  }

  this_hist->Fill(value_x, value_y, weight);

}

void AnalyzerCore::WriteHist(){

  outfile->cd();
  for(std::map< TString, TH1D* >::iterator mapit = maphist_TH1D.begin(); mapit!=maphist_TH1D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
  for(std::map< TString, TH2D* >::iterator mapit = maphist_TH2D.begin(); mapit!=maphist_TH2D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }
  for(std::map< TString, TH3D* >::iterator mapit = maphist_TH3D.begin(); mapit!=maphist_TH3D.end(); mapit++){
    TString this_fullname=mapit->second->GetName();
    TString this_name=this_fullname(this_fullname.Last('/')+1,this_fullname.Length());
    TString this_suffix=this_fullname(0,this_fullname.Last('/'));
    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);
    mapit->second->Write(this_name);
    outfile->cd();
  }

  outfile->cd();
  for(std::map< TString, std::map<TString, TH1D*> >::iterator mapit=JSmaphist_TH1D.begin(); mapit!=JSmaphist_TH1D.end(); mapit++){

    TString this_suffix = mapit->first;
    std::map< TString, TH1D* > this_maphist = mapit->second;


    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);

    for(std::map< TString, TH1D* >::iterator mapit = this_maphist.begin(); mapit!=this_maphist.end(); mapit++){
      mapit->second->Write();
    }

    outfile->cd();

  }

  for(std::map< TString, std::map<TString, TH2D*> >::iterator mapit=JSmaphist_TH2D.begin(); mapit!=JSmaphist_TH2D.end(); mapit++){

    TString this_suffix = mapit->first;
    std::map< TString, TH2D* > this_maphist = mapit->second;

    TDirectory *dir = outfile->GetDirectory(this_suffix);
    if(!dir){
      outfile->mkdir(this_suffix);
    }
    outfile->cd(this_suffix);

    for(std::map< TString, TH2D* >::iterator mapit = this_maphist.begin(); mapit!=this_maphist.end(); mapit++){
      mapit->second->Write();
    }

    outfile->cd();

  }

  //outfile->cd();
  if(fEventTree != nullptr) fEventTree->Write();

  outfile->cd();
  outfile->mkdir("Fitter");
  outfile->cd("Fitter");
  for(std::map<TString, TGraph2D*>::iterator mapit = Fitter -> map_TGraph2D.begin(); mapit!=Fitter -> map_TGraph2D.end(); mapit++){
    mapit->second->Write();
  }

}
