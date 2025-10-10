#ifndef Daughter_h
#define Daughter_h

#include <string>
#include <vector>

//#include <TROOT.h>


using namespace std;

class Daughter{
public:

  Daughter();
  virtual ~Daughter();

  // == Empty
  void SetIsEmpty(bool i_IsEmpty);
  inline bool IsEmpty() const { return j_IsEmpty; }

  // == Variables
  void Set_PFP_true_byHits_PDG(int i_PFP_true_byHits_PDG);
  void Set_PFP_true_byHits_ID(int i_PFP_true_byHits_ID);
  void Set_PFP_true_byHits_origin(int i_PFP_true_byHits_origin);
  void Set_PFP_true_byHits_parID(int i_PFP_true_byHits_parID);
  void Set_PFP_true_byHits_parPDG(int i_PFP_true_byHits_parPDG);
  void Set_PFP_true_byHits_process(string i_PFP_true_byHits_process);
  void Set_PFP_true_byHits_sharedHits(unsigned long i_PFP_true_byHits_sharedHits);
  void Set_PFP_true_byHits_emHits(unsigned long i_PFP_true_byHits_emHits);
  void Set_PFP_true_byHits_len(double i_PFP_true_byHits_len);
  void Set_PFP_true_byHits_startX(double i_PFP_true_byHits_startX);
  void Set_PFP_true_byHits_startY(double i_PFP_true_byHits_startY);
  void Set_PFP_true_byHits_startZ(double i_PFP_true_byHits_startZ);
  void Set_PFP_true_byHits_endX(double i_PFP_true_byHits_endX);
  void Set_PFP_true_byHits_endY(double i_PFP_true_byHits_endY);
  void Set_PFP_true_byHits_endZ(double i_PFP_true_byHits_endZ);
  void Set_PFP_true_byHits_startPx(double i_PFP_true_byHits_startPx);
  void Set_PFP_true_byHits_startPy(double i_PFP_true_byHits_startPy);
  void Set_PFP_true_byHits_startPz(double i_PFP_true_byHits_startPz);
  void Set_PFP_true_byHits_startP(double i_PFP_true_byHits_startP);
  void Set_PFP_true_byHits_startE(double i_PFP_true_byHits_startE);
  void Set_PFP_true_byHits_endProcess(string i_PFP_true_byHits_endProcess);
  void Set_PFP_true_byHits_purity(double i_PFP_true_byHits_purity);
  void Set_PFP_true_byHits_completeness(double i_PFP_true_byHits_completeness);
  void Set_PFP_true_byE_PDG(int i_PFP_true_byE_PDG);
  void Set_PFP_true_byE_len(double i_PFP_true_byE_len);
  void Set_PFP_ID(int i_PFP_ID);
  void Set_PFP_nHits(int i_PFP_nHits);
  void Set_PFP_nHits_collection(int i_PFP_nHits_collection);
  void Set_PFP_trackScore(double i_PFP_trackScore);
  void Set_PFP_emScore(double i_PFP_emScore);
  void Set_PFP_michelScore(double i_PFP_michelScore);
  void Set_PFP_trackScore_collection(double i_PFP_trackScore_collection);
  void Set_PFP_emScore_collection(double i_PFP_emScore_collection);
  void Set_PFP_michelScore_collection(double i_PFP_michelScore_collection);
  void Set_spacePts_X(vector<double> i_spacePts_X);
  void Set_spacePts_Y(vector<double> i_spacePts_Y);
  void Set_spacePts_Z(vector<double> i_spacePts_Z);
  void Set_allTrack_ID(int i_allTrack_ID);
  void Set_allTrack_EField_SCE(vector<double> i_allTrack_EField_SCE);
  void Set_allTrack_calo_X(vector<double> i_allTrack_calo_X);
  void Set_allTrack_calo_Y(vector<double> i_allTrack_calo_Y);
  void Set_allTrack_calo_Z(vector<double> i_allTrack_calo_Z);
  void Set_allTrack_resRange_SCE(vector<double> i_allTrack_resRange_SCE);
  void Set_allTrack_resRange_SCE_plane0(vector<double> i_allTrack_resRange_SCE_plane0);
  void Set_allTrack_resRange_SCE_plane1(vector<double> i_allTrack_resRange_SCE_plane1);
  void Set_allTrack_calibrated_dEdX_SCE(vector<double> i_allTrack_calibrated_dEdX_SCE);
  void Set_allTrack_calibrated_dEdX_SCE_plane0(vector<double> i_allTrack_calibrated_dEdX_SCE_plane0);
  void Set_allTrack_calibrated_dEdX_SCE_plane1(vector<double> i_allTrack_calibrated_dEdX_SCE_plane1);
  void Set_allTrack_Chi2_proton(double i_allTrack_Chi2_proton);
  void Set_allTrack_Chi2_pion(double i_allTrack_Chi2_pion);
  void Set_allTrack_Chi2_muon(double i_allTrack_Chi2_muon);
  void Set_allTrack_Chi2_ndof(double i_allTrack_Chi2_ndof);
  void Set_allTrack_Chi2_ndof_pion(double i_allTrack_Chi2_ndof_pion);
  void Set_allTrack_Chi2_ndof_muon(double i_allTrack_Chi2_ndof_muon);
  void Set_allTrack_Theta(double i_allTrack_Theta);
  void Set_allTrack_Phi(double i_allTrack_Phi);
  void Set_allTrack_startDirX(double i_allTrack_startDirX);
  void Set_allTrack_startDirY(double i_allTrack_startDirY);
  void Set_allTrack_startDirZ(double i_allTrack_startDirZ);
  void Set_allTrack_alt_len(double i_allTrack_alt_len);
  void Set_allTrack_startX(double i_allTrack_startX);
  void Set_allTrack_startY(double i_allTrack_startY);
  void Set_allTrack_startZ(double i_allTrack_startZ);
  void Set_allTrack_endX(double i_allTrack_endX);
  void Set_allTrack_endY(double i_allTrack_endY);
  void Set_allTrack_endZ(double i_allTrack_endZ);
  void Set_allTrack_vertex_michel_score(double i_allTrack_vertex_michel_score);
  void Set_allTrack_vertex_nHits(int i_allTrack_vertex_nHits);
  void Set_allShower_ID(int i_allShower_ID);
  void Set_allShower_len(double i_allShower_len);
  void Set_allShower_startX(double i_allShower_startX);
  void Set_allShower_startY(double i_allShower_startY);
  void Set_allShower_startZ(double i_allShower_startZ);
  void Set_allShower_dirX(double i_allShower_dirX);
  void Set_allShower_dirY(double i_allShower_dirY);
  void Set_allShower_dirZ(double i_allShower_dirZ);
  void Set_allShower_energy(double i_allShower_energy);
  void Set_pandora_type(int i_pandora_type);
  void Set_Beam_Cos(double i_Beam_Cos);
  void Set_Beam_Dist(double i_Beam_Dist);
  void Set_Beam_Dist_allShower(double i_Beam_Dist_allShower); 
  
  inline int PFP_true_byHits_PDG() const { return j_PFP_true_byHits_PDG; }
  inline int PFP_true_byHits_ID() const { return j_PFP_true_byHits_ID; }
  inline int PFP_true_byHits_origin() const { return j_PFP_true_byHits_origin; }
  inline int PFP_true_byHits_parID() const { return j_PFP_true_byHits_parID; }
  inline int PFP_true_byHits_parPDG() const  { return j_PFP_true_byHits_parPDG; }
  inline string PFP_true_byHits_process() const { return j_PFP_true_byHits_process; }
  inline unsigned long PFP_true_byHits_sharedHits() const { return j_PFP_true_byHits_sharedHits; }
  inline unsigned long PFP_true_byHits_emHits() const { return j_PFP_true_byHits_emHits; }
  inline double PFP_true_byHits_len() const { return j_PFP_true_byHits_len; }
  inline double PFP_true_byHits_startX() const { return j_PFP_true_byHits_startX; }
  inline double PFP_true_byHits_startY() const { return j_PFP_true_byHits_startY; }
  inline double PFP_true_byHits_startZ() const { return j_PFP_true_byHits_startZ; }
  inline double PFP_true_byHits_endX() const { return j_PFP_true_byHits_endX; }
  inline double PFP_true_byHits_endY() const { return j_PFP_true_byHits_endY; }
  inline double PFP_true_byHits_endZ() const { return j_PFP_true_byHits_endZ; }
  inline double PFP_true_byHits_startPx() const { return j_PFP_true_byHits_startPx; }
  inline double PFP_true_byHits_startPy() const { return j_PFP_true_byHits_startPy; }
  inline double PFP_true_byHits_startPz() const { return j_PFP_true_byHits_startPz; }
  inline double PFP_true_byHits_startP() const { return j_PFP_true_byHits_startP; }
  inline double PFP_true_byHits_startE() const { return j_PFP_true_byHits_startE; }
  inline string PFP_true_byHits_endProcess() const { return j_PFP_true_byHits_endProcess; }
  inline double PFP_true_byHits_purity() const { return j_PFP_true_byHits_purity; }
  inline double PFP_true_byHits_completeness() const { return j_PFP_true_byHits_completeness; }
  inline int PFP_true_byE_PDG() const { return j_PFP_true_byE_PDG; }
  inline double PFP_true_byE_len() const { return j_PFP_true_byE_len; }
  inline int PFP_ID() const { return j_PFP_ID; }
  inline int PFP_nHits() const { return j_PFP_nHits; }
  inline int PFP_nHits_collection() const { return j_PFP_nHits_collection; }
  inline double PFP_trackScore() const { return j_PFP_trackScore; }
  inline double PFP_emScore() const { return j_PFP_emScore; }
  inline double PFP_michelScore() const { return j_PFP_michelScore; }
  inline double PFP_trackScore_collection() const { return j_PFP_trackScore_collection; }
  inline double PFP_emScore_collection() const { return j_PFP_emScore_collection; }
  inline double PFP_michelScore_collection() const { return j_PFP_michelScore_collection; }
  inline vector<double> spacePts_X() const { return j_spacePts_X; }
  inline vector<double> spacePts_Y() const { return j_spacePts_Y; }
  inline vector<double> spacePts_Z() const { return j_spacePts_Z; }
  inline int allTrack_ID() const { return j_allTrack_ID; }
  inline vector<double> allTrack_EField_SCE() const { return j_allTrack_EField_SCE; }
  inline vector<double> allTrack_calo_X() const { return j_allTrack_calo_X; }
  inline vector<double> allTrack_calo_Y() const{ return j_allTrack_calo_Y; }
  inline vector<double> allTrack_calo_Z() const{ return j_allTrack_calo_Z; }
  inline vector<double> allTrack_resRange_SCE() const { return j_allTrack_resRange_SCE; }
  inline vector<double> allTrack_resRange_SCE_plane0() const { return j_allTrack_resRange_SCE_plane0; }
  inline vector<double> allTrack_resRange_SCE_plane1() const { return j_allTrack_resRange_SCE_plane1; }
  inline vector<double> allTrack_calibrated_dEdX_SCE() const { return j_allTrack_calibrated_dEdX_SCE; }
  inline vector<double> allTrack_calibrated_dEdX_SCE_plane0() const { return j_allTrack_calibrated_dEdX_SCE_plane0; }
  inline vector<double> allTrack_calibrated_dEdX_SCE_plane1() const { return j_allTrack_calibrated_dEdX_SCE_plane1; }
  inline double allTrack_Chi2_proton() const { return j_allTrack_Chi2_proton; }
  inline double allTrack_Chi2_pion() const { return j_allTrack_Chi2_pion; }
  inline double allTrack_Chi2_muon() const { return j_allTrack_Chi2_muon; }
  inline double allTrack_Chi2_ndof() const { return j_allTrack_Chi2_ndof; }
  inline double allTrack_Chi2_ndof_pion() const { return j_allTrack_Chi2_ndof_pion; }
  inline double allTrack_Chi2_ndof_muon() const { return j_allTrack_Chi2_ndof_muon; }
  inline double allTrack_Theta() const { return j_allTrack_Theta; }
  inline double allTrack_Phi() const { return j_allTrack_Phi; }
  inline double allTrack_startDirX() const { return j_allTrack_startDirX; }
  inline double allTrack_startDirY() const { return j_allTrack_startDirY; }
  inline double allTrack_startDirZ() const { return j_allTrack_startDirZ; }
  inline double allTrack_alt_len() const { return j_allTrack_alt_len; }
  inline double allTrack_startX() const { return j_allTrack_startX; }
  inline double allTrack_startY() const { return j_allTrack_startY; }
  inline double allTrack_startZ() const { return j_allTrack_startZ; }
  inline double allTrack_endX() const { return j_allTrack_endX; }
  inline double allTrack_endY() const { return j_allTrack_endY; }
  inline double allTrack_endZ() const { return j_allTrack_endZ; }
  inline double allTrack_vertex_michel_score() const { return j_allTrack_vertex_michel_score; }
  inline int allTrack_vertex_nHits() const { return j_allTrack_vertex_nHits; }
  inline int allShower_ID() const { return j_allShower_ID; }
  inline double allShower_len() const { return j_allShower_len; }
  inline double allShower_startX() const { return j_allShower_startX; }
  inline double	allShower_startY() const { return j_allShower_startY; }
  inline double	allShower_startZ() const { return j_allShower_startZ; }
  inline double	allShower_dirX() const { return j_allShower_dirX; }
  inline double allShower_dirY() const { return j_allShower_dirY; }
  inline double allShower_dirZ() const { return j_allShower_dirZ; }
  inline double allShower_energy() const { return j_allShower_energy; }
  inline int pandora_type() const { return j_pandora_type; }
  inline double Beam_Cos() const { return j_Beam_Cos; }
  inline double Beam_Dist() const { return j_Beam_Dist; }
  inline double Beam_Dist_allShower() const { return j_Beam_Dist_allShower; }

  
private:
  bool j_IsEmpty;

  int j_PFP_true_byHits_PDG;
  int j_PFP_true_byHits_ID;
  int j_PFP_true_byHits_origin;
  int j_PFP_true_byHits_parID;
  int j_PFP_true_byHits_parPDG;
  string j_PFP_true_byHits_process;
  unsigned long j_PFP_true_byHits_sharedHits;
  unsigned long j_PFP_true_byHits_emHits;
  double j_PFP_true_byHits_len;
  double j_PFP_true_byHits_startX;
  double j_PFP_true_byHits_startY;
  double j_PFP_true_byHits_startZ;
  double j_PFP_true_byHits_endX;
  double j_PFP_true_byHits_endY;
  double j_PFP_true_byHits_endZ;
  double j_PFP_true_byHits_startPx;
  double j_PFP_true_byHits_startPy;
  double j_PFP_true_byHits_startPz;
  double j_PFP_true_byHits_startP;
  double j_PFP_true_byHits_startE;
  string j_PFP_true_byHits_endProcess;
  double j_PFP_true_byHits_purity;
  double j_PFP_true_byHits_completeness;
  int j_PFP_true_byE_PDG;
  double j_PFP_true_byE_len;
  int j_PFP_ID;
  int j_PFP_nHits;
  int j_PFP_nHits_collection;
  double j_PFP_trackScore;
  double j_PFP_emScore;
  double j_PFP_michelScore;
  double j_PFP_trackScore_collection;
  double j_PFP_emScore_collection;
  double j_PFP_michelScore_collection;
  vector<double> j_spacePts_X;
  vector<double> j_spacePts_Y;
  vector<double> j_spacePts_Z;
  int j_allTrack_ID;
  vector<double> j_allTrack_EField_SCE;
  vector<double> j_allTrack_calo_X;
  vector<double> j_allTrack_calo_Y;
  vector<double> j_allTrack_calo_Z;
  vector<double> j_allTrack_resRange_SCE;
  vector<double> j_allTrack_resRange_SCE_plane0;
  vector<double> j_allTrack_resRange_SCE_plane1;
  vector<double> j_allTrack_calibrated_dEdX_SCE;
  vector<double> j_allTrack_calibrated_dEdX_SCE_plane0;
  vector<double> j_allTrack_calibrated_dEdX_SCE_plane1;
  double j_allTrack_Chi2_proton;
  double j_allTrack_Chi2_pion;
  double j_allTrack_Chi2_muon;
  double j_allTrack_Chi2_ndof;
  double j_allTrack_Chi2_ndof_pion;
  double j_allTrack_Chi2_ndof_muon;
  double j_allTrack_Theta;
  double j_allTrack_Phi;
  double j_allTrack_startDirX;
  double j_allTrack_startDirY;
  double j_allTrack_startDirZ;
  double j_allTrack_alt_len;
  double j_allTrack_startX;
  double j_allTrack_startY;
  double j_allTrack_startZ;
  double j_allTrack_endX;
  double j_allTrack_endY;
  double j_allTrack_endZ;
  double j_allTrack_vertex_michel_score;
  int j_allTrack_vertex_nHits;
  int j_allShower_ID;
  double j_allShower_len;
  double j_allShower_startX;
  double j_allShower_startY;
  double j_allShower_startZ;
  double j_allShower_dirX;
  double j_allShower_dirY;
  double j_allShower_dirZ;
  double j_allShower_energy;
  int j_pandora_type;
  double j_Beam_Cos;
  double j_Beam_Dist;
  double j_Beam_Dist_allShower;
};

#endif
