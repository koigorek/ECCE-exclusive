// Tell emacs that this is a C++ source
//  -*- C++ -*-.
#ifndef EXCL_ANA_H
#define EXCL_ANA_H

#include <fun4all/SubsysReco.h>

#include <string>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>

#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"


class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TTree;
class TNtuple;
class CaloEvalStack;
class CaloRawClusterEval;
class RawClusterContainer;
class SvtxTrackMap;
class SvtxEvalStack;
class SvtxTrackEval;
class PHG4TruthInfoContainer;


class excl_ana : public SubsysReco
{
 public:
CaloEvalStack *_caloevalstack;
CaloEvalStack* _caloevalstackFEMC;
CaloEvalStack* _caloevalstackEEMC;
CaloEvalStack* _caloevalstackCEMC;

SvtxEvalStack *_svtxEvalStack;
//  excl_ana(const std::string &name = "excl_ana");
  excl_ana(const std::string &name = "Diff_Tagg_ana", const std::string &fname = "MyNtuple.root");

  virtual ~excl_ana();

  /** Called during initialization.
      Typically this is where you can book histograms, and e.g.
      register them to Fun4AllServer (so they can be output to file
      using Fun4AllServer::dumpHistos() method).
   */
  int Init(PHCompositeNode *topNode) override;

  /** Called for first event when run number is known.
      Typically this is where you may want to fetch data from
      database, because you know the run number. A place
      to book histograms which have to know the run number.
   */
  int InitRun(PHCompositeNode *topNode) override;

  /** Called for each event.
      This is where you do the real work.
   */
  int process_event(PHCompositeNode *topNode) override;

  /// Clean up internals after each event.
  int ResetEvent(PHCompositeNode *topNode) override;

  /// Called at the end of each run.
  int EndRun(const int runnumber) override;

  /// Called at the end of all processing.
  int End(PHCompositeNode *topNode) override;

  /// Reset
  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  int process_g4hits(PHCompositeNode *);
  int process_RomanPots(PHCompositeNode *, int);
  int process_B0(PHCompositeNode *);
  int process_ClusterCalo(PHCompositeNode*, std::string);
  int process_tracks(PHCompositeNode *);

 private:

 TTree *tree;

float Epx;
float Epy;
float Epz;

float Ppx;
float Ppy;
float Ppz;

float Gpx;
float Gpy;
float Gpz;

    int RP1;
    int RP2;
    int RPhits;
    
    float RPx[10000];
    float RPy[10000];
    float RPz[10000];
    int RPind[10000];

    Int_t B0hits;
    
    Float_t B0x[10000];
    Float_t B0y[10000];
    Float_t B0z[10000];
    Int_t B0ind[10000];
    
int nHits;
int caloInd[10000];
float hitX[10000];
float hitY[10000];
float hitZ[10000];
float hitE[10000];
int hitPid[10000];
int hitsNtowers[10000];
int hitsEEMC;
int hitsFEMC;
int hitsCEMC;

int ntr;
float tr_px[10000];
float tr_py[10000];
float tr_pz[10000];
float tr_p[10000];
float tr_phi[10000];
float tr_eta[10000];
float charge[10000];
float tr_x[10000];
float tr_y[10000];
float tr_z[10000];
int tr_Pid[10000];

  void initializeVariables();
  void initializeTrees();

 protected:

  std::string detector;
  std::string outfilename;
  Fun4AllHistoManager *hm;

  TFile *outfile;
  TNtuple *g4hitntuple;

  unsigned long long int event_itt;
  gsl_rng* m_RandomGenerator;

  //*********************************
  // Energy and Position smearing

  float EMCAL_Smear(float E);
  float HCAL_Smear(float E);
  float PbWO4_Smear(float E);
  float Position_Smear(float E);

  //---------------------
  // From ejana

  double true_q2;
  double true_x;
  double true_s_e;
  double true_xpi;
  double true_ypi;
  double true_tpi;

  double have_true_dis_info = false;
  
  bool  HIT_IN_ZDC; 
  bool  HIT_IN_HEC;	

  double e_beam_energy;
  double ion_beam_energy;

  double crossing_angle;

  TLorentzVector r_lelectron;
//  TLorentzVector r_lproton;
//  TLorentzVector r_lproton;

  TLorentzVector r_lscatelec;
  TLorentzVector r_l_scat_nucleon;

  TLorentzVector lproton;

  Int_t ZDC_hit;

  TH2F* h2_ZDC_XY; 
  TH2F* h2_ZDC_XY_double; 

  TH1F* h1_E_dep_smeared;
  TH1F* h1_E_dep;

};

#endif // DIFF_TAGG_ANA_H
