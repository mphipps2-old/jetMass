#ifndef jetMass_jetMass_H
#define jetMass_jetMass_H

#include <AnaAlgorithm/AnaAlgorithm.h>
// GRL
#include <AsgAnalysisInterfaces/IGoodRunsListSelectionTool.h>
#include <AsgTools/AnaToolHandle.h>
#include <TrigConfInterfaces/ITrigConfigTool.h>
#include <TrigDecisionTool/TrigDecisionTool.h>
#include <MuonAnalysisInterfaces/IMuonSelectionTool.h>
#include <MuonAnalysisInterfaces/IMuonCalibrationAndSmearingTool.h>
#include "InDetTrackSelectionTool/InDetTrackSelectionTool.h"
#include <xAODJet/JetContainer.h>
#include <JetCalibTools/IJetCalibrationTool.h>
#include <TriggerMatchingTool/IMatchingTool.h>
#include "xAODTruth/TruthParticleContainer.h"

// Others
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <utility>

class jetMass : public EL::AnaAlgorithm
{
public:
  // this is a standard algorithm constructor
  jetMass (const std::string& name, ISvcLocator* pSvcLocator);
  
  // these are the functions inherited from Algorithm
  virtual StatusCode initialize () override;
  virtual StatusCode execute () override;
  virtual StatusCode finalize () override;

  void  clearVector();
  void  EventPlaneAnalysis();
  double dR(const double eta1,
            const double phi1,
            const double eta2,
            const double phi2);
  int    parent_classify(const xAOD::TruthParticle *theParticle);
  double mindRmuon(const xAOD::TrackParticle* trk);
  double mindRTrk(const xAOD::TrackParticle* metrk);
  int TruthMatching(double jetEta, double jetPhi,std::unique_ptr<xAOD::JetContainer> &truthJets);
private:
 
 
  asg::AnaToolHandle<IGoodRunsListSelectionTool> m_grl; //!
  asg::AnaToolHandle<Trig::TrigDecisionTool> m_trigDecisionTool; //!
  asg::AnaToolHandle<TrigConf::ITrigConfigTool> m_trigConfigTool; //!
  // MuonSelectionTool
  asg::AnaToolHandle<CP::IMuonSelectionTool> m_muonSelection; //!
  // MuonCalibrationAndSmearing
  asg::AnaToolHandle<CP::IMuonCalibrationAndSmearingTool> m_muonCalibrationAndSmearingTool; //!
  asg::AnaToolHandle<InDet::IInDetTrackSelectionTool> m_trkSelection; //!
  asg::AnaToolHandle<Trig::IMatchingTool> m_tmt; //!

  std::string m_outputName;
  std::string m_outputHist;
  float m_mcSumWight = 0;
  TF1 *m_d0_cut; //!
  TTree *m_tree; //!
  TTree *m_tree2; //!
  TTree *m_tree3; //!
  int m_EventNumber; //!
  int m_RunNumber; //!
  int m_LumiBlock; //!
  int m_mcEventNumber; //!
  int m_mcChannelNumber; //!
  float m_mcEventWeight; //!
  double m_FCalWeight; //!
  double m_jetWeight; //!
  double m_jzNorm; //!
  double m_FCal_Et; //!
  double m_primaryVertexZ; //!

  std::string m_grl_fileName; 
  std::string m_reco_jet_collection; 
  std::string m_truth_jet_collection; 
  std::string m_track_pt_cut;
  double m_jet_dR_truth_matching; 
  double m_jet_eta_cut; 
  double m_jet_pT_cut; 
  double m_truth_track_prob_cut; 
  double m_truth_jet_pT_cut; 
  double m_track_eta_cut; 
  bool m_isPP; 
  bool m_isMC; 

  std::vector<double>          m_Jet_pt; //!
  std::vector<double>          m_Jet_e; //!
  std::vector<double>          m_Jet_eta; //!
  std::vector<double>          m_Jet_y; //!
  std::vector<double>          m_Jet_phi; //!

  std::vector<double>          m_TruthJet_pt; //!
  std::vector<double>          m_TruthJet_e; //!
  std::vector<double>          m_TruthJet_eta; //!
  std::vector<double>          m_TruthJet_y; //!
  std::vector<double>          m_TruthJet_phi; //!
  std::vector<int>             m_TruthJet_pid; //!

  std::vector<double> m_Jet_sv1_pu; //!
  std::vector<double> m_Jet_sv1_pb; //!
  std::vector<double> m_Jet_sv1_pc; //!
  std::vector<float>  m_Jet_sv1_efrc; //!
  std::vector<float>  m_Jet_sv1_mass; //!
  std::vector<int>    m_Jet_sv1_n2t; //!
  std::vector<int>    m_Jet_sv1_ntrkv; //!
  std::vector<float>  m_Jet_sv1_LXY; //!
  std::vector<float>  m_Jet_sv1_L3d; //!
  std::vector<float>  m_Jet_sv1_sig3d; //!
  std::vector<float>  m_Jet_sv1_distmatlay; //!
  std::vector<float>  m_Jet_sv1_dR; //!

  std::vector<std::vector<double> > m_JetTrack_dR; //!
  std::vector<std::vector<double> > m_JetTrack_d0; //!
  std::vector<std::vector<double> > m_JetTrack_d0_cut; //!
  std::vector<std::vector<int> >    m_JetTrack_charge; //!
  std::vector<std::vector<double> > m_JetTrack_pt; //!
  std::vector<std::vector<double> > m_JetTrack_pz; //!
  std::vector<std::vector<double> > m_JetTrack_e; //!
  std::vector<std::vector<double> > m_JetTrack_eta; //!
  std::vector<std::vector<double> > m_JetTrack_phi; //!
  std::vector<std::vector<int> >    m_JetTrack_pid; //!
  std::vector<std::vector<int> >    m_JetTrack_status; //!
  std::vector<std::vector<int> >    m_JetTrack_type; //!
  std::vector<std::vector<double> > m_JetTrack_truth_dR; //!
  std::vector<std::vector<bool> >   m_JetTrack_truthMatched; //!

  std::vector<std::vector<int> >    m_JetTrackTruth_charge; //!
  std::vector<std::vector<double> > m_JetTrackTruth_pt; //!
  std::vector<std::vector<double> > m_JetTrackTruth_e; //!
  std::vector<std::vector<double> > m_JetTrackTruth_eta; //! 
  std::vector<std::vector<double> > m_JetTrackTruth_phi; //!

  std::vector<std::string> m_TriggerObject_Chain; //!
  std::vector<float>   m_TriggerObject_Ps; //!  
  // Trigger Tag and probe members
  std::vector<double> m_Muon_pt; //!
  std::vector<double> m_Muon_eta; //!
  std::vector<double> m_Muon_phi; //!
  std::vector<int>    m_Muon_charge; //!
  std::vector<double> m_Muon_parent; //!
  std::vector<double> m_Muon_e; //!
  std::vector<double> m_Muon_quality; //!
//
  // Configuration, and any other types of variables go here.
  //float m_cutValue;
  //TTree *m_myTree;
  //TH1 *m_myHist;
  };
  
 #endif
