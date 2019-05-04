#include <AsgTools/MessageCheck.h>
//#include <EventLoop/Job.h>
//#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <jetMass/jetMass.h>
#include <xAODEventInfo/EventInfo.h>
#include <TSystem.h>
#include <xAODJet/JetContainer.h>
#include <xAODEgamma/PhotonContainer.h>
#include "xAODHIEvent/HIEventShapeContainer.h"
#include <xAODMuon/MuonContainer.h>
#include "xAODTracking/TrackParticleContainer.h"
#include <TFile.h>
#include <PATInterfaces/CorrectionCode.h> // to check the return correction code status of tools
#include <xAODCore/ShallowAuxContainer.h>
#include <xAODCore/ShallowCopy.h>
#include <ParticleJetTools/JetFlavourInfo.h>
#include "xAODTruth/xAODTruthHelpers.h"
#include "xAODTruth/TruthParticleContainer.h"
#include "xAODTruth/TruthEventContainer.h"
#include <PathResolver/PathResolver.h>
#include <utility>
#include <TF1.h>
#include <TRandom3.h>
#include "data/jzWeight.h"
#include "fastjet/JetDefinition.hh"

using namespace fastjet;

jetMass :: jetMass (const std::string& name,
		    ISvcLocator *pSvcLocator)
  : EL::AnaAlgorithm (name, pSvcLocator),
  m_grl ("GoodRunsListSelectionTool/grl", this),
  m_trigDecisionTool ("Trig::TrigDecisionTool/TrigDecisionTool"),	
    m_trigConfigTool("TrigConf::xAODConfigTool/xAODConfigTool"),
  m_tmt ("Trig::MatchingTool/MyMatchingTool", this)
{
  m_grl.declarePropertyFor (this, "grlTool"); 
  declareProperty ("GRL_FileName",m_grl_fileName="jetMass/data17_5TeV.periodAllYear_DetStatus-v98-pro21-16_Unknown_PHYS_StandardGRL_All_Good_25ns_ignore_GLOBAL_LOWMU.xml","GRL file");
  declareProperty("RecoJetContainer",m_reco_jet_collection="AntiKt4HIJets","Reco jet container");
  declareProperty("TruthJetContainer",m_truth_jet_collection="AntiKt4TruthJets","Truth jet container");
  declareProperty("JetEtaCut",m_jet_eta_cut=2.1, "Max jet eta");
  declareProperty("JetPtCut",m_jet_pT_cut=30., "Min jet pt cut");
  declareProperty("Jet_dr_truthMatching",m_jet_dR_truth_matching=0.2, "Truth Matching radius");
  declareProperty("TruthTrackProbCut",m_truth_track_prob_cut=0.3,"MC truth track probability cut");
  declareProperty("TruthJetPtCut",m_truth_jet_pT_cut=30.,"Min truth jet pt cut");
  declareProperty("MC_Flag",m_isMC=true,"MC flag");
  declareProperty("pp_Flag",m_isPP=true,"pp flag");

}


StatusCode jetMass :: initialize ()
{
  ANA_MSG_INFO ("in initialize");

  ANA_CHECK (book (TH1F ("h_RejectionHisto","h_RejectionHisto",5,0,5)));

  TFile *outputFile = wk()->getOutputFile ("myOutput");
  // ANA_CHECK (book (TTree ("tree_bTag","jetMass ntuple")));
  //  TTree *myTree = tree ("tree_bTag");
  m_tree = new TTree ("tree_bTag","bTagging Tree");
  m_tree->SetDirectory (outputFile);

  m_tree->Branch ("EventNumber", &m_EventNumber);
  m_tree->Branch ("RunNumber", &m_RunNumber);
  m_tree->Branch ("LumiBlock", &m_LumiBlock);
  m_tree->Branch ("Fcal", &m_FCal_Et);
  m_tree->Branch ("FcalWeight", &m_FCalWeight);
  m_tree->Branch ("JzNorm", &m_jzNorm);
  m_tree->Branch ("PrimaryVertexZ", &m_primaryVertexZ);

  // Antikt4HI jets
  m_tree->Branch ("Jet_pt", &m_Jet_pt);
  m_tree->Branch ("Jet_pt2", &m_Jet_pt2);
  m_tree->Branch ("Jet_e",  &m_Jet_e);
  m_tree->Branch ("Jet_eta", &m_Jet_eta);
  m_tree->Branch ("Jet_y", &m_Jet_y);
  m_tree->Branch ("Jet_phi", &m_Jet_phi);
  m_tree->Branch ("Jet_mass2", &m_Jet_mass2);
  m_tree->Branch ("Jet_mass2OverPt2", &m_Jet_mass2OverPt2);

  // Antikt4HI truth jets
  m_tree->Branch ("TruthJet_pid", &m_TruthJet_pid);
  m_tree->Branch ("TruthJet_pt", &m_TruthJet_pt);
  m_tree->Branch ("TruthJet_pt2", &m_TruthJet_pt2);
  m_tree->Branch ("TruthJet_e",  &m_TruthJet_e);
  m_tree->Branch ("TruthJet_eta", &m_TruthJet_eta);
  m_tree->Branch ("TruthJet_y", &m_TruthJet_y);
  m_tree->Branch ("TruthJet_phi", &m_TruthJet_phi);
  m_tree->Branch ("TruthJet_mass2", &m_TruthJet_mass2);
  m_tree->Branch ("TruthJet_mass2OverPt2", &m_TruthJet_mass2OverPt2);
  
  // setting GRL 
  std::string GRLFilePath = PathResolverFindCalibFile(m_grl_fileName.c_str());
  std::vector<std::string> vecStringGRL;
  vecStringGRL.push_back(GRLFilePath);
  ANA_CHECK(m_grl.setProperty( "GoodRunsListVec", vecStringGRL));
  ANA_CHECK(m_grl.setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
  ANA_CHECK(m_grl.initialize());
  if (!m_isMC) {
    // Trigger configuration
    ANA_CHECK (m_trigConfigTool.initialize());
    // Trigger decision tool
    ANA_CHECK (m_trigDecisionTool.setProperty ("ConfigTool", m_trigConfigTool.getHandle()));
    ANA_CHECK (m_trigDecisionTool.setProperty ("TrigDecisionKey", "xTrigDecision"));
    ANA_CHECK (m_trigDecisionTool.initialize());
    // trigger matching tool
    ANA_CHECK (m_tmt.setProperty("TrigDecisionTool",m_trigDecisionTool.getHandle()));
    ANA_CHECK (m_tmt.initialize());
  }

  return StatusCode::SUCCESS;
}

StatusCode jetMass :: execute ()
{
  // retrieve the eventInfo object from the event store
  //  ANA_MSG_INFO ("execute");
  const xAOD::EventInfo *eventInfo = nullptr;
  ANA_CHECK (evtStore()->retrieve (eventInfo, "EventInfo"));
  m_EventNumber = eventInfo->eventNumber(); 
  m_RunNumber   = eventInfo->runNumber();
  m_LumiBlock   = eventInfo->lumiBlock();
  m_FCal_Et = -99; 
  m_FCalWeight = -1; 
  m_primaryVertexZ = -1;

  bool isMC = false;
  // check if the event is MC
  if (eventInfo->eventType (xAOD::EventInfo::IS_SIMULATION)) {
    isMC = true; // can do something with this later
  }

  if(isMC){ 
    m_mcEventNumber         = eventInfo->mcEventNumber();
    m_mcChannelNumber       = eventInfo->mcChannelNumber();
    m_mcEventWeight         = eventInfo->mcEventWeight();
    //    m_mcSumWeight += m_mcEventWeight;
  }

  if (m_mcChannelNumber == 420010 && m_isPP) m_jzNorm = perEvtWgtJZ0_PythDijet;
  else if (m_mcChannelNumber == 420011 && m_isPP) m_jzNorm = perEvtWgtJZ1_PythDijet;
  else if (m_mcChannelNumber == 420012 && m_isPP) m_jzNorm = perEvtWgtJZ2_PythDijet;
  else if (m_mcChannelNumber == 420013 && m_isPP) m_jzNorm = perEvtWgtJZ3_PythDijet;
  else if (m_mcChannelNumber == 420014 && m_isPP) m_jzNorm = perEvtWgtJZ4_PythDijet;
  else if (m_mcChannelNumber == 420015 && m_isPP) m_jzNorm = perEvtWgtJZ5_PythDijet;
  else {
    m_jzNorm = -1;
    ANA_MSG_INFO("warning: m_mcChannelNumber = " << m_mcChannelNumber << " eventNumber " << m_EventNumber);
  }

  //  ANA_MSG_INFO ("Event Cuts");
  // Event cuts
  hist ("h_RejectionHisto")->Fill(0.5); 
  bool keep = true; 

  // if data check if event passes GRL
  if (!isMC) { // it's data!
    if (!m_grl->passRunLB(*eventInfo)) {
      //ANA_MSG_INFO ("drop event: GRL");
      hist ("h_RejectionHisto")->Fill(1.5);
      keep = false; 
    }
  
    if ((eventInfo->errorState(xAOD::EventInfo::LAr) == xAOD::EventInfo::Error) ||
        (eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error) ||
        (eventInfo->errorState(xAOD::EventInfo::SCT) == xAOD::EventInfo::Error) ||
        (eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18))) {
      hist ("h_RejectionHisto")->Fill(2.5);
      keep = false; 
    } // end if event flags check
  } // end if not MC


  const xAOD::VertexContainer * vertices = 0;
  if ( !evtStore()->retrieve( vertices, "PrimaryVertices" ).isSuccess() ){
    ANA_MSG_ERROR("execute(): Failed to retrieve VertexContainer container. Exiting." );
    keep = false;
  }
  // at least one vertex -- note: first vertex is a dummy vertex (ie size < 2)
  if(vertices->size()<2) {
    hist("h_RejectionHisto")->Fill(3.5);
    keep = false;
  }
  // Find primary vertex
  xAOD::VertexContainer::const_iterator vtx_itr = vertices->begin();
  xAOD::VertexContainer::const_iterator vtx_end = vertices->end();
  const xAOD::Vertex* primaryVertex = 0;
  for(;vtx_itr!=vtx_end;++vtx_itr)	  {	
  if((*vtx_itr)->vertexType()==xAOD::VxType::PriVtx) {
      primaryVertex = (*vtx_itr);
      m_primaryVertexZ = fabs(primaryVertex->z());
      break;
    }
  }
 

  /*
  //Pileup rejection -- should be added for PbPb once tool is ready
  bool doPileupRejection = true;
  if (!m_isPP && doPileupRejection){
    bool m_is_pileup = false;
    if (!isMC) {
      const xAOD::ZdcModuleContainer* zdcMod = 0;
      EL_RETURN_CHECK("execute",evtStore()->retrieve( zdcMod, "ZdcModules"));
      // ZDC
      m_zdcTools->reprocessZdc();
      // is Pileup
      m_is_pileup = m_hiPileup->is_pileup( *calos, *zdcMod); // SAVE pileup Decision HERE 0 = NO pileup, 1 = pileup
    }
    else m_is_pileup = (FCalEt > 4.8); //Remove pileup in MC -- for MC data overlay since min bias can be pileup as well

    if (m_is_pileup){
      h_RejectionHisto->Fill(6.5);
      keep = false;
    }
  }
  */
  if (!keep) return StatusCode::SUCCESS; // go to the next event
  hist("h_RejectionHisto")->Fill(4.5);

  //  ANA_MSG_INFO ("begin muons ");

  if (!m_isPP && !m_isMC) {
    // only to be used for PbPb data
    const xAOD::HIEventShapeContainer *hiue(0);
    ANA_CHECK( evtStore()->retrieve(hiue, "CaloSums"));
    double m_fcalEt = hiue->at(5)->et();
    m_FCal_Et = m_fcalEt*1e-6;
  }
  
  if (!m_isMC) {
    bool dum_HLT_mu8 = false;
    bool dum_HLT_mu10 = false;
    bool dum_HLT_mu14 = false;
      
    auto chainGroup = m_trigDecisionTool->getChainGroup("HLT_.*");
    for (auto &trig : chainGroup->getListOfTriggers()) {
      auto cg = m_trigDecisionTool->getChainGroup(trig);
      if(cg->isPassed()) {
	std::string thisTrig = trig;
	m_TriggerObject_Chain.push_back(thisTrig);
	m_TriggerObject_Ps.push_back(cg->getPrescale());

	if(trig=="HLT_mu8")  dum_HLT_mu8 = true;
	if(trig=="HLT_mu10") dum_HLT_mu10 = true;
	if(trig=="HLT_mu14") dum_HLT_mu14 = true;

      }
    } 
  
  }

  //  ANA_MSG_INFO ("jet container ");
  const xAOD::JetContainer* jets = 0;
  ANA_CHECK(evtStore()->retrieve( jets, m_reco_jet_collection.c_str())); 
  //  ANA_MSG_INFO ("njets: " << jets->size());
  auto jets_shallowCopy = xAOD::shallowCopyContainer( *jets );
  std::unique_ptr<xAOD::JetContainer> jetsSC (jets_shallowCopy.first);
  std::unique_ptr<xAOD::ShallowAuxContainer> jetsAuxSC (jets_shallowCopy.second);
  //  ANA_MSG_INFO ("jet truth container ");
  const xAOD::JetContainer * truthJets = 0;
  ANA_CHECK(evtStore()->retrieve( truthJets, m_truth_jet_collection.c_str())); 
  auto truthJets_shallowCopy = xAOD::shallowCopyContainer( *truthJets );
  std::unique_ptr<xAOD::JetContainer> truthJetsSC (truthJets_shallowCopy.first);
  std::unique_ptr<xAOD::ShallowAuxContainer> truthJetsAuxSC (truthJets_shallowCopy.second);

  for (unsigned int i = 0; i < jetsSC->size(); i++) {
    xAOD::Jet* jet = (xAOD::Jet*) jetsSC->at(i);
    //    if (m_isPP && !m_jetCleaning->keep( *jet )) continue;

    xAOD::JetFourMom_t jet_4mom = jet->jetP4(); 
    fastjet::PseudoJet MomFourVec = fastjet::PseudoJet ( jet_4mom.px(), jet_4mom.py(), jet_4mom.pz(), jet_4mom.energy() );
    double jetPt = jet_4mom.pt()*1e-3;
    double jetPt2 = jetPt*jetPt;
    double jetEta = jet->eta();
    double jetPhi = jet->phi();
    double jetE = jet->e()*1e-3;
    double jetRap = MomFourVec.rapidity(); 
    double jetMass2 = MomFourVec.m2()*0.001*0.001;
    int truthJetIndex = -1;
    double truthJetPt = -1;
    int truthJet_pid = -1;
    double truthJetPt = -1.;
    double truthJetPt2 = -1.;
    double truthJetEta = -1.;
    double truthJetPhi = -1.;
    double truthJetE = -1.;
    double truthJetRap = -1.; 
    double truthJetMass2 = -1.;

    if ( fabs(jetEta) > m_jet_eta_cut ) 
      continue;    
    
    //if MC, match to truth
    if (m_isMC) {
      truthJetIndex = TruthMatching(jetEta,jetPhi,truthJetsSC);
      //      ANA_MSG_INFO ("truth jet index" << truthJetIndex << " truthjetvecsize " << jetsSC->size()); 
      if (truthJetIndex<0) continue; // truth matching -- cut background
      xAOD::Jet* truthJet = (xAOD::Jet*) truthJetsSC->at(truthJetIndex);
      xAOD::JetFourMom_t truthJet_4mom = truthJet->jetP4(); 
      truthJetPt = truthJet_4mom.pt()*1e-3;
      truthJetPt2 = truthJetPt*truthJetPt;
      truthJetEta = truthJet->eta();
      truthJetPhi = truthJet->phi();
      truthJetE = truthJet->e()*1e-3;
      truthJetRap = truthJet_4mom.rapidity(); 
      truthJetMass2 = truthJet_4mom.m2()*0.001*0.001;
      // LQ: 1-3; gluon: 9 or 21; c: 4; b: 5
      truthJet_pid = truthJet->getAttribute<int>("PartonTruthLabelID");
      if (truthJetPt < m_truth_jet_pT_cut) continue; // cut events with truth below cut

      //    if (!jetcorr->MCJetJERClean(truth_jet_pt_vector.at(truthindex),jet_pt,truth_jet_eta_vector.at(truthindex),cent_bin_fine) ) continue; //cut on JER balance -- cuts events reconstructed 5sigma (ie 5*JER) from truth
    }
  
    m_Jet_pt.push_back(jetPt);    
    m_Jet_pt2.push_back(jetPt2);    
    m_Jet_eta.push_back(jetEta);
    m_Jet_phi.push_back(jetPhi);
    m_Jet_e.push_back(jetE);
    m_Jet_y.push_back(jetRap);
    m_Jet_mass2.push_back(jetMass2);
    m_Jet_mass2OverPt2.push_back(jetMass2/jetPt2);
    m_TruthJet_pid.push_back(truthJet_pid);
    m_TruthJet_pt.push_back(truthJetPt);
    m_TruthJet_pt2.push_back(truthJetPt2);    
    m_TruthJet_eta.push_back(truthJetEta);
    m_TruthJet_phi.push_back(truthJetPhi);
    m_TruthJet_e.push_back(truthJetE);
    m_TruthJet_y.push_back(truthJetRap);
    m_TruthJet_mass2.push_back(truthJetMass2);
    m_TruthJet_mass2OverPt2.push_back(truthJetMass2/truthJetPt2);

    // add prescales if we're using lower trigger HI jets
  }
  
  m_tree->Fill();
  clearVector(); 

  return StatusCode::SUCCESS;
}

void jetMass :: clearVector () {

  m_Jet_pt.clear();
  m_Jet_pt2.clear();
  m_Jet_e.clear();
  m_Jet_eta.clear();
  m_Jet_y.clear();
  m_Jet_phi.clear();
  m_Jet_mass2.clear();
  m_Jet_mass2OverPt2.clear();
  m_TruthJet_pt.clear();
  m_TruthJet_pt2.clear();
  m_TruthJet_e.clear();
  m_TruthJet_eta.clear();
  m_TruthJet_y.clear();
  m_TruthJet_phi.clear();
  m_TruthJet_mass2.clear();
  m_TruthJet_mass2OverPt2.clear();
  
  if (!m_isMC) {
    m_TriggerObject_Chain.clear();
    m_TriggerObject_Ps.clear();
  }

}



double jetMass :: dR(const double eta1,
		     const double phi1,
		     const double eta2,
		     const double phi2)
{
  double deta = fabs(eta1 - eta2);
  double dphi = fabs(phi1 - phi2) < TMath::Pi() ? fabs(phi1 - phi2) : 2*TMath::Pi() - fabs(phi1 - phi2);
  return sqrt(deta*deta + dphi*dphi);
}


StatusCode jetMass :: finalize ()
{
  ANA_MSG_INFO ("in finalize");
  return StatusCode::SUCCESS;
}


int jetMass::TruthMatching(double jetEta, double jetPhi,std::unique_ptr<xAOD::JetContainer> &truthJets) {

  double dR_min=999;
  int truthJetIndex = -1;
  for (unsigned int i = 0; i < truthJets->size(); i++) {

    xAOD::Jet* truthJet = (xAOD::Jet*) truthJets->at(i);
    double truth_eta = truthJet->eta();
    double truth_phi = truthJet->phi();

    double R = dR(truth_eta,truth_phi,jetEta,jetPhi);
    //dR cut
    if(R < m_jet_dR_truth_matching && R < dR_min) {
      dR_min = R;
      truthJetIndex = i;
    }
  }

  return truthJetIndex;	
}

