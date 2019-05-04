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
  m_muonSelection ("CP::MuonSelectionTool", this),
  m_muonCalibrationAndSmearingTool ("CP::MuonCalibrationAndSmearingTool/MuonCorrectionTool",this),
    m_trkSelection ("InDet::InDetTrackSelectionTool/MyTrackTool", this),
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
  declareProperty("TrackEtaCut",m_track_eta_cut=2.5,"Max track eta cut");
  declareProperty("TrackPtCut",m_track_pt_cut="HITight","Track selector track pt cut");
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
  //  m_tree->Branch ("LumiBlock", &m_LumiBlock);
  //  m_tree->Branch ("Fcal", &m_FCal_Et);
  //  m_tree->Branch ("FcalWeight", &m_FCalWeight);
  m_tree->Branch ("JzNorm", &m_jzNorm);
  m_tree->Branch ("PrimaryVertexZ", &m_primaryVertexZ);

  // Antikt4HI jets
  m_tree->Branch ("Jet_pt", &m_Jet_pt);
  //  m_tree->Branch ("Jet_e",  &m_Jet_e);
  //  m_tree->Branch ("Jet_eta", &m_Jet_eta);
  m_tree->Branch ("Jet_y", &m_Jet_y);
  m_tree->Branch ("Jet_phi", &m_Jet_phi);

  // Antikt4HI truth jets
  m_tree->Branch ("TruthJet_pid", &m_TruthJet_pid);
  m_tree->Branch ("TruthJet_pt", &m_TruthJet_pt);
  //  m_tree->Branch ("TruthJet_e",  &m_TruthJet_e);
  //  m_tree->Branch ("TruthJet_eta", &m_TruthJet_eta);
  //  m_tree->Branch ("TruthJet_y", &m_TruthJet_y);
  //  m_tree->Branch ("TruthJet_phi", &m_TruthJet_phi);

  // Secondary vertex associated with jet
  m_tree->Branch ("JetSV1_pu", &m_Jet_sv1_pu);
  m_tree->Branch ("JetSV1_pb", &m_Jet_sv1_pb);
  m_tree->Branch ("JetSV1_pc", &m_Jet_sv1_pc);
  m_tree->Branch ("JetSV1_efrc", &m_Jet_sv1_efrc);
  m_tree->Branch ("JetSV1_mass", &m_Jet_sv1_mass);
  m_tree->Branch ("JetSV1_n2t", &m_Jet_sv1_n2t);
  m_tree->Branch ("JetSV1_ntrkv", &m_Jet_sv1_ntrkv);
  m_tree->Branch ("JetSV1_Lxy", &m_Jet_sv1_LXY);
  m_tree->Branch ("JetSV1_L3d", &m_Jet_sv1_L3d);
  m_tree->Branch ("JetSV1_sig3d", &m_Jet_sv1_sig3d);
  m_tree->Branch ("JetSV1_distmatlay", &m_Jet_sv1_distmatlay);
  m_tree->Branch ("JetSV1_dR", &m_Jet_sv1_dR);

  // Inner Detector Tracks associated with jet
  m_tree->Branch ("JetTrack_d0", &m_JetTrack_d0);
  m_tree->Branch ("JetTrack_d0_cut", &m_JetTrack_d0_cut);
  m_tree->Branch ("JetTrack_dR", &m_JetTrack_dR);
  m_tree->Branch ("JetTrack_charge", &m_JetTrack_charge);
  m_tree->Branch ("JetTrack_pt", &m_JetTrack_pt);
  m_tree->Branch ("JetTrack_pz", &m_JetTrack_pz);
  m_tree->Branch ("JetTrack_e", &m_JetTrack_e); 
  m_tree->Branch ("JetTrack_eta", &m_JetTrack_eta);
  m_tree->Branch ("JetTrack_phi", &m_JetTrack_phi);
  m_tree->Branch ("JetTrack_pid", &m_JetTrack_pid);
  m_tree->Branch ("JetTrack_truth_dR", &m_JetTrack_truth_dR);
  m_tree->Branch ("JetTrack_truthMatched", &m_JetTrack_truthMatched);

  /*
  m_tree->Branch ("JetTrack_status", &m_JetTrack_status);
  m_tree->Branch ("JetTrack_type", &m_JetTrack_type);
  m_tree->Branch ("JetTrackTruth_charge", &m_JetTrackTruth_charge);
  m_tree->Branch ("JetTrackTruth_pt", &m_JetTrackTruth_pt);
  m_tree->Branch ("JetTrackTruth_pz", &m_JetTrackTruth_pz);
  m_tree->Branch ("JetTrackTruth_e", &m_JetTrackTruth_e); 
  m_tree->Branch ("JetTrackTruth_eta", &m_JetTrackTruth_eta);
  m_tree->Branch ("JetTrackTruth_phi", &m_JetTrackTruth_phi);
  */

  // MUON
  m_tree->Branch ("Muon_pt",     &m_Muon_pt);
  m_tree->Branch ("Muon_eta",    &m_Muon_eta);
  m_tree->Branch ("Muon_phi",    &m_Muon_phi);
  m_tree->Branch ("Muon_charge", &m_Muon_charge);
  m_tree->Branch ("Muon_parent",  &m_Muon_parent);
  m_tree->Branch ("Muon_e",      &m_Muon_e);
  m_tree->Branch ("Muon_quality",&m_Muon_quality);
 
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
  // muon selector tool
  ANA_CHECK (m_muonSelection.setProperty("TrtCutOff",true));
  ANA_CHECK (m_muonSelection.setProperty("MaxEta", m_track_eta_cut)); 
  ANA_CHECK (m_muonSelection.setProperty("MuQuality", 1));
  ANA_CHECK (m_muonSelection.initialize());
  // muon calibration and smearing
  ANA_CHECK (m_muonCalibrationAndSmearingTool.initialize());
  // tracking selection tool
  ANA_CHECK (m_trkSelection.setProperty("CutLevel",m_track_pt_cut.c_str()));
  // ANA_CHECK (m_trkSelection.setProperty("maxZ0SinTheta",1.0));
  //  ANA_CHECK (m_trkSelection.setProperty("minPt",m_reco_track_pT_cut));
  // ANA_CHECK (m_trkSelection.setProperty("maxNSiSharedModules",100));
  ANA_CHECK (m_trkSelection.initialize());

  m_d0_cut = new TF1("f1", "[0]*exp([1]*x)+[2]*exp([3]*x)", 0.4, 500);
  m_d0_cut->SetParameters(0.472367, -0.149934, 0.193095, 0.000337765);

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
  

    if( !(dum_HLT_mu8 || dum_HLT_mu10 || dum_HLT_mu14)  ) return StatusCode::SUCCESS;
  }
  //  ANA_MSG_INFO ("Muon container");
  const xAOD::MuonContainer* muons = 0;
  ANA_CHECK(evtStore()->retrieve( muons, "Muons" ));
 
  auto muons_shallowCopy = xAOD::shallowCopyContainer( *muons );
  std::unique_ptr<xAOD::MuonContainer> muonsSC (muons_shallowCopy.first);
  std::unique_ptr<xAOD::ShallowAuxContainer> muonsAuxSC (muons_shallowCopy.second);

  for (unsigned int i = 0; i < muonsSC->size(); i++) {
    xAOD::Muon* muonSC0 = (xAOD::Muon*) muonsSC->at(i);
    if(m_muonCalibrationAndSmearingTool->applyCorrection(*muonSC0)!= CP::CorrectionCode::Ok){
      ANA_MSG_INFO ("execute(): Problem with Muon Calibration And Smearing Tool (Error or OutOfValidityRange) ");
    }
    if ( !m_muonSelection->accept(*muonSC0) ) continue;

    // start probe loop
    for (unsigned int j = i + 1; j < muonsSC->size(); j++) {
      xAOD::Muon* muonSC1 = (xAOD::Muon*) muonsSC->at(j);
      if(m_muonCalibrationAndSmearingTool->applyCorrection(*muonSC1)!= CP::CorrectionCode::Ok){
	ANA_MSG_INFO ("execute(): Problem with Muon Calibration And Smearing Tool (Error or OutOfValidityRange) ");
      }	
      if ( !m_muonSelection->accept(*muonSC1) ) continue;
      if (muonSC1->charge() == muonSC0->charge()) continue;

      TLorentzVector daughter0;
      TLorentzVector daughter1;
      TLorentzVector mom;
      daughter0.SetPtEtaPhiM(muonSC0->pt()/1.e3, muonSC0->eta(), muonSC0->phi(), 0.106);
      daughter1.SetPtEtaPhiM(muonSC1->pt()/1.e3, muonSC1->eta(), muonSC1->phi(), 0.106);
      mom = daughter0 + daughter1;
      if (mom.M() < 2.7 or mom.M() > 3.5) continue;
      //Probe info
      m_Muon_pt       .push_back( muonSC1->pt()*1e-3);
      m_Muon_eta      .push_back( muonSC1->eta() );
      m_Muon_phi      .push_back( muonSC1->phi() );
      m_Muon_charge   .push_back( muonSC1->charge() );
      m_Muon_e        .push_back( muonSC1->e()*1e-3 );
      m_Muon_quality  .push_back( m_muonSelection->getQuality(*muonSC1) );
    } // end loop probe muon
	
  } // end loop Tag


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

  bool isCJet = 0;
  bool isBJet = 0; 
  //  ANA_MSG_INFO ("loop through jets: " << jetsSC->size());
  //  for (xAOD::Jet *jet : *jets) {
  for (unsigned int i = 0; i < jetsSC->size(); i++) {
    xAOD::Jet* jet = (xAOD::Jet*) jetsSC->at(i);
    //    ANA_MSG_INFO ("new jet ");
    //    if (m_isPP && !m_jetCleaning->keep( *jet )) continue;

    xAOD::JetFourMom_t jet_4mom = jet->jetP4(); 
    fastjet::PseudoJet MomFourVec = fastjet::PseudoJet ( jet_4mom.px(), jet_4mom.py(), jet_4mom.pz(), jet_4mom.energy() );
    double jetPt = jet_4mom.pt()*1e-3;
    double jetEta = jet->eta();
    double jetPhi = jet->phi();
    double jetE = jet->e()*1e-3;
    double jetRap = MomFourVec.rapidity(); 
    int truthJetIndex = -1;
    double truthJetPt = -1;
    int truthJet_pid = -1;
    if ( fabs(jetEta) > m_jet_eta_cut ) 
      continue;

    
    float sv1_Lxy = -99.; float sv1_L3d = -99.; float sv1_mass = -99.; float sv1_efrc = -99.; 
    int sv1_n2t = -99; int sv1_ntrkv = -99; 
    float sv1_sig3d = -99.; float sv1_distmatlay = -99.; float sv1_dR = -99.;
    double sv1_pu = -99.; double sv1_pb = -99.; double sv1_pc = -99.; 

    const xAOD::BTagging *bjet(nullptr);
    bjet = jet->btagging();
    bool sv1_ok(false);
    std::vector< ElementLink< xAOD::VertexContainer > > myVertices_SV1;
    bjet->variable<std::vector<ElementLink<xAOD::VertexContainer> > >("SV1", "vertices", myVertices_SV1);
    if ( myVertices_SV1.size() > 0 && myVertices_SV1[0].isValid() ) {
      // if we found a vertex, then sv1 is okay to use
      sv1_ok = true;
    }
    if (sv1_ok) {
      bjet->variable<double>("SV1", "pu", sv1_pu);
      bjet->variable<double>("SV1", "pb", sv1_pb);
      bjet->variable<double>("SV1", "pc", sv1_pc);
            
      bjet->variable<float>("SV1", "masssvx",  sv1_mass);
      bjet->variable<float>("SV1", "efracsvx", sv1_efrc);
          bjet->variable<int>("SV1",   "N2Tpair",  sv1_n2t);
           bjet->variable<int>("SV1",   "NGTinSvx", sv1_ntrkv);
      bjet->variable<float>("SV1", "normdist", sv1_sig3d);
        
      bjet->variable<float>("SV1", "dstToMatLay" , sv1_distmatlay);
      bjet->variable<float>("SV1", "deltaR", sv1_dR);
      bjet->variable<float>("SV1", "Lxy",    sv1_Lxy);
      bjet->variable<float>("SV1", "L3d",    sv1_L3d);
    }
    
    //if MC, match to truth
    if (m_isMC) {
      truthJetIndex = TruthMatching(jetEta,jetPhi,truthJetsSC);
      //      ANA_MSG_INFO ("truth jet index" << truthJetIndex << " truthjetvecsize " << jetsSC->size()); 
      if (truthJetIndex<0) continue; // truth matching -- cut background
      xAOD::Jet* truthJet = (xAOD::Jet*) truthJetsSC->at(truthJetIndex);
      xAOD::JetFourMom_t truthJet_4mom = truthJet->jetP4(); 
      truthJetPt = truthJet_4mom.pt()*1e-3;
      truthJet_pid = truthJet->getAttribute<int>("PartonTruthLabelID");
      // LQ: 1-3; gluon: 9 or 21; c: 4; b: 5
      if (truthJet_pid == 1 || truthJet_pid == 2 || truthJet_pid == 3 || truthJet_pid == 9 || truthJet_pid == 21) {
	// keep 5% of total inclusive events
      }
      else if (truthJet_pid == 4) {
	// keep 50% of c jets (10% of total spectrum)
	isCJet = 1;
      }
      else if (truthJet_pid == 5) {
	// keep all b jets (around 5% of total spectrum)
	isBJet = 1; 
      }

      if (truthJetPt < m_truth_jet_pT_cut) continue; // cut events with truth below cut

      //    if (!jetcorr->MCJetJERClean(truth_jet_pt_vector.at(truthindex),jet_pt,truth_jet_eta_vector.at(truthindex),cent_bin_fine) ) continue; //cut on JER balance -- cuts events reconstructed 5sigma (ie 5*JER) from truth
    }
    //    ANA_MSG_INFO ("truth matched jet"); 
  
    ////////////////////// ASSOCIATED TRACKS PER JET //////////////////////
    const xAOD::TrackParticleContainer* recoTracks = 0;
    std::vector<double> jetTrack_dR;
    std::vector<double> jetTrack_d0;
    std::vector<double> jetTrack_d0_cut;
    std::vector<double> jetTrack_pt;
    std::vector<double> jetTrack_pz;
    std::vector<double> jetTrack_e;
    std::vector<double> jetTrack_eta;
    std::vector<double> jetTrack_phi;
    std::vector<int>    jetTrack_charge;
    std::vector<double> jetTrack_truth_dR;
    std::vector<bool>   jetTrack_truthMatched;
    std::vector<int> jetTrack_pid;
    ANA_CHECK(evtStore()->retrieve( recoTracks, "InDetTrackParticles" ));
    //  ANA_MSG_INFO ("looping through tracks "); 
    for (const auto& trk : *recoTracks) {
      //get the tracks....
      double pt = trk->pt()/1000.;
      double pz = trk->p4().Pz()/1000.;
      double eta = trk->eta();
      double e = trk->e();
      double phi = trk->phi();       
      int charge = trk->charge();
      if ( fabs(eta) > m_track_eta_cut ) continue;
      //      ANA_MSG_INFO ("survived track eta cut");
      double dR_track = dR(eta,phi,jetEta,jetPhi);
      double d0 = trk->d0();
      double d0_cut = m_d0_cut->Eval(pt);
      //       if(fabs(d0) > d0_cut) continue; //pT dependant d0 cut
      if(!m_trkSelection->accept(*trk)) continue; //track selector tool

      //Additional track selection: this would happen if the a track is misreconstructed high or a jet misreconstructed low -- note this cut would only cut a single track which would be wrong for my recNN case -- better solution would be to cut the entire event -- removing this for now but study this effect later offline!
      //       if (!trkcorr->PassTracktoJetBalance(pt, jet_pt, eta, jet_eta,cent_bin)) continue;

       // truth particle matching 
       bool isTruthMatched=false;
       fastjet::PseudoJet matchPar = fastjet::PseudoJet (0,0,0,0);
       int pid = 0.;
       //       int status = 0.; 
       //       int trktype = -1;
       //       int charge_truth = 0 ;
       //       double pt_truth = 0.;
       //       double pz_truth = 0.;
       //       double m_truth = 139.570; // pion mass in MeV
       //       double e_truth = 0.;
       double eta_truth = 0.;
       double phi_truth = 0.;
       double dR_truth = 0.;
       if ( m_isMC) { 
	 //Truth matching
	 ElementLink< xAOD::TruthParticleContainer > truthLink = trk->auxdata<ElementLink< xAOD::TruthParticleContainer > >("truthParticleLink");
	 float mcprob;
	 if(truthLink.isValid()) {
	   /*
	     trktype = getTypeReco((*truthLink)->barcode(),(*truthLink)->pdgId(),(*truthLink)->status(),(*truthLink)->charge(),mcprob,_mcProbCut);
	     // 0 - fake, 1 - primary, 2 - secondary, 3 - primary out-of-phase-space, 4 - secondary neutral, 5 - truth primary strange baryons

	     // status 1 is stable particle -- only save status 1
	     status = (*truthLink)->status();
	     pt_truth = (*truthLink)->pt()/1000.;
	     pz_truth = (*truthLink)->p4().Pz()/1000.;
	     e_truth = (*truthLink)->e();
	     charge_truth = (*truthLink)->charge();
	   */
	     pid = (*truthLink)->pdgId();
	     eta_truth = (*truthLink)->eta();
	     phi_truth = (*truthLink)->phi();
	     dR_truth = dR(eta,phi,eta_truth,phi_truth);
	     mcprob = trk->auxdata<float>("truthMatchProbability");
	     // pythia has barcode 0-10k, hijing 10k-200k, secondaries 200k+ 
	     if ( ( mcprob > m_truth_track_prob_cut ) && ( (*truthLink)->barcode() > 0 ) &&  ( (*truthLink)->barcode() < 200000 ) ) {  
	       isTruthMatched = true;
	       matchPar = fastjet::PseudoJet( (*truthLink)->p4() );
	     }
	 }
       }
       jetTrack_dR.push_back(dR_track);
       jetTrack_d0.push_back(d0);
       jetTrack_d0_cut.push_back(d0_cut);
       jetTrack_pt.push_back(pt);
       jetTrack_pz.push_back(pz);
       jetTrack_e.push_back(e) ;
       jetTrack_eta.push_back(eta);
       jetTrack_phi.push_back(phi);       
       jetTrack_charge.push_back(charge);
       jetTrack_truth_dR.push_back(dR_truth);
       jetTrack_truthMatched.push_back(isTruthMatched);
       jetTrack_pid.push_back(pid);
       /*
       jetTrackTruth_pt.push_back(pt_truth);
       jetTrackTruth_pt.push_back(pz_truth);		 
       jetTrackTruth_e.push_back(e_truth) ;
       jetTrackTruth_eta.push_back(eta_truth);
       jetTrackTruth_phi.push_back(phi_truth);
       jetTrack_type.push_back(trktype);
       jetTrack_status.push_back(status);
       */
    } // end reco track loop
    m_JetTrack_dR.push_back(jetTrack_dR);
    m_JetTrack_d0.push_back(jetTrack_d0);
    m_JetTrack_d0_cut.push_back(jetTrack_d0_cut);
    m_JetTrack_charge.push_back(jetTrack_charge);
    m_JetTrack_pt.push_back(jetTrack_pt);    
    m_JetTrack_pz.push_back(jetTrack_pz);
    m_JetTrack_e.push_back(jetTrack_e);
    m_JetTrack_eta.push_back(jetTrack_eta);
    m_JetTrack_phi.push_back(jetTrack_phi);
    m_JetTrack_pid.push_back(jetTrack_pid);
    m_JetTrack_truth_dR.push_back(jetTrack_truth_dR);
    m_JetTrack_truthMatched.push_back(jetTrack_truthMatched);
    /*
    m_JetTrack_status.push_back(jetTrack_status);
    m_JetTrack_type.push_back(jetTrack_type);    
    m_JetTrackTruth_charge.push_back(jetTrackTruth_charge);
    m_JetTrackTruth_pt.push_back(jetTrackTruth_pt);
    m_JetTrackTruth_e.push_back(jetTrackTruth_e);
    m_JetTrackTruth_eta.push_back(jetTrackTruth_eta);
    m_JetTrackTruth_phi.push_back(jetTrackTruth_phi);
    */
    m_Jet_pt.push_back(jetPt);    
    m_Jet_eta.push_back(jetEta);
    m_Jet_phi.push_back(jetPhi);
    m_Jet_e.push_back(jetE);
    m_Jet_y.push_back(jetRap);
    m_TruthJet_pid.push_back(truthJet_pid);
    m_TruthJet_pt.push_back(truthJetPt);

    m_Jet_sv1_pu.push_back(sv1_pu);
    m_Jet_sv1_pb.push_back(sv1_pb);
    m_Jet_sv1_pc.push_back(sv1_pc);
    m_Jet_sv1_mass.push_back(sv1_mass);
    m_Jet_sv1_efrc.push_back(sv1_efrc);
    m_Jet_sv1_n2t.push_back(sv1_n2t);
    m_Jet_sv1_ntrkv.push_back(sv1_ntrkv);
    m_Jet_sv1_sig3d.push_back(sv1_sig3d);
    m_Jet_sv1_distmatlay.push_back(sv1_distmatlay);
    m_Jet_sv1_dR.push_back(sv1_dR);
    m_Jet_sv1_LXY.push_back(sv1_Lxy);
    m_Jet_sv1_L3d.push_back(sv1_L3d);
    // add prescales if we're using lower trigger HI jets
  }
  //  srand(time(NULL));
  //  TRandom3 *rand = new TRandom3(rand());
  TRandom3 *rand = new TRandom3(0);
  double randRatio;
  bool keepInclusiveEvent = 0;
  bool keepCEvent = 0;
  // keep 5% of inclusive events
  if (!isBJet && !isCJet) {
    randRatio = rand->Uniform(0,1);  // returns uniform random # between [0,1]
    if (randRatio > 0.4 && randRatio < 0.45) {
      keepInclusiveEvent = 1;
    }
  }
  // keep 50% of c jet events
  if (isCJet) {
    randRatio = rand->Uniform(0,1);
    if (randRatio < 0.5) {
      keepCEvent = 1;
    }
  }
  if (isBJet || keepCEvent || keepInclusiveEvent) {
    m_tree->Fill();
  }
  clearVector(); 

  return StatusCode::SUCCESS;
}

void jetMass :: clearVector () {
  m_JetTrack_dR.clear();
  m_JetTrack_d0.clear();
  m_JetTrack_d0_cut.clear();
  m_JetTrack_charge.clear();
  m_JetTrack_pt.clear();
  m_JetTrack_e.clear();
  m_JetTrack_eta.clear();
  m_JetTrack_phi.clear();
  m_JetTrack_pid.clear();
  m_JetTrack_truth_dR.clear();
  m_JetTrack_truthMatched.clear();
  m_Jet_pt.clear();
  m_Jet_e.clear();
  m_Jet_eta.clear();
  m_Jet_y.clear();
  m_Jet_phi.clear();
  m_Muon_pt.clear();
  m_Muon_eta.clear();
  m_Muon_phi.clear();
  m_Muon_charge.clear();
  m_Muon_parent.clear();
  m_Muon_quality.clear();
  m_Muon_e.clear();

  m_Jet_sv1_pu.clear();
  m_Jet_sv1_pb.clear();
  m_Jet_sv1_pc.clear();
  m_Jet_sv1_mass.clear();
  m_Jet_sv1_efrc.clear();
  m_Jet_sv1_n2t.clear();
  m_Jet_sv1_ntrkv.clear();
  m_Jet_sv1_sig3d.clear();
  m_Jet_sv1_distmatlay.clear();
  m_Jet_sv1_dR.clear();
  m_Jet_sv1_LXY.clear();
  m_Jet_sv1_L3d.clear();
  
  if (!m_isMC) {
    m_TriggerObject_Chain.clear();
    m_TriggerObject_Ps.clear();
  }

}


double jetMass :: mindRmuon(const xAOD::TrackParticle* trk)
{
  double dr_min = 999;
  const xAOD::MuonContainer* muons2 = 0;
  if(evtStore()->retrieve( muons2, "Muons" ).isFailure()){
    std::cout << "Could not retrieve MuonContainer with key " << "Muons" << std::endl;
  }

  auto muons_shallowCopy2 = xAOD::shallowCopyContainer( *muons2 );
  std::unique_ptr<xAOD::MuonContainer> muonsSC2 (muons_shallowCopy2.first);
  std::unique_ptr<xAOD::ShallowAuxContainer> muonsAuxSC2 (muons_shallowCopy2.second);

  for (unsigned int i = 0; i < muonsSC2->size(); i++) {

    xAOD::Muon* muon2 = (xAOD::Muon*) muonsSC2->at(i);
    if(m_muonCalibrationAndSmearingTool->applyCorrection(*muon2)!= CP::CorrectionCode::Ok){
      ANA_MSG_INFO ("execute(): Problem with Muon Calibration And Smearing Tool (Error or OutOfValidityRange) ");
    }
    if ( !m_muonSelection->accept(*muon2) ) continue;
    if( muon2->charge() != trk->charge()) continue;

    double dr = dR(trk->eta(), trk->phi(), muon2->eta(), muon2->phi());

    if(dr<dr_min) dr_min = dr;
  }

  return dr_min;
}

double jetMass :: mindRTrk(const xAOD::TrackParticle* metrk)
{  
  double dr_min = 999;
  const xAOD::TrackParticleContainer* tracks = 0;
  if(evtStore()->retrieve( tracks, "InDetTrackParticles" ).isFailure()){
    std::cout << "Could not retrieve MuonContainer with key " << "InDetTrackParticles" << std::endl;
  }

  for (auto track: *tracks) {
    if ( !m_muonSelection->passedIDCuts(*track) ) continue;
    if( track->charge() != metrk->charge()) continue;

    double dr = dR(metrk->eta(), metrk->phi(), track->eta(), track->phi());
    if(dr<dr_min) dr_min = dr;
  }
  
  return dr_min;
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



int jetMass :: parent_classify(const xAOD::TruthParticle *theParticle) {
  const xAOD::TruthParticle *parent = 0; // the parent object
  Int_t particle_id = 999;
  Int_t parent_id = 999;

  if (theParticle == NULL) return parent_id;

  particle_id = theParticle->pdgId();
  parent = theParticle->parent(0);
  if (parent) parent_id = parent->pdgId();
  else return parent_id;

  while (fabs(parent_id) == fabs(particle_id) && fabs(parent_id) < 400 && fabs(parent_id) != 0) {
    parent = parent->parent(0);
    if (parent) parent_id = parent->pdgId();
    else break;
  }
  return parent_id;
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

