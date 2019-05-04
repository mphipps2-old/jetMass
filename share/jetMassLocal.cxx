#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>

void jetMassLocal (const std::string& submitDir)
{
  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  int nEvents=-1;
  // create a new sample handler to describe the data files we use
  SH::SampleHandler sh;

  std::string grlFileName="jetMass/data17_5TeV.periodAllYear_DetStatus-v98-pro21-16_Unknown_PHYS_StandardGRL_All_Good_25ns_ignore_GLOBAL_LOWMU.xml";
  std::string reco_jet_collection="AntiKt4HIJets";
  std::string truth_jet_collection="AntiKt4TruthJets";
  double jet_eta_cut = 2.1;
  double jet_pT_cut=40.;
  double jet_dR_truth_matching=0.2;
  double truth_jet_pT_cut=30.;
  int isMC=1;
  int isPP=1;  

  //  const char* inputFilePath = gSystem->ExpandPathName ("/usatlas/u/mphipps2/ML/analysisBase_jetMass/athena/jetMass/data/");
  const char* inputFilePath = gSystem->ExpandPathName ("/usatlas/u/mphipps2/usatlasdata/ML_data/AOD/mc16_5TeV");

  SH::ScanDir().filePattern("*AOD*").scan(sh,inputFilePath);

  // set the name of the tree in our files
  // in the xAOD the TTree containing the EDM containers is "CollectionTree"
  sh.setMetaString ("nc_tree", "CollectionTree");
// further sample handler configuration may go here

// print out the samples we found
   sh.print ();

  // this is the basic description of our job
  EL::Job job;
  job.sampleHandler (sh); // use SampleHandler in this job
  job.options()->setDouble (EL::Job::optMaxEvents, nEvents); // for testing purposes, limit to run over the first 500 events only!
    // define an output and an ntuple associated to that output  
  EL::OutputStream output  ("myOutput");
  job.outputAdd (output);
  //EL::NTupleSvc *ntuple = new EL::NTupleSvc ("myOutput");
  //job.algsAdd (ntuple);

  // add our algorithm to the job
  //job.useXAOD ();  
  EL::AnaAlgorithmConfig config;
  config.setType ("jetMass");
  // set the name of the algorithm (this is the name use with
  // messages)
  config.setName ("jetMass");
  config.setProperty("RecoJetContainer",reco_jet_collection.c_str()).ignore();
  config.setProperty("RecoJetContainer",reco_jet_collection.c_str()).ignore();
  config.setProperty("GRL_FileName",grlFileName.c_str()).ignore();
  config.setProperty("TruthJetContainer",truth_jet_collection.c_str()).ignore();
  config.setProperty("JetEtaCut",jet_eta_cut).ignore();
  config.setProperty("JetPtCut",jet_pT_cut).ignore();
  config.setProperty("Jet_dr_truthMatching",jet_dR_truth_matching).ignore();
  config.setProperty("TruthJetPtCut",truth_jet_pT_cut).ignore();
  config.setProperty("MC_Flag",isMC).ignore();
  config.setProperty("pp_Flag",isPP).ignore();

  // later on we'll add some configuration options for our algorithm that go here
  
  job.algsAdd (config);
  std::cout << " alg added " << std::endl;
  //config.m_outputName = "myOutput"; // give the name of the output to our algorithm

  // make the driver we want to use:
  // this one works by running the algorithm directly:
  EL::DirectDriver driver;
  // we can use other drivers to run things on the Grid, with PROOF, etc.
  // process the job using the driver
  driver.submit (job, submitDir);

}



