#include <EventLoopAlgs/NTupleSvc.h>
#include <EventLoopAlgs/AlgSelect.h>


void jetMassGrid (const std::string& submitDir)
{
  // Set up the job for xAOD access:
  xAOD::Init().ignore();

  // create a new sample handler to describe the data files we use
  SH::SampleHandler sh;
  int nEvents=-1;

  std::string grlFileName="jetMass/data17_5TeV.periodAllYear_DetStatus-v98-pro21-16_Unknown_PHYS_StandardGRL_All_Good_25ns_ignore_GLOBAL_LOWMU.xml";
  std::string reco_jet_collection="AntiKt4HIJets";
  std::string truth_jet_collection="AntiKt4TruthJets";
  double jet_eta_cut = 2.1;
  double jet_pT_cut=40.;
  double jet_dR_truth_matching=0.2;
  double truth_jet_pT_cut=30.;
  int isMC=1;
  int isPP=1;  


  //  string fileName = "mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ0R04.recon.AOD.e4108_s3238_r11199/";
//  string fileName = "mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e4108_s3238_r11199/";
//  string fileName = "mc16_5TeV.420011.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1R04.recon.AOD.e6608_s3238_r11199/";
  string fileName = "mc16_5TeV:mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e4108_s3238_r11199/";
//  string fileName = "mc16_5TeV.420012.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04.recon.AOD.e6608_s3238_r11199/";
//  string fileName = "mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e4108_s3238_r11199/";
//  string fileName = "mc16_5TeV.420013.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04.recon.AOD.e6608_s3238_r11199/";
// // HERWIG
//  string fileName = "mc16_5TeV:mc16_5TeV.364444.Herwig7EvtGen_H7UE_NNPDF30nlo_jetjet_JZ1.recon.AOD.e7288_s3238_r11199/";
//  string fileName = "mc16_5TeV:mc16_5TeV.364445.Herwig7EvtGen_H7UE_NNPDF30nlo_jetjet_JZ2.recon.AOD.e7288_s3238_r11199/";
//  string fileName = "mc16_5TeV:mc16_5TeV.364446.Herwig7EvtGen_H7UE_NNPDF30nlo_jetjet_JZ3.recon.AOD.e7288_s3238_r11199/";
//  // MC pp muon filter //
//  string fileName = "mc16_5TeV:mc16_5TeV.420267.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ3R04Opt_mufilter.merge.AOD.e7190_s3238_r10441_r10210/";
//  string fileName = "mc16_5TeV:mc16_5TeV.420265.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1BR04Opt_mufilter.merge.AOD.e7190_s3238_r10441_r10210"); 
//  string fileName = "mc16_5TeV:mc16_5TeV.420264.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1AR04_mufilter.merge.AOD.e7189_s3238_r10441_r10210/";  
//  string fileName = "mc16_5TeV:mc16_5TeV.420266.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ2R04Opt_mufilter.merge.AOD.e7190_s3238_r10441_r10210/";


//  string fileName = "data17_5TeV:data17_5TeV.00*.physics_Main.merge.AOD.f911_m1917/";
//  const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/s/stapiaar/Jetdata/ttbar_signal/group.phys-hi.mc16_valid.410503.PowhegPythia8EvtGen_A14_ttbar_hdamp258p75_dil_r20180318T2154_HI_mc15c_02_EXT1/");
 //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/s/stapiaar/JetB_tag/run/data/");
  //const char* inputFilePath = gSystem->ExpandPathName ("/afs/cern.ch/work/s/stapiaar/Jetdata/Hion2015/data15_hi.00287931.physics_HardProbes.merge.AOD.r7874_p2580");
  //SH::ScanDir().filePattern("*AOD*").scan(sh,inputFilePath);
  //SH::scanRucio (sh, "data18_hi:data18_hi.*.physics_HardProbes.merge.AOD.f1021_m2037");
  //SH::scanRucio (sh, "data18_hi:data18_hi.*.physics_MinBias.merge.AOD.f1021_m2037");
 // SH::scanRucio (sh, "data18_hi:data18_hi.00367233.physics_HardProbes.merge.AOD.f1030_m2048");
  //SH::scanRucio (sh, "data18_hi:data18_hi.00366805.physics_CC.merge.AOD.f1030_m2048");
  //SH::scanRucio (sh, "data18_hi:data18_hi.00366805.physics_PC.merge.AOD.f1030_m2048");

  SH::scanRucio (sh, fileName.c_str());  
  // set the name of the tree in our files
  // in the xAOD the TTree containing the EDM containers is "CollectionTree"
  sh.setMetaString ("nc_tree", "CollectionTree");
  sh.setMetaString( "nc_grid_filter", "*AOD*");
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

  // set the name of the algorithm (this is the name use with messages)
  config.setName ("jetMass");
  config.setProperty("GRL_FileName",grlFileName.c_str()).ignore();
  config.setProperty("RecoJetContainer",reco_jet_collection.c_str()).ignore();
  config.setProperty("RecoJetContainer",reco_jet_collection.c_str()).ignore();
  config.setProperty("TruthJetContainer",truth_jet_collection.c_str()).ignore();
  config.setProperty("JetEtaCut",jet_eta_cut).ignore();
  config.setProperty("JetPtCut",jet_pT_cut).ignore();
  config.setProperty("Jet_dr_truthMatching",jet_dR_truth_matching).ignore();
  config.setProperty("TruthJetPtCut",truth_jet_pT_cut).ignore();
  config.setProperty("MC_Flag",isMC).ignore();
  config.setProperty("pp_Flag",isPP).ignore();
  // later on we'll add some configuration options for our algorithm that go here
  
  job.algsAdd (config);

  //config.m_outputName = "myOutput"; // give the name of the output to our algorithm

// make the driver we want to use:
// this one works by running the algorithm directly:
  //EL::DirectDriver driver;
  EL::PrunDriver driver;
 driver.options()->setString("nc_outputSampleName", "user.mphipps.bTagging.V1.%in:name[2]%.%in:name[6]%");
  // we can use other drivers to run things on the Grid, with PROOF, etc.
    // process the job using the driver
 //driver.submit (job, submitDir);
 driver.submitOnly (job, submitDir);

}



