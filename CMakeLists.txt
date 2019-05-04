# The name of the package:
atlas_subdir (jetMass)

#atlas_depends_on_subdirs( PUBLIC
#                          InnerDetector/InDetRecTools/InDetTrackSelectionTool
#                          Tools/PathResolver
#                          )


# Add the shared library:
atlas_add_library (jetMassLib
  jetMass/*.h Root/*.cxx
  PUBLIC_HEADERS jetMass
  LINK_LIBRARIES JetUncertaintiesLib PathResolver ParticleJetToolsLib JetCalibToolsLib InDetTrackSelectionToolLib MuonAnalysisInterfacesLib TrigDecisionToolLib TriggerMatchingToolLib TrigConfInterfaces AnaAlgorithmLib EventLoopAlgs EventLoop xAODHIEvent xAODEgamma xAODMuon xAODJet xAODEventInfo AsgAnalysisInterfaces)

# Add the dictionary:
atlas_add_dictionary (jetMassDict
  jetMass/jetMassDict.h
  jetMass/selection.xml
  LINK_LIBRARIES jetMassLib)

if (NOT XAOD_STANDALONE)
  # Add a component library for AthAnalysis only:
  atlas_add_component (jetMass
    src/components/*.cxx
    LINK_LIBRARIES jetMassLib)
endif ()

# Install files from the package:
atlas_install_joboptions( share/*.py )
atlas_install_data( data/* )
