import FWCore.ParameterSet.Config as cms

process = cms.Process("TagProbe")

process.load('FWCore.MessageService.MessageLogger_cfi')
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource", 
    fileNames = cms.untracked.vstring(
      ""
                                      )
)
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )    

TightMuon_Cut = "isGlobalMuon && isPFMuon" + \
" && globalTrack().normalizedChi2 < 10.0" + \
" && globalTrack().hitPattern().numberOfValidMuonHits > 0" + \
" && globalTrack().hitPattern().numberOfValidPixelHits > 0" + \
" && numberOfMatchedStations>1" + \
" && globalTrack().hitPattern().trackerLayersWithMeasurement > 5" + \
" && dB() < 0.2"

print "TightMuon Definition: " + TightMuon_Cut

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string('MCRUN2_74_V9A::All')

### load default PAT sequence
process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.patMuons.muonSource = 'muons'
process.patMuons.embedTrack = True
process.patMuons.embedCombinedMuon = True
process.patMuons.embedStandAloneMuon = True
process.patMuons.embedPFCandidate = True
process.patMuons.embedCaloMETMuonCorrs = cms.bool(False)
process.patMuons.embedTcMETMuonCorrs = cms.bool(False)

# then switch off some features we don't need
process.patMuons.embedPickyMuon = False
process.patMuons.embedTpfmsMuon = False
process.patMuons.userIsolation = cms.PSet()
process.patMuons.isoDeposits = cms.PSet()
process.patMuons.addGenMatch = False

process.load("Configuration.StandardSequences.Reconstruction_cff") #Needed for patMuonsWithTrigger_cff

#Implement trigger information to PAT muons
process.load("MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff")
from MuonAnalysis.MuonAssociators.patMuonsWithTrigger_cff import useExistingPATMuons
useExistingPATMuons(process, "patMuons")

###Implement PFweighted & PUPPI Isolation on patMuons###

# -- PF-Weighted: 
  ##Generate 'pfWeightedPhotons' & 'pfWeightedNeutralHadrons'.
  ##Sequence: pfDeltaBetaWeightingSequence
process.load('CommonTools.ParticleFlow.deltaBetaWeights_cff')

# -- PUPPI
  ##Generate 'pfPUPPIChargedHadrons', 'pfPUPPINeutralHadrons', 'pfPUPPIPhotons'
  ##Sequence: pfPUPPISequence
from pfPUPPISequence_cff import *
load_pfPUPPI_sequence(process, 'pfPUPPISequence', algo = 'PUPPI',
  src_puppi = 'pfAllHadronsAndPhotonsForPUPPI',
  cone_puppi_central = 0.5
)

# load user-defined muon PF-isolation values
muon_src, cone_size = 'patMuonsWithTrigger', 0.4

from MuonPFIsolationSequence_cff import *

  ## Generate muPFIsoValueCHR04STAND, muPFIsoValueNHR04STAND, muPFIsoValuePhR04STAND and muPFIsoValuePUR04STAND
  ## Sequence: MuonPFIsoSequenceSTAND
load_muonPFiso_sequence(process, 'MuonPFIsoSequenceSTAND', algo = 'R04STAND',
  src = muon_src,
  src_charged_hadron = 'pfAllChargedHadronsPFBRECO',
  src_neutral_hadron = 'pfAllNeutralHadronsPFBRECO',
  src_photon         = 'pfAllPhotonsPFBRECO',
  src_charged_pileup = 'pfPileUpAllChargedParticlesPFBRECO',
  coneR = cone_size
)

  ## Generate muPFIsoValueNHR04PFWGT and muPFIsoValuePhR04PFWGT
  ## Sequence: MuonPFIsoSequencePFWGT
load_muonPFiso_sequence(process, 'MuonPFIsoSequencePFWGT', algo = 'R04PFWGT',
  src = muon_src,
  src_neutral_hadron = 'pfWeightedNeutralHadrons',
  src_photon         = 'pfWeightedPhotons',
  coneR = cone_size
)

  ## Generate muPFIsoValueCHR04PUPPI, muPFIsoValueNHR04PUPPI and muPFIsoValuePhR04PUPPI
  ## Sequence: MuonPFIsoSequencePUPPI
load_muonPFiso_sequence(process, 'MuonPFIsoSequencePUPPI', algo = 'R04PUPPI',
  src = muon_src,
  src_charged_hadron = 'pfPUPPIChargedHadrons',
  src_neutral_hadron = 'pfPUPPINeutralHadrons',
  src_photon         = 'pfPUPPIPhotons',
  coneR = cone_size
)

#Calculate 'Relative' PFWeighted isolation
process.RelIsoPFWeighted = cms.EDProducer("CalculateRelIso",
  muons = cms.InputTag("patMuonsWithTrigger"),
  SumCH = cms.InputTag("muPFIsoValueCHR04STAND"), #Charged Hadrons are not weigthed in PFWeighted isolation
  SumNH = cms.InputTag("muPFIsoValueNHR04PFWGT"),
  SumPh = cms.InputTag("muPFIsoValuePhR04PFWGT")
  )

#Calculate 'Relative' PUPPI isolation
process.RelIsoPUPPI = cms.EDProducer("CalculateRelIso",
  muons = cms.InputTag("patMuonsWithTrigger"),
  SumCH = cms.InputTag("muPFIsoValueCHR04PUPPI"),
  SumNH = cms.InputTag("muPFIsoValueNHR04PUPPI"),
  SumPh = cms.InputTag("muPFIsoValuePhR04PUPPI")
  )

process.MuonPFIsoSequences = cms.Sequence(
  process.MuonPFIsoSequenceSTAND *
  process.MuonPFIsoSequencePFWGT *
  process.MuonPFIsoSequencePUPPI *
  (process.RelIsoPFWeighted + process.RelIsoPUPPI)
)

###Tag&Probe Tree Part###

#Implement dxyPV, dzPV value on patMuons
process.MuonDxyDzPV = cms.EDProducer("MuonDxyDzPV",
  probes = cms.InputTag("patMuonsWithTrigger")
  )

# Tags.
process.tagMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string( "pt>20 && abs(eta)<2.4 && !triggerObjectMatchesByFilter('hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2').empty() &&" + TightMuon_Cut ), 
)

# Probes.
process.probeMuons = cms.EDFilter("PATMuonRefSelector",
    src = cms.InputTag("patMuonsWithTrigger"),
    cut = cms.string( "isTrackerMuon" ), 
)

# Combine Tags and Probes into Z candidates, applying a mass cut
process.tpPairs = cms.EDProducer("CandViewShallowCloneCombiner",
    decay = cms.string("tagMuons@+ probeMuons@-"), # charge coniugate states are implied
    cut   = cms.string("60 < mass < 120"),
)

# Match muons to MC
process.muMcMatch = cms.EDProducer("MCTruthDeltaRMatcherNew",
    pdgId = cms.vint32(13),
    src = cms.InputTag("patMuonsWithTrigger"),
    distMin = cms.double(0.3),
    matched = cms.InputTag("genParticles")
)

# Make the tree
process.muonEffs = cms.EDAnalyzer("TagProbeFitTreeProducer",
    # pairs
    tagProbePairs = cms.InputTag("tpPairs"),
    arbitration   = cms.string("OneProbe"),
    # variables to use
    variables = cms.PSet(
        eta = cms.string("eta"),
        pt  = cms.string("pt"),
        dxy = cms.string("muonBestTrack().dxy()"),
        dz = cms.string("muonBestTrack().dz()"),
        RelPFIso = cms.string("(pfIsolationR04.sumChargedHadronPt+pfIsolationR04.sumNeutralHadronEt+pfIsolationR04.sumPhotonEt) / pt"),
        RelPFIso_PFWeighted = cms.InputTag("RelIsoPFWeighted"),
        RelPFIso_PUPPI = cms.InputTag("RelIsoPUPPI"),
        RelPFIso_dBeta = cms.string("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt"),
        RelTrkIso = cms.string("isolationR03().sumPt()/pt"),
        dzPV = cms.InputTag("MuonDxyDzPV", "dzPV"),
        dxyPV = cms.InputTag("MuonDxyDzPV", "dxyPV")
    ),
    # choice of what defines a 'passing' probe
    flags = cms.PSet(
        ## two defined by simple string cuts
        passingHLT_Mu17Mu8dzIsoVVL = cms.string("!triggerObjectMatchesByFilter('hltDiMuonGlb17Glb8RelTrkIsoFiltered0p4DzFiltered0p2').empty()"),
        passingHLT_Mu17Mu8dzIsoVVL_Filter8 = cms.string("!triggerObjectMatchesByFilter('hltL3pfL1sDoubleMu103p5L1f0L2pf0L3PreFiltered8').empty()"),
        passingHLT_Mu17Mu8dzIsoVVL_Filter17 = cms.string("!triggerObjectMatchesByFilter('hltL3fL1sDoubleMu103p5L1f0L2f10OneMuL3Filtered17').empty()"),
        passingTight = cms.string( TightMuon_Cut ),
        passingGlb = cms.string( "isGlobalMuon" ),
    ),
    # mc-truth info
    isMC = cms.bool(True),
    motherPdgId = cms.vint32(22,23),
    makeMCUnbiasTree = cms.bool(True),
    checkMotherInUnbiasEff = cms.bool(True),
    tagMatches = cms.InputTag("muMcMatch"),
    probeMatches  = cms.InputTag("muMcMatch"),
    allProbes     = cms.InputTag("probeMuons"),
)

##    ____       _   _     
##   |  _ \ __ _| |_| |__  
##   | |_) / _` | __| '_ \ 
##   |  __/ (_| | |_| | | |
##   |_|   \__,_|\__|_| |_|
##                         
process.tagAndProbe = cms.Path(
    process.patCandidates * 
    process.patMuonsWithTriggerSequence *
    process.MuonDxyDzPV *
    (process.tagMuons + process.probeMuons) *   # 'A*B' means 'B needs output of A';
    (process.tpPairs + process.muMcMatch) *
    process.muonEffs
)

#Delete unnecessary sequences
process.tagAndProbe.remove(process.makePatElectrons)
process.tagAndProbe.remove(process.makePatPhotons)
process.tagAndProbe.remove(process.makePatJets)
process.tagAndProbe.remove(process.makePatTaus)
process.tagAndProbe.remove(process.makePatMETs)
process.tagAndProbe.remove(process.patCandidateSummary)

process.tagAndProbe.replace(
  process.pfParticleSelectionPFBRECOSequence,
  process.pfParticleSelectionPFBRECOSequence  *
  process.pfDeltaBetaWeightingSequence *
  process.pfPUPPISequence
)

process.tagAndProbe.replace(
  process.patMuonsWithTriggerSequence,
  process.patMuonsWithTriggerSequence *
  process.MuonPFIsoSequences
)

process.TFileService = cms.Service("TFileService", fileName = cms.string("ROOTFile_TPTree_ZMuMu_RecoID_MC.root"))