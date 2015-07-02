import FWCore.ParameterSet.Config as cms


pfWeightedPhotons = cms.EDProducer('DeltaBetaWeights',
                                   src = cms.InputTag('pfAllPhotonsPFBRECO'),
                                   chargedFromPV = cms.InputTag("pfAllChargedParticlesPFBRECO"),
                                   chargedFromPU = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
)
pfWeightedNeutralHadrons = cms.EDProducer('DeltaBetaWeights',
                                   src = cms.InputTag('pfAllNeutralHadronsPFBRECO'),
                                   chargedFromPV = cms.InputTag("pfAllChargedParticlesPFBRECO"),
                                   chargedFromPU = cms.InputTag("pfPileUpAllChargedParticlesPFBRECO")
)



