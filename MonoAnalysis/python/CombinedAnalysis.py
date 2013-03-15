import FWCore.ParameterSet.Config as cms
 
process = cms.Process("Mpl")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')


## Fitter-smoother: loosen outlier rejection as for first data-taking with LHC "collisions"
process.KFFittingSmootherWithOutliersRejectionAndRK.BreakTrajWith2ConsecutiveMissing = False
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 1000



### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_53_V8::All'

### Track refitter specific stuff
process.load("RecoTracker.TrackProducer.TrackRefitters_cff") #the correct one


process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(10)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'file:/data1/swilbur/Pythia_534_monopole_1000GeV_RECO_RAW.root'
	'file:/data1/swilbur/Pythia_526_monopole_1000GeV_RECO_NoRadiation.root'
    )
)

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process.Monopoler = cms.EDAnalyzer('MonoNtupleDumper'
  ,isData = cms.bool(False)
  ,Output = cms.string("Combined.test.root")
  ,EcalEBRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB") 
  ,ElectronTag = cms.InputTag("gsfElectrons","")
  ,StripSeedLength = cms.uint32(3)
  ,ClusterLength = cms.uint32(5)
  ,SeedThreshold = cms.double(50.)
  ,TrackSource=cms.string("TrackRefitter")
  ,TrackChi2Cut=cms.untracked.double(7.5)
  ,TrackPtCut=cms.untracked.double(3.0)
  ,TrackDeDxCut=cms.untracked.double(0)
  ,TrackDefaultError=cms.untracked.double(0.05)
  ,TrackErrorFudge=cms.untracked.double(0.02)
  ,TrackHitOutput=cms.untracked.bool(True)
)


process.refit_step = cms.Path(process.TrackRefitter)
process.mpl_step = cms.Path(process.Monopoler)

process.options = cms.untracked.PSet(     wantSummary = cms.untracked.bool(True) )


process.p1 = cms.Schedule(process.refit_step,
			process.mpl_step
)
#process.outpath = cms.EndPath(process.TRACKS)
