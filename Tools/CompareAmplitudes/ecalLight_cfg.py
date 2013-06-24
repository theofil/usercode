import FWCore.ParameterSet.Config as cms

process = cms.Process("chi2analysis")

# StandardSequences configurations
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# Global Tag: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Global_Run_data
process.GlobalTag.globaltag = 'GR_P_V17::All' # validity 311X - 41X   Prompt reconstruction Global Tags 

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = "WARNING"

#process.MessageLogger = cms.Service("MessageLogger",
#       destinations   = cms.untracked.vstring(
#                                             'detailedInfo'
#                                               ,'critical'
#                                               ,'cerr'
#                    ),
#       critical       = cms.untracked.PSet(
#                        threshold = cms.untracked.string('ERROR')
#        ),
#       detailedInfo   = cms.untracked.PSet(
#                        threshold = cms.untracked.string('INFO')
#       ),
#       cerr           = cms.untracked.PSet(
#                       threshold = cms.untracked.string('WARNING')
#        )
#)


#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    useJobReport = cms.untracked.bool(True)
#)

process.TFileService = cms.Service("TFileService",
        fileName = cms.string("beam7TeV.root")
        # closeFileFast = cms.untracked.bool(True)
)


#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(2000))
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
 ###   'file:/data/theofil08/Data/Beam10/bit40or41skim_startingFromSkim2.root'
##"/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/716/547F5FB2-7E42-DF11-B366-003048D46296.root"
#"/store/data/Commissioning10/MinimumBias/RAW-RECO/v8/000/132/601/C4A9EF9C-F240-DF11-A5F3-002481E150FC.root"
"/store/data/Run2011A/DoubleElectron/RAW/v1/000/161/217/50F7D1B6-E954-E011-A438-000423D98950.root"
#'file:~/skimmed/r132601eve2864282_4_1.root'
  )
)
process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
process.source.inputCommands = cms.untracked.vstring('drop *','keep *FEDRawData*_*_*_*')


process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

process.load('RecoLocalCalo/EcalRecProducers/ecalRatioUncalibRecHit_cfi')
import RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi
process.ecalRecHitRatio = RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi.ecalRecHit.clone()
process.ecalRecHitRatio.EBuncalibRecHitCollection = 'ecalRatioUncalibRecHit:EcalUncalibRecHitsEB'
process.ecalRecHitRatio.EEuncalibRecHitCollection = 'ecalRatioUncalibRecHit:EcalUncalibRecHitsEE'
process.ecalRecHitRatio.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]


process.load('RecoLocalCalo/EcalRecProducers/ecalGlobalUncalibRecHit_cfi')
process.load('RecoLocalCalo/EcalRecProducers/ecalRecHit_cfi')
process.ecalRecHit.EBuncalibRecHitCollection = 'ecalGlobalUncalibRecHit:EcalUncalibRecHitsEB'
process.ecalRecHit.EEuncalibRecHitCollection = 'ecalGlobalUncalibRecHit:EcalUncalibRecHitsEE'
process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]




process.demo = cms.EDAnalyzer('CompareAmplitudes')
process.demo.EBURecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEB")
process.demo.EEURecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEE")
process.demo.EBRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB")
process.demo.EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE")


process.demo.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
process.demo.EEDigiCollection = cms.InputTag("ecalDigis","eeDigis")
process.demo.ebEnergyMin = cms.untracked.double(15);  # in GeV
process.demo.eeEnergyMin = cms.untracked.double(15);  # in GeV
process.demo.gncMinEnergy = cms.untracked.double(0.8); # in GeV




process.p = cms.Path(
      process.ecalDigis
      *process.ecalGlobalUncalibRecHit
      *process.ecalDetIdToBeRecovered
      *process.ecalRecHit
      *process.demo
)

