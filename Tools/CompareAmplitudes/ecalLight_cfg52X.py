import FWCore.ParameterSet.Config as cms

process = cms.Process("chi2analysis")

# StandardSequences configurations
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("CalibCalorimetry.EcalLaserCorrection.ecalLaserCorrectionService_cfi")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")


# Global Tag: see https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideFrontierConditions#Global_Tags_for_Global_Run_data
process.GlobalTag.globaltag = 'GR_P_V32::All' 

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.threshold = "WARNING"

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
"/store/data/Run2011B/Photon/RAW/v1/000/176/548/C49FEB43-AEE0-E011-811F-BCAEC5329701.root"

#'file:rawFiveEventsComplete.root'
  )
)
process.source.dropDescendantsOfDroppedBranches = cms.untracked.bool(False)
process.source.inputCommands = cms.untracked.vstring('drop *','keep *FEDRawData*_*_*_*')


process.load("Configuration.StandardSequences.RawToDigi_Data_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")



process.load('RecoLocalCalo/EcalRecProducers/ecalGlobalUncalibRecHit_cfi')
process.load('RecoLocalCalo/EcalRecProducers/ecalRecHit_cfi')
process.ecalRecHit.EBuncalibRecHitCollection = 'ecalGlobalUncalibRecHit:EcalUncalibRecHitsEB'
process.ecalRecHit.EEuncalibRecHitCollection = 'ecalGlobalUncalibRecHit:EcalUncalibRecHitsEE'
#process.ecalRecHit.ChannelStatusToBeExcluded = [ 1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 78, 142 ]




process.demo = cms.EDAnalyzer('CompareAmplitudes')
process.demo.EBURecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEB")
process.demo.EEURecHitCollection = cms.InputTag("ecalGlobalUncalibRecHit","EcalUncalibRecHitsEE")
process.demo.EBRecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEB")
process.demo.EERecHitCollection = cms.InputTag("ecalRecHit","EcalRecHitsEE")


process.demo.EBDigiCollection = cms.InputTag("ecalDigis","ebDigis")
process.demo.EEDigiCollection = cms.InputTag("ecalDigis","eeDigis")
process.demo.ebEnergyMin = cms.untracked.double(15);  # in GeV
#process.demo.ebEnergyMin = cms.untracked.double(-1.e+9);  # in GeV
process.demo.eeEnergyMin = cms.untracked.double(15);  # in GeV
#process.demo.eeEnergyMin = cms.untracked.double(-1.e+9);  # in GeV
process.demo.gncMinEnergy = cms.untracked.double(0.8); # in GeV

'''
process.output = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *'),
#    outputCommands = cms.untracked.vstring('drop *','keep *_*ecal*RecHit*_*_*','keep *FEDRawData*_*_*_*'),
#    outputCommands = cms.untracked.vstring('drop *','keep *_*ecal*UncalibRecHit*_*_*'),
    fileName = cms.untracked.string('edmRECO.root'),
)
'''


process.p = cms.Path(
      process.ecalDigis
      *process.ecalGlobalUncalibRecHit
      *process.ecalDetIdToBeRecovered
      *process.ecalRecHit
      *process.demo
)


#process.e = cms.EndPath(process.output)

