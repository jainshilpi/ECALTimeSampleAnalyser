import FWCore.ParameterSet.Config as cms
#from Configuration.StandardSequences.Eras import eras
#process = cms.Process("EcalTimeSample", eras.Run2_2017)

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('Analyse',Run3)

process.load("FWCore.MessageService.MessageLogger_cfi")
#######################
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.RecoSim_cff')
#process.load('CommonTools.ParticleFlow.EITopPAG_cff')
#process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
#process.load('Configuration.StandardSequences.PATMC_cff')
#process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")

#############################

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v8', '') ###UL 2017 BH 
#process.GlobalTag = GlobalTag(process.GlobalTag, '105X_upgrade2018_realistic_IdealEcalIC_v4', '') ###UL 2017 BH 
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')  
#process.GlobalTag = GlobalTag(process.GlobalTag, '113X_mcRun3_2021_realistic_v10', '')  
process.GlobalTag = GlobalTag(process.GlobalTag, '122X_mcRun3_2021_realistic_v9', '')  

#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(500) )
#process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50) )
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )

#

process.source = cms.Source("PoolSource",
                                # replace 'myfile.root' with the source file you want to use
                                fileNames = cms.untracked.vstring(
                                    #'/store/user/shilpi/V21_AODSIM_BEAM1ON_PU/BeamHalo_2017_Beam1ON/BeamHalo_AODSIM_V21_AODSIM_BEAM1ON_PU/200215_230732/0000/BH_3_1.root'
                                    #'/store/user/shilpi/V21_AODSIM_BEAM1ON_PU/BeamHalo_2017_Beam1ON/BeamHalo_AODSIM_V21_AODSIM_BEAM1ON_PU/200215_230732/0000/BH_3_907.root'
                                    #'/store/mc/RunIIWinter19PFCalibDR/DoubleElectron_FlatPt-1To300/AODSIM/2018ConditionsFlatPU0to70ECALGT_105X_upgrade2018_realistic_IdealEcalIC_v4-v1/40000/0038F7DC-41D2-7A4B-BE93-5F20220400C5.root'
                                    #'/store/mc/Run3Winter20DRMiniAOD/DoubleElectron_FlatPt-1To300/AODSIM/FlatPU0to80_110X_mcRun3_2021_realistic_v6-v3/10000/030B2A09-E305-7A48-AEE6-2BB967D727B9.root'
                                    #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/030B2A09-E305-7A48-AEE6-2BB967D727B9.root'

                                    ### using this till 26th July, 24
                                    'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/045196a3-6768-412f-bef9-96e07364fee6_alcareco.root'
                                    ### trial to run on the mother gensimraw
                                    #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/0ce3d318-6865-476a-8aa0-2fe02f38d6d5.root'

                                    
                                    #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/1fff896e-d966-4b44-8ccc-8d9ae08e36fb.root'
                                    #'file:/store/mc/Run3Winter22DRPremix/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/ALCARECO/EcalUncalZElectron-ALCARECO_122X_mcRun3_2021_realistic_v9_ext1-v3/60000/64ea8a9d-bd0b-45be-9455-7f94b298abe8.root'
                                ),
                            
                            secondaryFileNames = cms.untracked.vstring(

                                ### mother of 045196a3-6768-412f-bef9-96e07364fee6.root
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/0098d95d-c5c2-42ff-ab50-523853ec99ec.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/0ce3d318-6865-476a-8aa0-2fe02f38d6d5.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/13db54b9-f99f-4fba-b09f-f215260a6aef.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/142534c7-0663-4f75-b1e7-2aef2dee099e.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/15c8ea48-75e3-4ca3-a444-196a91ab6452.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/16b7246b-37a7-45db-83c3-db2d09dea9bb.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/1d29b2e1-6390-4bec-b234-2dec450e0dca.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/1ecf8049-7d31-4a99-9a01-c34f10d8b5d2.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/1f6e1262-b126-462b-a860-9e11f3c1a781.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/2edeaa23-53db-4412-a73b-d4625444d103.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/300c77f3-6038-4ac4-b3e7-54049cb8da76.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/35d231f7-b0c1-4b56-8fa0-546c32310965.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/37250065-9c31-4e2b-9a58-2e17f33acb9b.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/3ce937cb-2fbc-4ec9-95cf-7a5dac8dc7c6.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/461937a9-f6d4-44ca-9532-27977e5d930b.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/507aa847-6ae3-4382-87b8-48f294b48cd7.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/5099e05f-b2c3-4173-b89e-13e93891a930.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/53dfd967-9289-4268-aa2e-5d53191f5b6f.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/5512a607-a8d2-4623-bff9-372e491c9f25.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/55d17799-f70e-49cc-929b-089b0c7abe6c.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/560c3659-1dcc-4564-ae59-b5578390a138.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/582b8b83-7d9c-491a-b119-fd23c9b7fead.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/5a238420-90a9-41cc-bd97-dd26d91b26e2.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/5df0cc75-4cc3-4454-b9a0-a79c715a7ab1.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/5eafa783-fe27-4aae-be38-14563ebbb7c3.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/5f872965-d2da-4a43-acd9-2320527977f2.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/62f9b728-6d27-419a-b2d2-3e3363e46e52.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/63262eb1-79e7-49e1-8495-f1321f7d247b.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/691ea470-bae8-49a2-ad75-74ebf34b0474.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/71617ce2-86cc-47fe-be4e-106526ea1fcc.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/77928c0d-5174-4d2b-84b7-b0e7913bbe9e.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/7f0c5be2-484c-42c3-80ad-da4fc5fa3d0a.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/800a4547-4249-4bb4-9001-ff1929c0ec4f.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/8c92c573-9db8-44a4-917d-e610e1b5d0da.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/8cdf8e14-659a-4ac1-917f-a731e0618607.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/8cf02e9a-28e6-4b5d-ac03-31861a9b2b0e.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/8f44b15f-e2d1-40e9-ab3a-b6e6fb15cb1f.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/9156f536-9da8-4b53-895c-5d098cfb666e.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/9982b7e8-8d99-4927-ae23-10182552a524.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/9bdcfc77-0654-4aa5-bfc7-6baace79f9fe.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/9ea64793-a05e-4aae-81af-c53cfd11a317.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/a5c8b0c7-5f15-439e-86fa-d0c70389e2d1.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/b00043f1-f6f5-4095-899e-29011bb8de93.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/b08d02a2-75fd-4b48-9bae-57d383f51671.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/cda239e4-2bf4-4df5-89b7-f4ed79e269ae.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/cdd7098c-280d-4ae2-ac90-985f7699b7a1.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/da4d4862-07a4-4775-a83d-b4eb15553f12.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/ddd8db4c-5d6d-4f39-89b3-219652e7231b.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/e0bc49b3-2932-4695-a29c-09cdc0eab781.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/e2e1469a-d695-4a02-8536-4f4459759dd2.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/e682c8cb-eab5-434f-8536-f0216029a2a0.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/e716eff7-4549-4790-b703-c493b2b80f63.root',
                                #'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/ea8b8650-be12-4f88-9b92-8f4a22e02aed.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/ef0dc043-7a81-464b-89e0-1f882c69cd72.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/efa24374-a5dd-4c05-b4f6-56d8098439f9.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/f7dcc73d-399c-4f7a-a63f-78273a013ef0.root',
                                'file:/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/fe37489b-f8f2-426e-9c06-683dd74a930e.root'
                            )
                            

                        )


process.InterimOutput = cms.OutputModule("PoolOutputModule",
                                         fileName = cms.untracked.string('myOutputFile.root'), 
                                          SelectEvents = cms.untracked.PSet(
                                              SelectEvents = cms.vstring("p")
                                              ),
                                         outputCommands = cms.untracked.vstring('keep *')
                                         
                                     )

                                                        


process.TFileService = cms.Service("TFileService", fileName = cms.string('/tmp/shilpi/timeSampleTree.root'))

###https://cmssdt.cern.ch/lxr/source/RecoEcal/EgammaClusterProducers/python/ecalDigiSelector_cfi.py
process.timeSample = cms.EDAnalyzer('ECALTimeSampleAnalyser',
                                    #EBdigiCollection = cms.InputTag("selectDigi","selectedEcalEBDigiCollection"),
                                    #EEdigiCollection = cms.InputTag("selectDigi","selectedEcalEEDigiCollection"),
                                    #ebRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
                                    #eeRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
                                    
                                    simHitEBCollection = cms.InputTag("g4SimHits","EcalHitsEB"),
                                    simHitEECollection = cms.InputTag("g4SimHits","EcalHitsEE"),
                                    #EBdigiCollection = cms.InputTag("ecalDigis","ebDigis", "RECO"),
                                    #EEdigiCollection = cms.InputTag("ecalDigis","eeDigis", "RECO"),
                                    EBdigiCollection = cms.InputTag("ecalDigis","ebDigis", "Analyse"),
                                    EEdigiCollection = cms.InputTag("ecalDigis","eeDigis", "Analyse"),

                                    ebRecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEB"),
                                    eeRecHitCollection = cms.InputTag("ecalRecHit", "EcalRecHitsEE"),
                                    #electronSrc          = cms.InputTag("selectedPatElectrons"),
                                    #photonSrc            = cms.InputTag("selectedPatPhotons")
                                    #recoEleSrc            = cms.InputTag("gedGsfElectrons") ,
                                    recoEleSrc            = cms.InputTag("electronRecalibSCAssociator") ,
                                    genParticleSrc        = cms.InputTag("genParticles"),
                                    #recoPhotonSrc            = cms.InputTag("gedPhotons") ,
                                    runMinBias            = cms.bool(False)
)

process.load('RecoLuminosity.LumiProducer.bunchSpacingProducer_cfi')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load("RecoEcal.EgammaClusterProducers.reducedRecHitsSequence_cff")
#process.load("RecoLocalCalo.EcalRecProducers.ecalLocalRecoSequence_cff")
process.load("RecoLocalCalo.Configuration.ecalLocalRecoSequence_cff")
#process.load("RecoEcal.EgammaClusterProducers.particleFlowSuperClusterECAL_cfi")
process.load("RecoEcal.Configuration.RecoEcal_cff")
process.load("Calibration.EcalCalibAlgos.electronRecalibSCAssociator_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterECALUncorrected_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowClusterPS_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi")
process.load("RecoParticleFlow.PFClusterProducer.particleFlowRecHitPS_cfi")

###electron-SC matching needs to be run since I am re-running the SC collection and hence when ele-->SC, it doesnot recognize. So I recreate another electron collection. SC algo needed to be run because it didnt exist and needed for ele-->SC. 
###https://cmssdt.cern.ch/lxr/source/Calibration/EcalCalibAlgos/src/ElectronRecalibSuperClusterAssociator.cc
### https://cmssdt.cern.ch/lxr/source/Calibration/EcalCalibAlgos/python/electronRecalibSCAssociator_cfi.py

### inputs needed for generating EBSr and EESr
#process.load('L1Trigger.Configuration.SimL1Emulator_cff')
#process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

### https://cmssdt.cern.ch/lxr/source/SimCalorimetry/Configuration/python/ecalDigiSequenceComplete_cff.py#0008
# Selective Readout Processor producer
#process.load("SimCalorimetry.EcalSelectiveReadoutProducers.ecalDigis_cff")

process.load('Configuration.StandardSequences.Digi_cff')

#process.ecalDigis.cpu.forceToKeepFRData = cms.bool(True)


process.p = cms.Path(
    process.bunchSpacingProducer *
    process.ecalDigis *
    process.ecalPreshowerDigis *
    process.ecalPreshowerRecHit *
    process.ecalLocalRecoSequence *
    process.ecalClusters *
    process.particleFlowRecHitPS *
    process.particleFlowRecHitECAL *
    process.particleFlowClusterECALUncorrected *
    process.particleFlowClusterPS *
    process.particleFlowClusterECAL *
    process.electronRecalibSCAssociator 
    *process.timeSample
)

#process.e = cms.EndPath(process.InterimOutput)  


print(process.dumpPython())
