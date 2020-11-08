import FWCore.ParameterSet.Config as cms

process = cms.Process("TEST")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('/store/mc/RunIIFall17MiniAODv2/WGToLNuG_01J_5f_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v3/30000/002D72F8-662F-E911-A6B3-A4BADB1E67B8.root'),
    secondaryFileNames = cms.untracked.vstring()
)
process.ChargeSignificanceTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('ChargeSignificanceTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0)
)

process.CkfBaseTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.CkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('CkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    estimator = cms.string('Chi2'),
    intermediateCleaning = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    updator = cms.string('KFUpdator')
)

process.CompositeTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('CompositeTrajectoryFilter'),
    filters = cms.VPSet()
)

process.GroupedCkfTrajectoryBuilder = cms.PSet(
    ComponentType = cms.string('GroupedCkfTrajectoryBuilder'),
    MeasurementTrackerName = cms.string(''),
    TTRHBuilder = cms.string('WithTrackAngle'),
    alwaysUseInvalidHits = cms.bool(True),
    bestHitOnly = cms.bool(True),
    estimator = cms.string('Chi2'),
    foundHitBonus = cms.double(10.0),
    inOutTrajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    intermediateCleaning = cms.bool(True),
    keepOriginalIfRebuildFails = cms.bool(False),
    lockHits = cms.bool(True),
    lostHitPenalty = cms.double(30.0),
    maxCand = cms.int32(5),
    minNrOfHitsForRebuild = cms.int32(5),
    propagatorAlong = cms.string('PropagatorWithMaterial'),
    propagatorOpposite = cms.string('PropagatorWithMaterialOpposite'),
    requireSeedHitsInRebuild = cms.bool(True),
    trajectoryFilter = cms.PSet(
        refToPSet_ = cms.string('CkfBaseTrajectoryFilter_block')
    ),
    updator = cms.string('KFUpdator'),
    useSameTrajFilter = cms.bool(True)
)

process.HFRecalParameterBlock = cms.PSet(
    HFdepthOneParameterA = cms.vdouble(
        0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
        0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
        0.058939, 0.125497
    ),
    HFdepthOneParameterB = cms.vdouble(
        -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
        2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
        0.000425, 0.000209
    ),
    HFdepthTwoParameterA = cms.vdouble(
        0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
        0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
        0.051579, 0.086593
    ),
    HFdepthTwoParameterB = cms.vdouble(
        -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
        1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
        0.000157, -3e-06
    )
)

process.MaxCCCLostHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxCCCLostHitsTrajectoryFilter'),
    maxCCCLostHits = cms.int32(3),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutLoose')
    )
)

process.MaxConsecLostHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxConsecLostHitsTrajectoryFilter'),
    maxConsecLostHits = cms.int32(1)
)

process.MaxHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxHitsTrajectoryFilter'),
    maxNumberOfHits = cms.int32(100)
)

process.MaxLostHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MaxLostHitsTrajectoryFilter'),
    maxLostHits = cms.int32(2)
)

process.MinHitsTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MinHitsTrajectoryFilter'),
    minimumNumberOfHits = cms.int32(5)
)

process.MinPtTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('MinPtTrajectoryFilter'),
    minHitsMinPt = cms.int32(3),
    minPt = cms.double(1.0),
    nSigmaMinPt = cms.double(5.0)
)

process.SiStripClusterChargeCutLoose = cms.PSet(
    value = cms.double(1620.0)
)

process.SiStripClusterChargeCutNone = cms.PSet(
    value = cms.double(-1.0)
)

process.SiStripClusterChargeCutTight = cms.PSet(
    value = cms.double(1945.0)
)

process.SiStripClusterChargeCutTiny = cms.PSet(
    value = cms.double(800.0)
)

process.ThresholdPtTrajectoryFilter_block = cms.PSet(
    ComponentType = cms.string('ThresholdPtTrajectoryFilter'),
    minHitsThresholdPt = cms.int32(3),
    nSigmaThresholdPt = cms.double(5.0),
    thresholdPt = cms.double(10.0)
)

process.calibratedEgammaPatSettings = cms.PSet(
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    semiDeterministic = cms.bool(True)
)

process.calibratedEgammaSettings = cms.PSet(
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    semiDeterministic = cms.bool(True)
)

process.ckfBaseInOutTrajectoryFilter = cms.PSet(
    ComponentType = cms.string('CkfBaseTrajectoryFilter'),
    chargeSignificance = cms.double(-1.0),
    constantValueForLostHitsFractionFilter = cms.double(2.0),
    extraNumberOfHitsBeforeTheFirstLoop = cms.int32(4),
    maxCCCLostHits = cms.int32(9999),
    maxConsecLostHits = cms.int32(1),
    maxLostHits = cms.int32(999),
    maxLostHitsFraction = cms.double(0.1),
    maxNumberOfHits = cms.int32(100),
    minGoodStripCharge = cms.PSet(
        refToPSet_ = cms.string('SiStripClusterChargeCutNone')
    ),
    minHitsMinPt = cms.int32(3),
    minNumberOfHitsForLoopers = cms.int32(13),
    minNumberOfHitsPerLoop = cms.int32(4),
    minPt = cms.double(0.9),
    minimumNumberOfHits = cms.int32(5),
    nSigmaMinPt = cms.double(5.0),
    pixelSeedExtension = cms.bool(False),
    seedExtension = cms.int32(0),
    seedPairPenalty = cms.int32(0),
    strictSeedExtension = cms.bool(False)
)

process.ecalTrkCombinationRegression = cms.PSet(
    ecalTrkRegressionConfig = cms.PSet(
        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
        eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
        forceHighEnergyTrainingIfSaturated = cms.bool(False),
        lowEtHighEtBoundary = cms.double(50.0),
        rangeMax = cms.double(3.0),
        rangeMin = cms.double(-1.0)
    ),
    ecalTrkRegressionUncertConfig = cms.PSet(
        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
        forceHighEnergyTrainingIfSaturated = cms.bool(False),
        lowEtHighEtBoundary = cms.double(50.0),
        rangeMax = cms.double(0.5),
        rangeMin = cms.double(0.0002)
    ),
    maxEPDiffInSigmaForComb = cms.double(15.0),
    maxEcalEnergyForComb = cms.double(200.0),
    maxRelTrkMomErrForComb = cms.double(10.0),
    minEOverPForComb = cms.double(0.025)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(30000)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    wantSummary = cms.untracked.bool(True)
)

process.pfJetIDSelector = cms.PSet(
    quality = cms.string('TIGHT'),
    version = cms.string('WINTER17')
)

process.JetUserData = cms.EDProducer("JetUserData",
    candSVTagInfos = cms.string('pfInclusiveSecondaryVertexFinder'),
    coneSize = cms.double(0.4),
    getJERFromTxt = cms.bool(False),
    hlt2reco_deltaRmax = cms.double(0.2),
    hltJetFilter = cms.InputTag("hltPFHT"),
    hltPath = cms.string('HLT_PFHT800'),
    jecAK4chsPayloadNames_jetUserdata = cms.vstring(
        'Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt', 
        'Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt', 
        'Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt'
    ),
    jerLabel = cms.string('AK4PFchs'),
    jetCorrLabel = cms.string('AK4PFchs'),
    jetLabel = cms.InputTag("slimmedJets"),
    resolutionsFile = cms.string('Fall17_17Nov2017_V32_MC_PtResolution_AK4PFchs.txt'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    scaleFactorsFile = cms.string('Fall17_17Nov2017_V32_MC_SF_AK4PFchs.txt'),
    triggerResults = cms.InputTag("TriggerResults","","HLT"),
    triggerSummary = cms.InputTag("hltTriggerSummaryAOD","","HLT"),
    vertex_jetUserdata = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.Wtoenu = cms.EDProducer("PKUWLepProducer",
    MET = cms.InputTag("slimmedMETs"),
    cut = cms.string(''),
    leptons = cms.InputTag("goodElectrons")
)


process.Wtomunu = cms.EDProducer("PKUWLepProducer",
    MET = cms.InputTag("slimmedMETs"),
    cut = cms.string(''),
    leptons = cms.InputTag("goodMuons")
)


process.calibratedElectrons = cms.EDProducer("CalibratedElectronProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    epCombConfig = cms.PSet(
        ecalTrkRegressionConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMax = cms.double(3.0),
            rangeMin = cms.double(-1.0)
        ),
        ecalTrkRegressionUncertConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMax = cms.double(0.5),
            rangeMin = cms.double(0.0002)
        ),
        maxEPDiffInSigmaForComb = cms.double(15.0),
        maxEcalEnergyForComb = cms.double(200.0),
        maxRelTrkMomErrForComb = cms.double(10.0),
        minEOverPForComb = cms.double(0.025)
    ),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("gedGsfElectrons")
)


process.calibratedPatElectrons = cms.EDProducer("CalibratedPatElectronProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    epCombConfig = cms.PSet(
        ecalTrkRegressionConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMax = cms.double(3.0),
            rangeMin = cms.double(-1.0)
        ),
        ecalTrkRegressionUncertConfig = cms.PSet(
            ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
            ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
            eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
            eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
            forceHighEnergyTrainingIfSaturated = cms.bool(False),
            lowEtHighEtBoundary = cms.double(50.0),
            rangeMax = cms.double(0.5),
            rangeMin = cms.double(0.0002)
        ),
        maxEPDiffInSigmaForComb = cms.double(15.0),
        maxEcalEnergyForComb = cms.double(200.0),
        maxRelTrkMomErrForComb = cms.double(10.0),
        minEOverPForComb = cms.double(0.025)
    ),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(False),
    recHitCollectionEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("slimmedElectrons","","@skipCurrentProcess")
)


process.calibratedPatPhotons = cms.EDProducer("CalibratedPatPhotonProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(False),
    recHitCollectionEB = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    recHitCollectionEE = cms.InputTag("reducedEgamma","reducedEERecHits"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
)


process.calibratedPhotons = cms.EDProducer("CalibratedPhotonProducer",
    correctionFile = cms.string('EgammaAnalysis/ElectronTools/data/ScalesSmearings/Run2017_17Nov2017_v1_ele_unc'),
    minEtToCalibrate = cms.double(5.0),
    produceCalibratedObjs = cms.bool(True),
    recHitCollectionEB = cms.InputTag("reducedEcalRecHitsEB"),
    recHitCollectionEE = cms.InputTag("reducedEcalRecHitsEE"),
    semiDeterministic = cms.bool(True),
    src = cms.InputTag("gedPhotons")
)


process.ckfTrackCandidates = cms.EDProducer("CkfTrackCandidateMaker",
    MeasurementTrackerEvent = cms.InputTag("MeasurementTrackerEvent"),
    NavigationSchool = cms.string('SimpleNavigationSchool'),
    RedundantSeedCleaner = cms.string('CachingSeedCleanerBySharedInput'),
    SimpleMagneticField = cms.string(''),
    TrajectoryBuilder = cms.string('GroupedCkfTrajectoryBuilder'),
    TrajectoryBuilderPSet = cms.PSet(
        refToPSet_ = cms.string('GroupedCkfTrajectoryBuilder')
    ),
    TrajectoryCleaner = cms.string('TrajectoryCleanerBySharedHits'),
    TransientInitialStateEstimatorParameters = cms.PSet(
        numberMeasurementsForFit = cms.int32(4),
        propagatorAlongTISE = cms.string('PropagatorWithMaterial'),
        propagatorOppositeTISE = cms.string('PropagatorWithMaterialOpposite')
    ),
    cleanTrajectoryAfterInOut = cms.bool(True),
    doSeedingRegionRebuilding = cms.bool(True),
    maxNSeeds = cms.uint32(500000),
    maxSeedsBeforeCleaning = cms.uint32(5000),
    reverseTrajectories = cms.bool(False),
    src = cms.InputTag("globalMixedSeeds"),
    useHitsSplitting = cms.bool(True)
)


process.cleanAK4Jets = cms.EDProducer("PATJetCleaner",
    checkOverlaps = cms.PSet(
        electrons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.4),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(True),
            src = cms.InputTag("goodElectrons")
        ),
        muons = cms.PSet(
            algorithm = cms.string('byDeltaR'),
            checkRecoComponents = cms.bool(False),
            deltaR = cms.double(0.4),
            pairCut = cms.string(''),
            preselection = cms.string(''),
            requireNoOverlaps = cms.bool(True),
            src = cms.InputTag("goodMuons")
        ),
        photons = cms.PSet(

        ),
        taus = cms.PSet(

        ),
        tkIsoElectrons = cms.PSet(

        )
    ),
    finalCut = cms.string('pt > 20 && abs(eta) < 4.7'),
    preselection = cms.string(''),
    src = cms.InputTag("goodAK4Jets")
)


process.goodElectrons = cms.EDProducer("PATElectronIdSelector",
    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'),
    idLabel = cms.string('tight'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedElectrons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.goodMuons = cms.EDProducer("PATMuonIdSelector",
    idLabel = cms.string('tight'),
    src = cms.InputTag("slimmedMuons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.leptonicV = cms.EDProducer("CandViewMerger",
    cut = cms.string(''),
    src = cms.VInputTag("Wtoenu", "Wtomunu")
)


process.looseMuons = cms.EDProducer("PATMuonIdSelector",
    idLabel = cms.string('loose'),
    src = cms.InputTag("slimmedMuons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.photonIDValueMapProducer = cms.EDProducer("PhotonIDValueMapProducer",
    ebReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    ebReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    eeReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
    eeReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedEERecHits"),
    esReducedRecHitCollection = cms.InputTag("reducedEcalRecHitsES"),
    esReducedRecHitCollectionMiniAOD = cms.InputTag("reducedEgamma","reducedESRecHits"),
    particleBasedIsolation = cms.InputTag("particleBasedIsolation","gedPhotons"),
    pfCandidates = cms.InputTag("particleFlow"),
    pfCandidatesMiniAOD = cms.InputTag("packedPFCandidates"),
    src = cms.InputTag("gedPhotons"),
    srcMiniAOD = cms.InputTag("slimmedPhotons","","@skipCurrentProcess"),
    vertices = cms.InputTag("offlinePrimaryVertices"),
    verticesMiniAOD = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.prefiringweight = cms.EDProducer("L1ECALPrefiringWeightProducer",
    DataEra = cms.string('2017BtoF'),
    L1Maps = cms.string('L1PrefiringMaps.root'),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = cms.bool(False),
    TheJets = cms.InputTag("slimmedJets"),
    ThePhotons = cms.InputTag("slimmedPhotons"),
    UseJetEMPt = cms.bool(False)
)


process.slimmedElectrons = cms.EDProducer("ModifiedElectronProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(
            cms.PSet(
                electron_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPreCorr"),
                    ecalTrkEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPostCorr"),
                    ecalTrkEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPreCorr"),
                    ecalTrkEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPostCorr"),
                    ecalTrkEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPreCorr"),
                    electronSrc = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
                    energyScaleDown = cms.InputTag("calibratedPatElectrons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatElectrons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatElectrons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatElectrons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatElectrons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatElectrons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatElectrons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatElectrons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatElectrons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatElectrons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatElectrons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatElectrons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatElectrons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatElectrons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatElectrons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatElectrons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatElectrons","energySmearNrSigma")
                ),
                modifierName = cms.string('EGExtraInfoModifierFromFloatValueMaps'),
                overrideExistingValues = cms.bool(True),
                photon_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPreCorr"),
                    energyScaleDown = cms.InputTag("calibratedPatPhotons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatPhotons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatPhotons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatPhotons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatPhotons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatPhotons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatPhotons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatPhotons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatPhotons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatPhotons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatPhotons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatPhotons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatPhotons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatPhotons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatPhotons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatPhotons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatPhotons","energySmearNrSigma"),
                    photonSrc = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
                )
            ), 
            cms.PSet(
                epCombConfig = cms.PSet(
                    ecalTrkRegressionConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMax = cms.double(3.0),
                        rangeMin = cms.double(-1.0)
                    ),
                    ecalTrkRegressionUncertConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMax = cms.double(0.5),
                        rangeMin = cms.double(0.0002)
                    ),
                    maxEPDiffInSigmaForComb = cms.double(15.0),
                    maxEcalEnergyForComb = cms.double(200.0),
                    maxRelTrkMomErrForComb = cms.double(10.0),
                    minEOverPForComb = cms.double(0.025)
                ),
                modifierName = cms.string('EGEtScaleSysModifier'),
                overrideExistingValues = cms.bool(True),
                uncertFunc = cms.PSet(
                    highEt = cms.double(46.5),
                    highEtUncert = cms.double(-0.002),
                    lowEt = cms.double(43.5),
                    lowEtUncert = cms.double(0.002),
                    name = cms.string('UncertFuncV1')
                )
            )
        )
    ),
    src = cms.InputTag("slimmedElectrons","","@skipCurrentProcess")
)


process.slimmedPhotons = cms.EDProducer("ModifiedPhotonProducer",
    modifierConfig = cms.PSet(
        modifications = cms.VPSet(
            cms.PSet(
                electron_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalEnergyPreCorr"),
                    ecalTrkEnergyErrPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPostCorr"),
                    ecalTrkEnergyErrPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyErrPreCorr"),
                    ecalTrkEnergyPostCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPostCorr"),
                    ecalTrkEnergyPreCorr = cms.InputTag("calibratedPatElectrons","ecalTrkEnergyPreCorr"),
                    electronSrc = cms.InputTag("slimmedElectrons","","@skipCurrentProcess"),
                    energyScaleDown = cms.InputTag("calibratedPatElectrons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatElectrons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatElectrons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatElectrons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatElectrons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatElectrons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatElectrons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatElectrons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatElectrons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatElectrons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatElectrons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatElectrons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatElectrons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatElectrons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatElectrons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatElectrons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatElectrons","energySmearNrSigma")
                ),
                modifierName = cms.string('EGExtraInfoModifierFromFloatValueMaps'),
                overrideExistingValues = cms.bool(True),
                photon_config = cms.PSet(
                    ecalEnergyErrPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPostCorr"),
                    ecalEnergyErrPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyErrPreCorr"),
                    ecalEnergyPostCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPostCorr"),
                    ecalEnergyPreCorr = cms.InputTag("calibratedPatPhotons","ecalEnergyPreCorr"),
                    energyScaleDown = cms.InputTag("calibratedPatPhotons","energyScaleDown"),
                    energyScaleGainDown = cms.InputTag("calibratedPatPhotons","energyScaleGainDown"),
                    energyScaleGainUp = cms.InputTag("calibratedPatPhotons","energyScaleGainUp"),
                    energyScaleStatDown = cms.InputTag("calibratedPatPhotons","energyScaleStatDown"),
                    energyScaleStatUp = cms.InputTag("calibratedPatPhotons","energyScaleStatUp"),
                    energyScaleSystDown = cms.InputTag("calibratedPatPhotons","energyScaleSystDown"),
                    energyScaleSystUp = cms.InputTag("calibratedPatPhotons","energyScaleSystUp"),
                    energyScaleUp = cms.InputTag("calibratedPatPhotons","energyScaleUp"),
                    energyScaleValue = cms.InputTag("calibratedPatPhotons","energyScaleValue"),
                    energySigmaDown = cms.InputTag("calibratedPatPhotons","energySigmaDown"),
                    energySigmaPhiDown = cms.InputTag("calibratedPatPhotons","energySigmaPhiDown"),
                    energySigmaPhiUp = cms.InputTag("calibratedPatPhotons","energySigmaPhiUp"),
                    energySigmaRhoDown = cms.InputTag("calibratedPatPhotons","energySigmaRhoDown"),
                    energySigmaRhoUp = cms.InputTag("calibratedPatPhotons","energySigmaRhoUp"),
                    energySigmaUp = cms.InputTag("calibratedPatPhotons","energySigmaUp"),
                    energySigmaValue = cms.InputTag("calibratedPatPhotons","energySigmaValue"),
                    energySmearNrSigma = cms.InputTag("calibratedPatPhotons","energySmearNrSigma"),
                    photonSrc = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
                )
            ), 
            cms.PSet(
                epCombConfig = cms.PSet(
                    ecalTrkRegressionConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMax = cms.double(3.0),
                        rangeMin = cms.double(-1.0)
                    ),
                    ecalTrkRegressionUncertConfig = cms.PSet(
                        ebHighEtForestName = cms.string('electron_eb_ECALTRK_var'),
                        ebLowEtForestName = cms.string('electron_eb_ECALTRK_lowpt_var'),
                        eeHighEtForestName = cms.string('electron_ee_ECALTRK_var'),
                        eeLowEtForestName = cms.string('electron_ee_ECALTRK_lowpt_var'),
                        forceHighEnergyTrainingIfSaturated = cms.bool(False),
                        lowEtHighEtBoundary = cms.double(50.0),
                        rangeMax = cms.double(0.5),
                        rangeMin = cms.double(0.0002)
                    ),
                    maxEPDiffInSigmaForComb = cms.double(15.0),
                    maxEcalEnergyForComb = cms.double(200.0),
                    maxRelTrkMomErrForComb = cms.double(10.0),
                    minEOverPForComb = cms.double(0.025)
                ),
                modifierName = cms.string('EGEtScaleSysModifier'),
                overrideExistingValues = cms.bool(True),
                uncertFunc = cms.PSet(
                    highEt = cms.double(46.5),
                    highEtUncert = cms.double(-0.002),
                    lowEt = cms.double(43.5),
                    lowEtUncert = cms.double(0.002),
                    name = cms.string('UncertFuncV1')
                )
            )
        )
    ),
    src = cms.InputTag("slimmedPhotons","","@skipCurrentProcess")
)


process.vetoElectrons = cms.EDProducer("PATElectronIdSelector",
    effAreasConfigFile = cms.FileInPath('RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'),
    idLabel = cms.string('veto'),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    src = cms.InputTag("slimmedElectrons"),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.BadChargedCandidateFilter = cms.EDFilter("BadParticleFilter",
    PFCandidates = cms.InputTag("packedPFCandidates"),
    filterType = cms.string('BadChargedCandidate'),
    innerTrackRelErr = cms.double(1.0),
    maxDR = cms.double(1e-05),
    minMuonPt = cms.double(100.0),
    minMuonTrackRelErr = cms.double(2.0),
    minPtDiffRel = cms.double(1e-05),
    muons = cms.InputTag("slimmedMuons"),
    segmentCompatibility = cms.double(0.3),
    taggingMode = cms.bool(False)
)


process.BadPFMuonFilter = cms.EDFilter("BadParticleFilter",
    PFCandidates = cms.InputTag("packedPFCandidates"),
    algo = cms.int32(14),
    filterType = cms.string('BadPFMuon'),
    innerTrackRelErr = cms.double(1.0),
    maxDR = cms.double(0.001),
    minMuonPt = cms.double(100),
    minMuonTrackRelErr = cms.double(2.0),
    minPtDiffRel = cms.double(0.0),
    muons = cms.InputTag("slimmedMuons"),
    segmentCompatibility = cms.double(0.3),
    taggingMode = cms.bool(False)
)


process.ecalBadCalibFilter = cms.EDFilter("EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEcalRecHitsEE"),
    baddetEcal = cms.vuint32(),
    debug = cms.bool(False),
    ecalMinEt = cms.double(50.0),
    taggingMode = cms.bool(False)
)


process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter("EcalBadCalibFilter",
    EcalRecHitSource = cms.InputTag("reducedEgamma","reducedEERecHits"),
    baddetEcal = cms.vuint32(
        872439604, 872422825, 872420274, 872423218, 872423215, 
        872416066, 872435036, 872439336, 872420273, 872436907, 
        872420147, 872439731, 872436657, 872420397, 872439732, 
        872439339, 872439603, 872422436, 872439861, 872437051, 
        872437052, 872420649, 872422436, 872421950, 872437185, 
        872422564, 872421566, 872421695, 872421955, 872421567, 
        872437184, 872421951, 872421694, 872437056, 872437057, 
        872437313
    ),
    debug = cms.bool(False),
    ecalMinEt = cms.double(50.0),
    taggingMode = cms.bool(True)
)


process.goodAK4Jets = cms.EDFilter("PFJetIDSelectionFunctorFilter",
    filterParams = cms.PSet(
        quality = cms.string('TIGHT'),
        version = cms.string('WINTER17')
    ),
    src = cms.InputTag("JetUserData")
)


process.goodPhotons = cms.EDFilter("PATPhotonSelector",
    cut = cms.string('pt > 15 && fabs(eta) < 2.5'),
    src = cms.InputTag("slimmedPhotons")
)


process.leptonicVFilter = cms.EDFilter("CandViewCountFilter",
    minNumber = cms.uint32(0),
    src = cms.InputTag("leptonicV")
)


process.leptonicVSelector = cms.EDFilter("CandViewSelector",
    cut = cms.string('pt > 0.0'),
    filter = cms.bool(False),
    src = cms.InputTag("leptonicV")
)


process.treeDumper = cms.EDAnalyzer("PKUTreeMaker",
    PKUChannel = cms.string('VW_CHANNEL'),
    RunOnMC = cms.bool(True),
    ak4jetsSrc = cms.InputTag("cleanAK4Jets"),
    beamSpot = cms.InputTag("offlineBeamSpot","","RECO"),
    conversions = cms.InputTag("reducedEgamma","reducedConversions","PAT"),
    crossSectionPb = cms.double(1),
    effAreaChHadFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt'),
    effAreaNeuHadFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfNeutralHadrons_90percentBased.txt'),
    effAreaPhoFile = cms.FileInPath('RecoEgamma/PhotonIdentification/data/Fall17/effAreaPhotons_cone03_pfPhotons_90percentBased.txt'),
    elPaths1 = cms.vstring('HLT_Ele23_WPTight_Gsf_v*'),
    elPaths2 = cms.vstring('HLT_Ele32 WPTight Gsf L1DoubleEG_v*'),
    electrons = cms.InputTag("slimmedElectrons"),
    full5x5SigmaIEtaIEtaMap = cms.InputTag("photonIDValueMapProducer","phoFull5x5SigmaIEtaIEta"),
    genJet = cms.InputTag("slimmedGenJets"),
    genSrc = cms.InputTag("prunedGenParticles"),
    generator = cms.InputTag("generator"),
    goodeleSrc = cms.InputTag("goodElectrons"),
    goodmuonSrc = cms.InputTag("goodMuons"),
    hltToken = cms.InputTag("TriggerResults","","HLT"),
    isGen = cms.bool(False),
    jecAK4PayloadNames = cms.vstring(
        'Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt', 
        'Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt', 
        'Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt'
    ),
    jecAK4chsPayloadNames = cms.vstring(
        'Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs.txt', 
        'Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs.txt', 
        'Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs.txt'
    ),
    leptonicVSrc = cms.InputTag("leptonicV"),
    lhe = cms.InputTag("externalLHEProducer"),
    looseelectronSrc = cms.InputTag("vetoElectrons"),
    loosemuonSrc = cms.InputTag("looseMuons"),
    metSrc = cms.InputTag("slimmedMETs"),
    muPaths1 = cms.vstring(
        'HLT_IsoMu20_v*', 
        'HLT_IsoTkMu20_v*'
    ),
    muPaths2 = cms.vstring(
        'HLT_IsoMu24_v*', 
        'HLT_IsoTkMu24_v*'
    ),
    muPaths3 = cms.vstring('HLT_IsoMu27_v*'),
    noiseFilter = cms.InputTag("TriggerResults","","PAT"),
    noiseFilterSelection_EcalDeadCellTriggerPrimitiveFilter = cms.string('Flag_EcalDeadCellTriggerPrimitiveFilter'),
    noiseFilterSelection_HBHENoiseFilter = cms.string('Flag_HBHENoiseFilter'),
    noiseFilterSelection_HBHENoiseIsoFilter = cms.string('Flag_HBHENoiseIsoFilter'),
    noiseFilterSelection_badChargedHadron = cms.InputTag("BadChargedCandidateFilter"),
    noiseFilterSelection_badMuon = cms.InputTag("BadPFMuonFilter"),
    noiseFilterSelection_eeBadScFilter = cms.string('Flag_eeBadScFilter'),
    noiseFilterSelection_globalTightHaloFilter = cms.string('Flag_globalTightHalo2016Filter'),
    noiseFilterSelection_goodVertices = cms.string('Flag_goodVertices'),
    originalNEvents = cms.int32(1),
    phoChargedIsolation = cms.InputTag("photonIDValueMapProducer","phoChargedIsolation"),
    phoNeutralHadronIsolation = cms.InputTag("photonIDValueMapProducer","phoNeutralHadronIsolation"),
    phoPhotonIsolation = cms.InputTag("photonIDValueMapProducer","phoPhotonIsolation"),
    photonSrc = cms.InputTag("slimmedPhotons"),
    pileup = cms.InputTag("slimmedAddPileupInfo"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    t1jetSrc = cms.InputTag("slimmedJets"),
    t1jetSrc_user = cms.InputTag("JetUserData"),
    t1muSrc = cms.InputTag("slimmedMuons"),
    targetLumiInvPb = cms.double(1.0),
    vertex = cms.InputTag("offlineSlimmedPrimaryVertices")
)


process.MessageLogger = cms.Service("MessageLogger",
    FrameworkJobReport = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        optionalPSet = cms.untracked.bool(True)
    ),
    categories = cms.untracked.vstring(
        'FwkJob', 
        'FwkReport', 
        'FwkSummary', 
        'Root_NoDictionary'
    ),
    cerr = cms.untracked.PSet(
        FwkJob = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        FwkReport = cms.untracked.PSet(
            limit = cms.untracked.int32(99999999),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(10)
        ),
        FwkSummary = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000),
            optionalPSet = cms.untracked.bool(True),
            reportEvery = cms.untracked.int32(1)
        ),
        INFO = cms.untracked.PSet(
            limit = cms.untracked.int32(0)
        ),
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        default = cms.untracked.PSet(
            limit = cms.untracked.int32(10000000)
        ),
        noTimeStamps = cms.untracked.bool(False),
        optionalPSet = cms.untracked.bool(True),
        threshold = cms.untracked.string('INFO')
    ),
    cerr_stats = cms.untracked.PSet(
        optionalPSet = cms.untracked.bool(True),
        output = cms.untracked.string('cerr'),
        threshold = cms.untracked.string('WARNING')
    ),
    cout = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring(),
    debugs = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    default = cms.untracked.PSet(

    ),
    destinations = cms.untracked.vstring(
        'warnings', 
        'errors', 
        'infos', 
        'debugs', 
        'cout', 
        'cerr'
    ),
    errors = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    ),
    fwkJobReports = cms.untracked.vstring('FrameworkJobReport'),
    infos = cms.untracked.PSet(
        Root_NoDictionary = cms.untracked.PSet(
            limit = cms.untracked.int32(0),
            optionalPSet = cms.untracked.bool(True)
        ),
        optionalPSet = cms.untracked.bool(True),
        placeholder = cms.untracked.bool(True)
    ),
    statistics = cms.untracked.vstring('cerr_stats'),
    suppressDebug = cms.untracked.vstring(),
    suppressInfo = cms.untracked.vstring(),
    suppressWarning = cms.untracked.vstring(),
    warnings = cms.untracked.PSet(
        placeholder = cms.untracked.bool(True)
    )
)


process.TFileService = cms.Service("TFileService",
    fileName = cms.string('treePKU.root')
)


process.CSCGeometryESModule = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.CaloGeometryBuilder = cms.ESProducer("CaloGeometryBuilder",
    SelectedCalos = cms.vstring(
        'HCAL', 
        'ZDC', 
        'CASTOR', 
        'EcalBarrel', 
        'EcalEndcap', 
        'EcalPreshower', 
        'TOWER'
    )
)


process.CaloTopologyBuilder = cms.ESProducer("CaloTopologyBuilder")


process.CaloTowerGeometryFromDBEP = cms.ESProducer("CaloTowerGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.CaloTowerTopologyEP = cms.ESProducer("CaloTowerTopologyEP")


process.CastorDbProducer = cms.ESProducer("CastorDbProducer",
    appendToDataLabel = cms.string('')
)


process.CastorGeometryFromDBEP = cms.ESProducer("CastorGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.Chi2MeasurementEstimator = cms.ESProducer("Chi2MeasurementEstimatorESProducer",
    ComponentName = cms.string('Chi2'),
    MaxChi2 = cms.double(30),
    MaxDisplacement = cms.double(0.5),
    MaxSagitta = cms.double(2),
    MinPtForHitRecoveryInGluedDet = cms.double(1000000000000),
    MinimalTolerance = cms.double(0.5),
    appendToDataLabel = cms.string(''),
    nSigma = cms.double(3)
)


process.DTGeometryESModule = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.EcalBarrelGeometryFromDBEP = cms.ESProducer("EcalBarrelGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalElectronicsMappingBuilder = cms.ESProducer("EcalElectronicsMappingBuilder")


process.EcalEndcapGeometryFromDBEP = cms.ESProducer("EcalEndcapGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalLaserCorrectionService = cms.ESProducer("EcalLaserCorrectionService")


process.EcalPreshowerGeometryFromDBEP = cms.ESProducer("EcalPreshowerGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.EcalTrigTowerConstituentsMapBuilder = cms.ESProducer("EcalTrigTowerConstituentsMapBuilder",
    MapFile = cms.untracked.string('Geometry/EcalMapping/data/EndCap_TTMap.txt')
)


process.GlobalTrackingGeometryESProducer = cms.ESProducer("GlobalTrackingGeometryESProducer")


process.HcalAlignmentEP = cms.ESProducer("HcalAlignmentEP")


process.HcalGeometryFromDBEP = cms.ESProducer("HcalGeometryFromDBEP",
    applyAlignment = cms.bool(True)
)


process.KFUpdatorESProducer = cms.ESProducer("KFUpdatorESProducer",
    ComponentName = cms.string('KFUpdator')
)


process.MaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterial'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('alongMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.MeasurementTracker = cms.ESProducer("MeasurementTrackerESProducer",
    ComponentName = cms.string(''),
    DebugPixelModuleQualityDB = cms.untracked.bool(False),
    DebugPixelROCQualityDB = cms.untracked.bool(False),
    DebugStripAPVFiberQualityDB = cms.untracked.bool(False),
    DebugStripModuleQualityDB = cms.untracked.bool(False),
    DebugStripStripQualityDB = cms.untracked.bool(False),
    HitMatcher = cms.string('StandardMatcher'),
    MaskBadAPVFibers = cms.bool(True),
    PixelCPE = cms.string('PixelCPEGeneric'),
    SiStripQualityLabel = cms.string(''),
    StripCPE = cms.string('StripCPEfromTrackAngle'),
    UsePixelModuleQualityDB = cms.bool(True),
    UsePixelROCQualityDB = cms.bool(True),
    UseStripAPVFiberQualityDB = cms.bool(True),
    UseStripModuleQualityDB = cms.bool(True),
    UseStripStripQualityDB = cms.bool(True),
    badStripCuts = cms.PSet(
        TEC = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TIB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TID = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        ),
        TOB = cms.PSet(
            maxBad = cms.uint32(4),
            maxConsecutiveBad = cms.uint32(2)
        )
    )
)


process.MuonDetLayerGeometryESProducer = cms.ESProducer("MuonDetLayerGeometryESProducer")


process.OppositeMaterialPropagator = cms.ESProducer("PropagatorWithMaterialESProducer",
    ComponentName = cms.string('PropagatorWithMaterialOpposite'),
    Mass = cms.double(0.105),
    MaxDPhi = cms.double(1.6),
    PropagationDirection = cms.string('oppositeToMomentum'),
    SimpleMagneticField = cms.string(''),
    ptMin = cms.double(-1.0),
    useRungeKutta = cms.bool(False)
)


process.ParabolicParametrizedMagneticFieldProducer = cms.ESProducer("AutoParametrizedMagneticFieldProducer",
    label = cms.untracked.string('ParabolicMf'),
    valueOverride = cms.int32(18268),
    version = cms.string('Parabolic')
)


process.RPCGeometryESModule = cms.ESProducer("RPCGeometryESModule",
    compatibiltyWith11 = cms.untracked.bool(True),
    useDDD = cms.untracked.bool(False)
)


process.SiStripRecHitMatcherESProducer = cms.ESProducer("SiStripRecHitMatcherESProducer",
    ComponentName = cms.string('StandardMatcher'),
    NSigmaInside = cms.double(3.0),
    PreFilter = cms.bool(False)
)


process.SteppingHelixPropagatorAlong = cms.ESProducer("SteppingHelixPropagatorESProducer",
    ApplyRadX0Correction = cms.bool(True),
    AssumeNoMaterial = cms.bool(False),
    ComponentName = cms.string('SteppingHelixPropagatorAlong'),
    NoErrorPropagation = cms.bool(False),
    PropagationDirection = cms.string('alongMomentum'),
    SetVBFPointer = cms.bool(False),
    VBFName = cms.string('VolumeBasedMagneticField'),
    debug = cms.bool(False),
    endcapShiftInZNeg = cms.double(0.0),
    endcapShiftInZPos = cms.double(0.0),
    returnTangentPlane = cms.bool(True),
    sendLogWarning = cms.bool(False),
    useEndcapShiftsInZ = cms.bool(False),
    useInTeslaFromMagField = cms.bool(False),
    useIsYokeFlag = cms.bool(True),
    useMagVolumes = cms.bool(True),
    useMatVolumes = cms.bool(True),
    useTuningForL2Speed = cms.bool(False)
)


process.StripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('SimpleStripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.StripCPEfromTrackAngleESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('StripCPEfromTrackAngle'),
    ComponentType = cms.string('StripCPEfromTrackAngle'),
    parameters = cms.PSet(
        mLC_P0 = cms.double(-0.326),
        mLC_P1 = cms.double(0.618),
        mLC_P2 = cms.double(0.3),
        mTEC_P0 = cms.double(-1.885),
        mTEC_P1 = cms.double(0.471),
        mTIB_P0 = cms.double(-0.742),
        mTIB_P1 = cms.double(0.202),
        mTID_P0 = cms.double(-1.427),
        mTID_P1 = cms.double(0.433),
        mTOB_P0 = cms.double(-1.026),
        mTOB_P1 = cms.double(0.253),
        maxChgOneMIP = cms.double(6000.0),
        useLegacyError = cms.bool(False)
    )
)


process.TrackerRecoGeometryESProducer = cms.ESProducer("TrackerRecoGeometryESProducer")


process.VolumeBasedMagneticFieldESProducer = cms.ESProducer("VolumeBasedMagneticFieldESProducerFromDB",
    debugBuilder = cms.untracked.bool(False),
    label = cms.untracked.string(''),
    valueOverride = cms.int32(18268)
)


process.ZdcGeometryFromDBEP = cms.ESProducer("ZdcGeometryFromDBEP",
    applyAlignment = cms.bool(False)
)


process.beamHaloNavigationSchoolESProducer = cms.ESProducer("NavigationSchoolESProducer",
    ComponentName = cms.string('BeamHaloNavigationSchool'),
    SimpleMagneticField = cms.string('')
)


process.cosmicsNavigationSchoolESProducer = cms.ESProducer("SkippingLayerCosmicNavigationSchoolESProducer",
    ComponentName = cms.string('CosmicNavigationSchool'),
    allSelf = cms.bool(True),
    noPXB = cms.bool(False),
    noPXF = cms.bool(False),
    noTEC = cms.bool(False),
    noTIB = cms.bool(False),
    noTID = cms.bool(False),
    noTOB = cms.bool(False),
    selfSearch = cms.bool(True)
)


process.fakeForIdealAlignment = cms.ESProducer("FakeAlignmentProducer",
    appendToDataLabel = cms.string('fakeForIdeal')
)


process.hcalDDDRecConstants = cms.ESProducer("HcalDDDRecConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalDDDSimConstants = cms.ESProducer("HcalDDDSimConstantsESModule",
    appendToDataLabel = cms.string('')
)


process.hcalTopologyIdeal = cms.ESProducer("HcalTopologyIdealEP",
    Exclude = cms.untracked.string(''),
    MergePosition = cms.untracked.bool(False),
    appendToDataLabel = cms.string('')
)


process.hcal_db_producer = cms.ESProducer("HcalDbProducer",
    dump = cms.untracked.vstring(''),
    file = cms.untracked.string('')
)


process.idealForDigiCSCGeometry = cms.ESProducer("CSCGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    debugV = cms.untracked.bool(False),
    useCentreTIOffsets = cms.bool(False),
    useDDD = cms.bool(False),
    useGangedStripsInME1a = cms.bool(True),
    useOnlyWiresInME1a = cms.bool(False),
    useRealWireGeometry = cms.bool(True)
)


process.idealForDigiDTGeometry = cms.ESProducer("DTGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.idealForDigiTrackerGeometry = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string('fakeForIdeal'),
    appendToDataLabel = cms.string('idealForDigi'),
    applyAlignment = cms.bool(False),
    fromDDD = cms.bool(False)
)


process.navigationSchoolESProducer = cms.ESProducer("NavigationSchoolESProducer",
    ComponentName = cms.string('SimpleNavigationSchool'),
    SimpleMagneticField = cms.string('')
)


process.siPixelQualityESProducer = cms.ESProducer("SiPixelQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelQualityFromDbRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiPixelDetVOffRcd'),
            tag = cms.string('')
        )
    )
)


process.siStripBackPlaneCorrectionDepESProducer = cms.ESProducer("SiStripBackPlaneCorrectionDepESProducer",
    BackPlaneCorrectionDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    BackPlaneCorrectionPeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripBackPlaneCorrectionRcd')
    ),
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    )
)


process.siStripGainESProducer = cms.ESProducer("SiStripGainESProducer",
    APVGain = cms.VPSet(
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGainRcd')
        ), 
        cms.PSet(
            Label = cms.untracked.string(''),
            NormalizationFactor = cms.untracked.double(1.0),
            Record = cms.string('SiStripApvGain2Rcd')
        )
    ),
    AutomaticNormalization = cms.bool(False),
    appendToDataLabel = cms.string(''),
    printDebug = cms.untracked.bool(False)
)


process.siStripLorentzAngleDepESProducer = cms.ESProducer("SiStripLorentzAngleDepESProducer",
    LatencyRecord = cms.PSet(
        label = cms.untracked.string(''),
        record = cms.string('SiStripLatencyRcd')
    ),
    LorentzAngleDeconvMode = cms.PSet(
        label = cms.untracked.string('deconvolution'),
        record = cms.string('SiStripLorentzAngleRcd')
    ),
    LorentzAnglePeakMode = cms.PSet(
        label = cms.untracked.string('peak'),
        record = cms.string('SiStripLorentzAngleRcd')
    )
)


process.siStripQualityESProducer = cms.ESProducer("SiStripQualityESProducer",
    ListOfRecordToMerge = cms.VPSet(
        cms.PSet(
            record = cms.string('SiStripDetVOffRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripDetCablingRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('RunInfoRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadChannelRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadFiberRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadModuleRcd'),
            tag = cms.string('')
        ), 
        cms.PSet(
            record = cms.string('SiStripBadStripRcd'),
            tag = cms.string('')
        )
    ),
    PrintDebugOutput = cms.bool(False),
    ReduceGranularity = cms.bool(False),
    ThresholdForReducedGranularity = cms.double(0.3),
    UseEmptyRunInfo = cms.bool(False),
    appendToDataLabel = cms.string('')
)


process.sistripconn = cms.ESProducer("SiStripConnectivity")


process.stripCPEESProducer = cms.ESProducer("StripCPEESProducer",
    ComponentName = cms.string('stripCPE'),
    ComponentType = cms.string('SimpleStripCPE'),
    parameters = cms.PSet(

    )
)


process.trackerGeometryDB = cms.ESProducer("TrackerDigiGeometryESModule",
    alignmentsLabel = cms.string(''),
    appendToDataLabel = cms.string(''),
    applyAlignment = cms.bool(True),
    fromDDD = cms.bool(False)
)


process.trackerNumberingGeometryDB = cms.ESProducer("TrackerGeometricDetESModule",
    appendToDataLabel = cms.string(''),
    fromDDD = cms.bool(False)
)


process.trackerTopology = cms.ESProducer("TrackerTopologyEP",
    appendToDataLabel = cms.string('')
)


process.trajectoryCleanerBySharedHits = cms.ESProducer("TrajectoryCleanerESProducer",
    ComponentName = cms.string('TrajectoryCleanerBySharedHits'),
    ComponentType = cms.string('TrajectoryCleanerBySharedHits'),
    MissingHitPenalty = cms.double(20.0),
    ValidHitBonus = cms.double(5.0),
    allowSharedFirstHit = cms.bool(True),
    fractionShared = cms.double(0.19)
)


process.ttrhbwr = cms.ESProducer("TkTransientTrackingRecHitBuilderESProducer",
    ComponentName = cms.string('WithTrackAngle'),
    ComputeCoarseLocalPositionFromDisk = cms.bool(False),
    Matcher = cms.string('StandardMatcher'),
    PixelCPE = cms.string('PixelCPEGeneric'),
    StripCPE = cms.string('StripCPEfromTrackAngle')
)


process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string(''),
        authenticationSystem = cms.untracked.int32(0),
        messageLevel = cms.untracked.int32(0),
        security = cms.untracked.string('')
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('94X_mc2017_realistic_v17'),
    pfnPostfix = cms.untracked.string(''),
    pfnPrefix = cms.untracked.string(''),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)


process.HcalTimeSlewEP = cms.ESSource("HcalTimeSlewEP",
    appendToDataLabel = cms.string('HBHE'),
    timeSlewParametersM2 = cms.VPSet(
        cms.PSet(
            slope = cms.double(-3.178648),
            tmax = cms.double(16.0),
            tzero = cms.double(23.960177)
        ), 
        cms.PSet(
            slope = cms.double(-1.556668),
            tmax = cms.double(10.0),
            tzero = cms.double(13.307784)
        ), 
        cms.PSet(
            slope = cms.double(-1.075824),
            tmax = cms.double(6.25),
            tzero = cms.double(9.109694)
        )
    ),
    timeSlewParametersM3 = cms.VPSet(
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(15.5),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-3.2),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(32.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        ), 
        cms.PSet(
            cap = cms.double(6.0),
            tspar0 = cms.double(12.2999),
            tspar0_siPM = cms.double(0.0),
            tspar1 = cms.double(-2.19142),
            tspar1_siPM = cms.double(0.0),
            tspar2 = cms.double(0.0),
            tspar2_siPM = cms.double(0.0)
        )
    )
)


process.eegeom = cms.ESSource("EmptyESSource",
    firstValid = cms.vuint32(1),
    iovIsRunNotTime = cms.bool(True),
    recordName = cms.string('EcalMappingRcd')
)


process.es_hardcode = cms.ESSource("HcalHardcodeCalibrations",
    GainWidthsForTrigPrims = cms.bool(False),
    HBRecalibration = cms.bool(False),
    HBmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHB.txt'),
    HBreCalibCutoff = cms.double(20.0),
    HERecalibration = cms.bool(False),
    HEmeanenergies = cms.FileInPath('CalibCalorimetry/HcalPlugins/data/meanenergiesHE.txt'),
    HEreCalibCutoff = cms.double(20.0),
    HFRecalParameterBlock = cms.PSet(
        HFdepthOneParameterA = cms.vdouble(
            0.004123, 0.00602, 0.008201, 0.010489, 0.013379, 
            0.016997, 0.021464, 0.027371, 0.034195, 0.044807, 
            0.058939, 0.125497
        ),
        HFdepthOneParameterB = cms.vdouble(
            -4e-06, -2e-06, 0.0, 4e-06, 1.5e-05, 
            2.6e-05, 6.3e-05, 8.4e-05, 0.00016, 0.000107, 
            0.000425, 0.000209
        ),
        HFdepthTwoParameterA = cms.vdouble(
            0.002861, 0.004168, 0.0064, 0.008388, 0.011601, 
            0.014425, 0.018633, 0.023232, 0.028274, 0.035447, 
            0.051579, 0.086593
        ),
        HFdepthTwoParameterB = cms.vdouble(
            -2e-06, -0.0, -7e-06, -6e-06, -2e-06, 
            1e-06, 1.9e-05, 3.1e-05, 6.7e-05, 1.2e-05, 
            0.000157, -3e-06
        )
    ),
    HFRecalibration = cms.bool(False),
    SiPMCharacteristics = cms.VPSet(
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(36000)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(2500)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.17),
            nonlin1 = cms.double(1.00985),
            nonlin2 = cms.double(7.84089e-06),
            nonlin3 = cms.double(2.86282e-10),
            pixels = cms.int32(27370)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.196),
            nonlin1 = cms.double(1.00546),
            nonlin2 = cms.double(6.40239e-06),
            nonlin3 = cms.double(1.27011e-10),
            pixels = cms.int32(38018)
        ), 
        cms.PSet(
            crosstalk = cms.double(0.0),
            nonlin1 = cms.double(1.0),
            nonlin2 = cms.double(0.0),
            nonlin3 = cms.double(0.0),
            pixels = cms.int32(0)
        )
    ),
    hb = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.19),
        gainWidth = cms.vdouble(0.0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.285),
        pedestalWidth = cms.double(0.809),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.49, 1.8, 7.2, 37.9),
        qieSlope = cms.vdouble(0.912, 0.917, 0.922, 0.923),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(8)
    ),
    hbUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(150),
            intlumiToNeutrons = cms.double(367000000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(-5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    he = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.23),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(125),
        pedestal = cms.double(3.163),
        pedestalWidth = cms.double(0.9698),
        photoelectronsToAnalog = cms.double(0.3305),
        qieOffset = cms.vdouble(-0.38, 2.0, 7.6, 39.6),
        qieSlope = cms.vdouble(0.912, 0.916, 0.92, 0.922),
        qieType = cms.int32(0),
        recoShape = cms.int32(105),
        zsThreshold = cms.int32(9)
    ),
    heUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.01, 0.015),
        doRadiationDamage = cms.bool(True),
        gain = cms.vdouble(0.0006252),
        gainWidth = cms.vdouble(0),
        mcShape = cms.int32(206),
        pedestal = cms.double(17.3),
        pedestalWidth = cms.double(1.5),
        photoelectronsToAnalog = cms.double(40.0),
        qieOffset = cms.vdouble(0.0, 0.0, 0.0, 0.0),
        qieSlope = cms.vdouble(0.05376, 0.05376, 0.05376, 0.05376),
        qieType = cms.int32(2),
        radiationDamage = cms.PSet(
            depVsNeutrons = cms.vdouble(5.543e-10, 8.012e-10),
            depVsTemp = cms.double(0.0631),
            intlumiOffset = cms.double(75),
            intlumiToNeutrons = cms.double(29200000.0),
            temperatureBase = cms.double(20),
            temperatureNew = cms.double(5)
        ),
        recoShape = cms.int32(206),
        zsThreshold = cms.int32(16)
    ),
    hf = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(9.354),
        pedestalWidth = cms.double(2.516),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(-0.87, 1.4, 7.8, -29.6),
        qieSlope = cms.vdouble(0.359, 0.358, 0.36, 0.367),
        qieType = cms.int32(0),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    hfUpgrade = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.14, 0.135),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(301),
        pedestal = cms.double(13.33),
        pedestalWidth = cms.double(3.33),
        photoelectronsToAnalog = cms.double(0.0),
        qieOffset = cms.vdouble(0.0697, -0.7405, 12.38, -671.9),
        qieSlope = cms.vdouble(0.297, 0.298, 0.298, 0.313),
        qieType = cms.int32(1),
        recoShape = cms.int32(301),
        zsThreshold = cms.int32(-9999)
    ),
    ho = cms.PSet(
        darkCurrent = cms.vdouble(0.0),
        doRadiationDamage = cms.bool(False),
        gain = cms.vdouble(0.006, 0.0087),
        gainWidth = cms.vdouble(0.0, 0.0),
        mcShape = cms.int32(201),
        pedestal = cms.double(12.06),
        pedestalWidth = cms.double(0.6285),
        photoelectronsToAnalog = cms.double(4.0),
        qieOffset = cms.vdouble(-0.44, 1.4, 7.1, 38.5),
        qieSlope = cms.vdouble(0.907, 0.915, 0.92, 0.921),
        qieType = cms.int32(0),
        recoShape = cms.int32(201),
        zsThreshold = cms.int32(24)
    ),
    iLumi = cms.double(-1.0),
    killHE = cms.bool(False),
    testHEPlan1 = cms.bool(False),
    testHFQIE10 = cms.bool(False),
    toGet = cms.untracked.vstring('GainWidths'),
    useHBUpgrade = cms.bool(False),
    useHEUpgrade = cms.bool(False),
    useHFUpgrade = cms.bool(False),
    useHOUpgrade = cms.bool(True),
    useIeta18depth1 = cms.bool(True),
    useLayer0Weight = cms.bool(False)
)


process.prefer("es_hardcode")

process.egammaUpdatorTask = cms.Task()


process.egammaPostRecoPatUpdatorTask = cms.Task(process.slimmedElectrons, process.slimmedPhotons)


process.egammaScaleSmearTask = cms.Task(process.calibratedPatElectrons, process.calibratedPatPhotons)


process.egammaVIDTask = cms.Task()


process.egammaUpdatorSeq = cms.Sequence(process.egammaUpdatorTask)


process.eleSequence = cms.Sequence(process.goodElectrons+process.vetoElectrons)


process.egammaVIDSeq = cms.Sequence(process.egammaVIDTask)


process.muSequence = cms.Sequence(process.goodMuons+process.looseMuons)


process.egammaScaleSmearSeq = cms.Sequence(process.egammaScaleSmearTask)


process.leptonicVSequence = cms.Sequence(process.Wtomunu+process.Wtoenu+process.leptonicV)


process.NJetsSequence = cms.Sequence(process.goodAK4Jets+process.cleanAK4Jets)


process.egammaPostRecoPatUpdatorSeq = cms.Sequence(process.egammaPostRecoPatUpdatorTask)


process.photonSequence = cms.Sequence(process.goodPhotons)


process.metfilterSequence = cms.Sequence(process.BadPFMuonFilter+process.BadChargedCandidateFilter)


process.jetSequence = cms.Sequence(process.NJetsSequence)


process.egammaPostRecoSeq = cms.Sequence(process.egammaUpdatorSeq+process.egammaScaleSmearSeq+process.egammaVIDSeq+process.egammaPostRecoPatUpdatorSeq)


process.leptonSequence = cms.Sequence(process.muSequence+process.eleSequence+process.leptonicVSequence+process.leptonicVSelector+process.leptonicVFilter)


process.analysis = cms.Path(process.JetUserData+process.leptonSequence+process.jetSequence+process.metfilterSequence+process.ecalBadCalibReducedMINIAODFilter+process.prefiringweight+process.treeDumper)


