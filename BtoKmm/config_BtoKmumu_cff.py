# Lets try to reconstruct the channel B -> K* Mu Mu

import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *



from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.nano_cff import *
from PhysicsTools.NanoAOD.vertices_cff import *
from PhysicsTools.NanoAOD.NanoAODEDMEventContent_cff import *
from PhysicsTools.BParkingNano.trgbits_cff import *

##for gen and trigger muon
from PhysicsTools.BParkingNano.genparticlesBPark_cff import *
from PhysicsTools.BParkingNano.particlelevelBPark_cff import *
from PhysicsTools.BParkingNano.triggerObjectsBPark_cff import *
from PhysicsTools.BParkingNano.muonsBPark_cff import muonTriggerMatchedTable, muonTriggerMatchedTables

## filtered input collections
from PhysicsTools.BParkingNano.electronsBPark_cff import * 
#from PhysicsTools.BParkingNano.tracksBPark_cff import *
from PhysicsTools.BParkingNano.common_cff import *




######################################################################################################
############################################ Muons ###################################################
######################################################################################################

muonTrgSelector = cms.EDProducer("MuonTriggerSelector",
                                 muonCollection = cms.InputTag("slimmedMuons"), #same collection as in NanoAOD                                                           
                                 bits = cms.InputTag("TriggerResults","","HLT"),
                                 prescales = cms.InputTag("patTrigger"),
                                 objects = cms.InputTag("slimmedPatTrigger"),
                                 vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                 
                                 ##for the output trigger matched collection
                                 maxdR_matching = cms.double(0.1),
                                 
                                 ## for the output selected collection (tag + all compatible in dZ)
                                 dzForCleaning_wrtTrgMuon = cms.double(1.),

                                 ptMin = cms.double(0.5),
                                 absEtaMax = cms.double(2.4),
                                 # keeps only muons with at soft Quality flag
                                 softMuonsOnly = cms.bool(False)
                             )

countTrgMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "trgMuons")
)

countSelectedMuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(2),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonTrgSelector", "SelectedMuons")
)


muonBParkTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("muonTrgSelector:SelectedMuons"),
    cut = cms.string(""), #we should not filter on cross linked collections
    name = cms.string("Muon"),
    doc  = cms.string("slimmedMuons for BPark after basic selection"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for the muons
    variables = cms.PSet(CandVars,
        # ptErr   = Var("bestTrack().ptError()", float, doc = "ptError of the muon track", precision=6),
        # dz = Var("dB('PVDZ')",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        # dzErr = Var("abs(edB('PVDZ'))",float,doc="dz uncertainty, in cm",precision=6),
        # dxy = Var("dB('PV2D')",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        # dxyErr = Var("edB('PV2D')",float,doc="dxy uncertainty, in cm",precision=6),
        # vx = Var("vx()",float,doc="x coordinate of vertex position, in cm",precision=6),
        # vy = Var("vy()",float,doc="y coordinate of vertex position, in cm",precision=6),
        # vz = Var("vz()",float,doc="z coordinate of vertex position, in cm",precision=6),
        # px = Var("px()",float,doc="x momentum",precision=6),
        # py = Var("py()",float,doc="y momentum",precision=6),
        # pz = Var("pz()",float,doc="z momentum",precision=6),
        # energy = Var("energy()", float, doc = "energy of the muon", precision=6),
        # ip3d = Var("abs(dB('PV3D'))",float,doc="3D impact parameter wrt first PV, in cm",precision=10),
        # sip3d = Var("abs(dB('PV3D')/edB('PV3D'))",float,doc="3D impact parameter significance wrt first PV",precision=10),
        # segmentComp   = Var("segm entCompatibility()", float, doc = "muon segment compatibility", precision=14), # keep higher precision since people have cuts with 3 digits on this
        # nStations = Var("numberOfMatchedStations", int, doc = "number of matched stations with default arbitration (segment & track)"),
        # nTrackerLayers = Var("innerTrack().hitPattern().trackerLayersWithMeasurement()", int, doc = "number of layers in the tracker"),
        # pfRelIso03_chg = Var("pfIsolationR03().sumChargedHadronPt/pt",float,doc="PF relative isolation dR=0.3, charged component"),
        # pfRelIso03_all = Var("(pfIsolationR03().sumChargedHadronPt + max(pfIsolationR03().sumNeutralHadronEt + pfIsolationR03().sumPhotonEt - pfIsolationR03().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.3, total (deltaBeta corrections)"),
        # pfRelIso04_all = Var("(pfIsolationR04().sumChargedHadronPt + max(pfIsolationR04().sumNeutralHadronEt + pfIsolationR04().sumPhotonEt - pfIsolationR04().sumPUPt/2,0.0))/pt",float,doc="PF relative isolation dR=0.4, total (deltaBeta corrections)"),
        # tightCharge = Var("?(muonBestTrack().ptError()/muonBestTrack().pt() < 0.2)?2:0",int,doc="Tight charge criterion using pterr/pt of muonBestTrack (0:fail, 2:pass)"),
        isPFcand = Var("isPFMuon",bool,doc="muon is PF candidate"),
        isGlobal = Var("isGlobalMuon",bool,doc="muon is global muon"),
        isTracker = Var("isTrackerMuon",bool,doc="muon is tracker muon"),
        # mediumId = Var("passed('CutBasedIdMedium')",bool,doc="cut-based ID, medium WP"),
        # # mediumPromptId = Var("passed('CutBasedIdMediumPrompt')",bool,doc="cut-based ID, medium prompt WP"),
        # tightId = Var("passed('CutBasedIdTight')",bool,doc="cut-based ID, tight WP"),
        # softId = Var("passed('SoftCutBasedId')",bool,doc="soft cut-based ID"),
        # # softMvaId = Var("passed('SoftMvaId')",bool,doc="soft MVA ID"),
        # # highPtId = Var("?passed('CutBasedIdGlobalHighPt')?2:passed('CutBasedIdTrkHighPt')","uint8",doc="high-pT cut-based ID (1 = tracker high pT, 2 = global high pT, which includes tracker high pT)"),
        # pfIsoId = Var("passed('PFIsoVeryLoose')+passed('PFIsoLoose')+passed('PFIsoMedium')+passed('PFIsoTight')+passed('PFIsoVeryTight')+passed('PFIsoVeryVeryTight')","uint8",doc="PFIso ID from miniAOD selector (1=PFIsoVeryLoose, 2=PFIsoLoose, 3=PFIsoMedium, 4=PFIsoTight, 5=PFIsoVeryTight, 6=PFIsoVeryVeryTight)"),
        # tkIsoId = Var("?passed('TkIsoTight')?2:passed('TkIsoLoose')","uint8",doc="TkIso ID (1=TkIsoLoose, 2=TkIsoTight)"),
        # # mvaId = Var("passed('MvaLoose')+passed('MvaMedium')+passed('MvaTight')","uint8",doc="Mva ID from miniAOD selector (1=MvaLoose, 2=MvaMedium, 3=MvaTight)"),
        # # miniIsoId = Var("passed('MiniIsoLoose')+passed('MiniIsoMedium')+passed('MiniIsoTight')+passed('MiniIsoVeryTight')","uint8",doc="MiniIso ID from miniAOD selector (1=MiniIsoLoose, 2=MiniIsoMedium, 3=MiniIsoTight, 4=MiniIsoVeryTight)"),
        # # multiIsoId = Var("?passed('MultiIsoMedium')?2:passed('MultiIsoLoose')","uint8",doc="MultiIsoId from miniAOD selector (1=MultiIsoLoose, 2=MultiIsoMedium)"),
        # triggerIdLoose = Var("passed('TriggerIdLoose')",bool,doc="TriggerIdLoose ID"),
        # # inTimeMuon = Var("passed('InTimeMuon')",bool,doc="inTimeMuon ID"),
        isTriggering = Var("userInt('isTriggering')", int,doc="flag the reco muon is also triggering")
    ),
)

muonBParkSequence = cms.Sequence(muonTrgSelector * countTrgMuons * countSelectedMuons)
# muonBParkMC = cms.Sequence(muonBParkSequence + muonsBParkMCMatchForTable + selectedMuonsMCMatchEmbedded + muonBParkMCTable)
muonBParkTables = cms.Sequence(muonBParkTable)
# muonTriggerMatchedTables = cms.Sequence(muonTriggerMatchedTable)




#######################################################################################################
############################################ Tracks ###################################################
#######################################################################################################


tracksBPark = cms.EDProducer('TrackMerger',
                             beamSpot   = cms.InputTag("offlineBeamSpot"),
                             trgMuon    = cms.InputTag("muonTrgSelector:trgMuons"),
                             tracks     = cms.InputTag("packedPFCandidates"),
                             lostTracks = cms.InputTag("lostTracks"),
                             trkPtCut = cms.double(1.0),
                             muons      = cms.InputTag("slimmedMuons"),
                             pfElectrons= cms.InputTag("slimmedElectrons"),
                             vertices   = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             trkEtaCut = cms.double(2.5),
                             dzTrg_cleaning = cms.double(1.),
                             drTrg_Cleaning = cms.double(0.03),
                             dcaSig = cms.double(-100000),
                             trkNormChiMin = cms.int32(-1),
                             trkNormChiMax = cms.int32(-1)
                             )

trackBParkTable = cms.EDProducer(
    "SimpleCompositeCandidateFlatTableProducer",
    src = cms.InputTag("tracksBPark:SelectedTracks"),
    cut = cms.string(""),
    name = cms.string("ProbeTracks"),
    doc  = cms.string("track collection probe side for BPark after basic selection"),
    singleton = cms.bool(False),
    extension = cms.bool(False), 
    variables = cms.PSet(
        P4Vars, 
        vx = Var("vx()", float, doc="x coordinate of vertex position, in cm", precision=10),
        vy = Var("vy()", float, doc="y coordinate of vertex position, in cm", precision=10),
        vz = Var("vz()", float, doc="z coordinate of vertex position, in cm", precision=10),
        # isPacked = Var("userInt('isPacked')",int,doc="track from packedCandidate collection", precision=10),
        # isLostTrk = Var("userInt('isLostTrk')",int,doc="track from lostTrack collection", precision=10),
        px = Var("userFloat('px')",float,doc="pz", precision=10),
        py = Var("userFloat('py')",float,doc="py", precision=10),
        pz = Var("userFloat('pz')",float,doc="pz", precision=10),
        energy = Var("userFloat('energy')",float,doc="energy", precision=10),
        #mass = Var("userFloat('mass')",float,doc="mass", precision=10),
        # dz = Var("userFloat('dz')",float,doc="dz (with sign) wrt first PV, in cm", precision=10),
        # dxy = Var("userFloat('dxy')",float,doc="dxy (with sign) wrt first PV, in cm", precision=10),
        # dzS = Var("userFloat('dzS')", float, doc="dz/err (with sign) wrt first PV, in cm", precision=10),
        # dxyS = Var("userFloat('dxyS')", float, doc="dxy/err (with sign) wrt first PV, in cm", precision=10),
        # DCASig=Var("userFloat('DCASig')", float,doc="significance of xy-distance of closest approach wrt beamspot", precision=10),
        # isMatchedToMuon = Var("userInt('isMatchedToMuon')",bool,doc="track was used to build a muon", precision=10),
        # isMatchedToLooseMuon = Var("userInt('isMatchedToLooseMuon')",bool,doc="track was used to build a muon passing LooseID", precision=10),
        # isMatchedToSoftMuon = Var("userInt('isMatchedToSoftMuon')",bool,doc="track was used to build a muon passing softID", precision=10),
        # isMatchedToMediumMuon = Var("userInt('isMatchedToMediumMuon')",bool,doc="track was used to build a muon passing mediumID", precision=10),
        # isMatchedToEle = Var("userInt('isMatchedToEle')",bool,doc="track was used to build a PF ele", precision=10),
        nValidHits = Var("userInt('nValidHits')", int,doc="Number of valid hits on track", precision=10),
        #dEdXStrip=Var("userFloat('dEdXStrip')", float,doc="dE/dX from strips of associated isolated track"),
        #dEdXPixel=Var("userFloat('dEdXPixel')", float,doc="dE/dX from pixels of associated isolated track"),
        ),
)
 

tracksBParkSequence = cms.Sequence(tracksBPark)
tracksBParkTables = cms.Sequence(trackBParkTable)

#######################################################################################################
############################################ Dimuons ##################################################
#######################################################################################################


muonPairsForKmumu = cms.EDProducer(
    'DiMuonBuilder',
    src = cms.InputTag('muonTrgSelector', 'SelectedMuons'),
    transientTracksSrc = cms.InputTag('muonTrgSelector', 'SelectedTransientMuons'),
    lep1Selection = cms.string('pt > 1.5'),
    lep2Selection = cms.string(''),
    preVtxSelection = cms.string('abs(userCand("l1").vz - userCand("l2").vz) <= 1. && mass() < 4. '
                                 '&& mass() > 2.8 && charge() == 0 && userFloat("lep_deltaR") > 0.03'),
    postVtxSelection = cms.string('userFloat("sv_chi2") < 998 && userFloat("sv_prob") > 1.e-5'),
)

diMuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("muonPairsForKmumu"),
    cut = cms.string(""),
    name = cms.string("dimuons"),
    doc = cms.string("Dimuons Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        #CandVars,
        FittedMass = ufloat("fitted_mass"),
        ErrorFittedMass = ufloat("fitted_massErr"),
        l1_idx = uint("l1_idx"),
        #l1_pt = Var('userCand("l1").pt()', float, doc ='pt of l1 muon'),
        l2_idx = uint("l2_idx"),
        lep_deltaR = Var('userFloat("lep_deltaR")', float,doc="delta R from the two muons"),
        sv_chi2 = Var('userFloat("sv_chi2")', float, doc='chi2 for secondary vertex'),
        sv_prob = Var('userFloat("sv_prob")', float, doc='probability for secondary vertex'),
        lep_vert = Var('abs(userCand("l1").vz - userCand("l2").vz)', float, doc='Distance in vertex z'),
        # fit_l1_pt = ufloat('fitted_l1_pt'),
        # fit_l1_eta = ufloat('fitted_l1_eta'),
        # fit_l1_phi = ufloat('fitted_l1_phi'),
        # fit_l2_pt = ufloat('fitted_l2_pt'),
        # fit_l2_eta = ufloat('fitted_l2_eta'),
        # fit_l2_phi = ufloat('fitted_l2_phi'),
        vtx_ex = ufloat('vtx_ex'),
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        #vtx_exy = ufloat('vtx_exy'),
        vtx_eyx = ufloat('vtx_eyx'),
        vtx_ezy = ufloat('vtx_ezy'),
        vtx_ezx = ufloat('vtx_ezx'),
           )
)

countDimuons = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("muonPairsForKmumu")
)
#######################################################################################################
#######################################  B -> K* Mu Mu  ###############################################
#######################################################################################################
BToKmumu = cms.EDProducer(
    'BToKMMBuilder',
    dileptons = cms.InputTag('muonPairsForKmumu'),
    leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    ##HCL
    trgMuon = cms.InputTag("muonTrgSelector:trgMuons"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"), 
    dzCut = cms.double(1.0),
    ##HCL
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& userFloat("fitted_mass") > 4.8 && userFloat("fitted_mass") < 5.8')
    # preVtxSelection = cms.string(''),
    # postVtxSelection = cms.string(''),
)


BToKmumuOriginal = cms.EDProducer(
    'BToKLLBuilder',
    dileptons = cms.InputTag('muonPairsForKmumu'),
    leptonTransientTracks = muonPairsForKmumu.transientTracksSrc,
    kaons = cms.InputTag('tracksBPark', 'SelectedTracks'),
    kaonsTransientTracks = cms.InputTag('tracksBPark', 'SelectedTransientTracks'),
    beamSpot = cms.InputTag("offlineBeamSpot"),
    tracks = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    kaonSelection = cms.string(''),
    isoTracksSelection = cms.string('pt > 0.7 && abs(eta)<2.5'),
    # This in principle can be different between electrons and muons
    preVtxSelection = cms.string(
        'pt > 3. && userFloat("min_dr") > 0.03'
        '&& mass < 7. && mass > 4.'
        ),
    postVtxSelection = cms.string(
        'userInt("sv_OK") == 1 && userFloat("sv_prob") > 0.001 '
        '&& userFloat("fitted_cos_theta_2D") >= 0'
        '&& userFloat("fitted_mass") > 4.5 && userFloat("fitted_mass") < 6.'
    )
)


BToKmumuTableOriginal = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKmumuOriginal"),
    cut = cms.string(""),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        minDR = ufloat('min_dr'),
        maxDR = ufloat('max_dr'),
        # fit and vtx info
        #chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
        mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),
        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),
        fit_l1_pt = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        fit_l2_pt = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        fit_k_pt = ufloat('fitted_k_pt'),
        fit_k_eta = ufloat('fitted_k_eta'),
        fit_k_phi = ufloat('fitted_k_phi'),
        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        k_iso03  = ufloat('k_iso03'),
        k_iso04  = ufloat('k_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        n_k_used = uint('n_k_used'),
        n_l1_used = uint('n_l1_used'),
        n_l2_used = uint('n_l2_used'),
    )
)

BToKmumuTable = cms.EDProducer(
    'SimpleCompositeCandidateFlatTableProducer',
    src = cms.InputTag("BToKmumu:Bcollection"),
    cut = cms.string(""),
    name = cms.string("BToKMuMu"),
    doc = cms.string("BToKMuMu Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
        l1Idx = uint('l1_idx'),
        l2Idx = uint('l2_idx'),
        kIdx = uint('k_idx'),
        # minDR = ufloat('min_dr'),
        # maxDR = ufloat('max_dr'),
        dr_l1 = ufloat('dr_l1'),
        dr_l2 = ufloat('dr_l2'),
        dz_l1 = ufloat('dz_l1'),
        dz_l2 = ufloat('dz_l2'),
        # #HCL
        # l1Vx = ufloat("l1_vx"),
        # l1Vy = ufloat("l1_vy"),
        # l1Vz = ufloat("l1_vz"),
        # l2Vx = ufloat("l2_vx"),
        # l2Vy = ufloat("l2_vy"),
        # l2Vz = ufloat("l2_vz"),
        # kVx = ufloat("k_vx"),
        # kVy = ufloat("k_vy"),
        # kVz = ufloat("k_vz"),
        # dmuonVx = ufloat("dmuon_vx"),
        # dmuonVy = ufloat("dmuon_vy"),
        # dmuonVz = ufloat("dmuon_vz"),
        # fit and vtx info
        #chi2 = ufloat('sv_chi2'),
        svprob = ufloat('sv_prob'),
        l_xy = ufloat('l_xy'),
        l_xy_unc = ufloat('l_xy_unc'),
        # beamspot_x = ufloat('beamSpot_x'),
        # beamspot_y = ufloat('beamSpot_y'),
        vtx_x = ufloat('vtx_x'),
        vtx_y = ufloat('vtx_y'),
        vtx_z = ufloat('vtx_z'),
        vtx_ex = ufloat('vtx_ex'), ## only saving diagonal elements of the cov matrix
        vtx_ey = ufloat('vtx_ey'),
        vtx_ez = ufloat('vtx_ez'),
        #vtx_exy = ufloat('vtx_exy'),
        vtx_eyx = ufloat('vtx_eyx'),
        vtx_ezy = ufloat('vtx_ezy'),
        vtx_ezx = ufloat('vtx_ezx'),

        # Mll
        mll_raw = Var('userCand("dilepton").mass()', float),
        mll_llfit = Var('userCand("dilepton").userFloat("fitted_mass")', float), # this might not work
        mllErr_llfit = Var('userCand("dilepton").userFloat("fitted_massErr")', float), # this might not work
        # mll_fullfit = ufloat('fitted_mll'),
        # Cos(theta)
        cos2D = ufloat('cos_theta_2D'),
        fit_cos2D = ufloat('fitted_cos_theta_2D'),
        # post-fit momentum
        fit_mass = ufloat('fitted_mass'),
        fit_massErr = ufloat('fitted_massErr'),

        fit_px = ufloat('fitted_px'),
        fit_py = ufloat('fitted_py'),
        fit_pz = ufloat('fitted_pz'),

        fit_pt = ufloat('fitted_pt'),
        fit_eta = ufloat('fitted_eta'),
        fit_phi = ufloat('fitted_phi'),

        fit_l1_pt = ufloat('fitted_l1_pt'),
        fit_l1_eta = ufloat('fitted_l1_eta'),
        fit_l1_phi = ufloat('fitted_l1_phi'),
        l1_charge = ufloat('l1_charge'),

        fit_l2_pt = ufloat('fitted_l2_pt'),
        fit_l2_eta = ufloat('fitted_l2_eta'),
        fit_l2_phi = ufloat('fitted_l2_phi'),
        l2_charge = ufloat('l2_charge'),

        fit_k_pt = ufloat('fitted_k_pt'),
        fit_k_eta = ufloat('fitted_k_eta'),
        fit_k_phi = ufloat('fitted_k_phi'),
        k_charge = ufloat('k_charge'),

        l1_iso03 = ufloat('l1_iso03'),
        l1_iso04 = ufloat('l1_iso04'),
        l2_iso03 = ufloat('l2_iso03'),
        l2_iso04 = ufloat('l2_iso04'),
        k_iso03  = ufloat('k_iso03'),
        k_iso04  = ufloat('k_iso04'),
        b_iso03  = ufloat('b_iso03'),
        b_iso04  = ufloat('b_iso04'),
        # n_k_used = uint('n_k_used'),
        # n_l1_used = uint('n_l1_used'),
        # n_l2_used = uint('n_l2_used'),
        cosThetaMuMu = ufloat("cosTheta_mm"),
        cosThetaKMu = ufloat("cosTheta_km"),

        cosAlpha_0 = ufloat("cosAlpha0"),
        cosAlpha_1 = ufloat("cosAlpha1"),
        cosAlpha_2 = ufloat("cosAlpha2"),

        lxy_pv0 = ufloat("lxy_pv0"),
        lxy_pv1 = ufloat("lxy_pv1"),
        lxy_pv2 = ufloat("lxy_pv2"),

        significance0 = ufloat("significance0"),
        significance1 = ufloat("significance1"),
        significance2 = ufloat("significance2"),
    )
)

CountBToKmumu = cms.EDFilter("PATCandViewCountFilter",
    minNumber = cms.uint32(1),
    maxNumber = cms.uint32(999999),
    src = cms.InputTag("BToKmumu:Bcollection")
)    


BToKMuMuSequence = cms.Sequence(
    (muonPairsForKmumu * countDimuons * diMuTable * BToKmumu)
)





vertexProducer = cms.EDProducer(
    'PrimaryVertexSelector',
    bMesons = cms.InputTag('BToKmumu'),
    trgMuon = cms.InputTag("muonTrgSelector:trgMuons"),
    vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"), 
    dzCut = cms.double(1.0),
)


VertexTable = cms.EDProducer(
    'SimpleCandidateFlatTableProducer',
    src = cms.InputTag("vertexProducer"),
    cut = cms.string(""),
    name = cms.string("PV"),
    doc = cms.string("Primary Vertex Variable"),
    singleton=cms.bool(False),
    extension=cms.bool(False),
    variables=cms.PSet(
        # pre-fit quantities
        CandVars,
    )
)




######################################################################################################
########################################### NanoAOD ##################################################
######################################################################################################


nanoBKMuMuSequence = cms.Sequence( BToKMuMuSequence + BToKmumuTable) #+vertexProducer)#+VertexTable)

nanoSequence = cms.Sequence(# Original
                            #nanoMetadata + 
                            # vertexSequence +           
                            # globalTables + vertexTables + 
                            # triggerObjectBParkTables + l1bits +

                            # # customizeMuonTriggerBPark
                            muonBParkSequence + muonBParkTables +#muonTriggerMatchedTables +
                            
                            # # customizeTrackFilteredBPark
                            tracksBParkSequence #+ tracksBParkTables

                            # #customizeTriggerBitsBPark
                            #trgTables 
                            )
                            # )
