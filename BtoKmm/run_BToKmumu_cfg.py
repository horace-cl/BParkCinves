from FWCore.ParameterSet.VarParsing import VarParsing
import FWCore.ParameterSet.Config as cms

####################################################################################################
#################################### PARSEO DE VARIABLES (ARGS) ####################################
####################################################################################################

options = VarParsing('python')

options.register('isMC', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('globalTag', 'NOTSET',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "Set global tag"
)
options.register('wantSummary', True,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('wantFullRECO', False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Run this on real data"
)
options.register('reportEvery', 100,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "report every N events"
)
options.register('skip', 0,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "skip first N events"
)

options.register('maxE', 3000,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Maximum number of events"
)

options.register('tg', '3k_parseArgs',
    VarParsing.multiplicity.singleton,
    VarParsing.varType.string,
    "tag for outputfile"
)

options.parseArguments()
options.setDefault('maxEvents', options.maxE)
options.setDefault('tag', options.tg)
 
extension = {False : 'data', True : 'mc'}
outputFileNANO = cms.untracked.string('_'.join(['BParkNANO', extension[options.isMC], options.tag])+'.root')
outputFileFEVT = cms.untracked.string('_'.join(['BParkFullEvt', extension[options.isMC], options.tag])+'.root')
if not options.inputFiles:
	if options.isMC:
		options.inputFiles = ['/store/cmst3/group/bpark/BToKmumu_1000Events_MINIAOD.root']
	else:
		options.inputFiles = [
		# Bloque  --->  /ParkingBPH4/Run2018B-05May2019-v2/MINIAOD#62e4110e-94d0-4a3b-b8e4-4df4863cd558   <---   con 637K eventos
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root', 
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/F7E7EF39-476F-1C48-95F7-74CB5C7A542C.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/C8242DD9-5C7C-304F-A652-366512704CAC.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/C366AD8C-D054-1641-B51C-A4EA23B849AC.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/BDB59105-A1F2-E947-A9A7-377DF4BF246D.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/B7A9D939-4DCF-1E4B-B51F-D0B8341E98C5.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/AFAE473F-CCF6-4744-9D31-43DC25B17699.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/902D0E02-4CBD-054C-885A-84806BB4A050.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/82B8B790-F38F-ED45-BDC3-C939E89E3382.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/81B2BC95-8E17-124D-A947-202576A42CB0.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root', #Default
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/50D2AEAA-A149-B44A-928D-59BC4A27DD02.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/4DCF5A74-7ECC-8148-A672-07DC7EDAC8A5.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/46B692A6-2186-A146-A1C0-D7C2EEFC590C.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/382A8ED8-B10F-C042-A965-93B0D223BA98.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/2D37F175-285B-4C4C-9F8A-187A41ED226A.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/2A239304-3909-1F4E-937A-4429CABB4072.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/12FCD72D-CD84-4542-AF68-DDAB511A224C.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/12A7CB03-BAA7-084E-979C-B98929C001E7.root',
		'/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/08C016D2-8DEA-DF44-BA83-05A3ABB22139.root',
				 			]
    # options.inputFiles = ['/store/data/Run2018B/ParkingBPH4/MINIAOD/05May2019-v2/230000/6B5A24B1-0E6E-504B-8331-BD899EB60110.root'] if not options.isMC else \
    #                      ['/store/cmst3/group/bpark/BToKmumu_1000Events_MINIAOD.root']
annotation = '%s nevts:%d' % (outputFileNANO, options.maxEvents)





####################################################################################################
################################# ERAS E INSTANCIACION DEL PROCESO #################################
####################################################################################################

from Configuration.StandardSequences.Eras import eras
process = cms.Process('BParkNANO',eras.Run2_2018)




#####################################################################################################
#################################     QUE PASA EXACTAMENTE AQUI     #################################
#####################################################################################################

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load('PhysicsTools.BParkingNano.nanoBPark_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.load('PhysicsTools.BParkingNano.BToKMuMu_cff')
process.load('config_BtoKmumu_cff')
 




#####################################################################################################
################################### INPUT AND OUTPUT CONFIGURATION ##################################
#####################################################################################################


################################ Input source
process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles),
    secondaryFileNames = cms.untracked.vstring(),
    skipEvents=cms.untracked.uint32(options.skip),
)

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary),
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
)

process.nanoMetadata.strings.tag = annotation




################################ Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string(annotation),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)




################################ Output definition

process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)


process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileFEVT,
    outputCommands = (cms.untracked.vstring('keep *',
                                            'drop *_*_SelectedTransient*_*',
                     )),
    splitLevel = cms.untracked.int32(0)
)

process.NANOAODoutput = cms.OutputModule("NanoAODOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(9),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('NANOAOD'),
        filterName = cms.untracked.string('')
    ),
    fileName = outputFileNANO,
    outputCommands = cms.untracked.vstring(
      # A primer vista no hay diferencia entre usar o no drop *
      'drop *',
      # Es en esta linea donde se indica que se almacenen los branches en Event a parte de
      #  b'run'
      #	 b'luminosityBlock'
      #	 b'event'
      "keep nanoaodFlatTable_*_*_*",     # event data
      #Esta linea se encarga de agregar el tag (Tree)
      "keep nanoaodUniqueString_nanoMetadata_*_*",   # basic metadata
    )

)
 



#####################################################################################################
######################################### GLOBAL TAG ################################################
#####################################################################################################

globaltag = '102X_dataRun2_v11' if not options.isMC else '102X_upgrade2018_realistic_v15'
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, globaltag, '')




#####################################################################################################
################################### PATH, ENDPAT Y SCHEDULE #########################################
#####################################################################################################

from config_BtoKmumu_cff import *

process.nanoAOD_KMuMu_step = cms.Path(nanoSequence  + nanoBKMuMuSequence + CountBToKmumu )

process.endjob_step = cms.EndPath(process.endOfProcess)

process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

process.NANOAODoutput_step = cms.EndPath(process.NANOAODoutput)

process.schedule = cms.Schedule(
                                process.nanoAOD_KMuMu_step,
                                process.endjob_step, 
                                process.NANOAODoutput_step
                               )


 
 

#####################################################################################################
###################################### WHAT IS THIS FOR? ############################################
#####################################################################################################



from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

process.NANOAODoutput.SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring(
                                   'nanoAOD_KMuMu_step', 
                                   )
)


### from https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3287/1/1/1/1/1.html
process.add_(cms.Service('InitRootHandlers', EnableIMT = cms.untracked.bool(False)))
process.NANOAODoutput.fakeNameForCrab=cms.untracked.bool(True)    

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
