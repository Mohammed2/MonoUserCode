import FWCore.ParameterSet.Config as cms


process = cms.Process('MPGen')


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedRealistic7TeV2011Collision_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')




process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)


# Input source
process.source = cms.Source("EmptySource")

process.options = cms.untracked.PSet(

)


process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('Pythia_monopole_GEN.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    ),
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    )
) 


# Other statements
process.GlobalTag.globaltag = 'MC_44_V7::All'




# Get a unique random number
import os, struct
myseed = abs(struct.unpack('i', os.urandom(4))[0]) % 900000000
process.RandomNumberGeneratorService.generator.initialSeed = myseed 


from Configuration.Generator.PythiaUEZ2Settings_cfi import *


process.generator = cms.EDFilter('Pythia6GeneratorFilter',
	comEnergy = cms.double(8000.0),
	crossSection = cms.untracked.double(1.562933e+05),
	filterEfficiency = cms.untracked.double(1),
	maxEventsToPrint = cms.untracked.int32(0),
	pythiaHepMCVerbosity = cms.untracked.bool(False),
	pythiaPylistVerbosity = cms.untracked.int32(0),

	PythiaParameters = cms.PSet(
		pythiaUESettingsBlock,
		processParameters = cms.vstring(
			"MSEL = 0"              #user defined process
			,"MSTP(1)=4"
			,"MSUB(1)=1"
			,"MDCY( 17,1)=0"
			,"PMAS(C17,1)=1000"
			,"CKIN(1)=600.0"
			,"KCHG(C17,1)=204"

		),
		parameterSets = cms.vstring(
			'pythiaUESettings',
			'processParameters',
		)
	)
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
	version = cms.untracked.string('\$Revision: 0.0 $'),
	name = cms.untracked.string('\$Source: src/Monopole/MonoGen/Pythia_monopole_gen_config.py'),
	annotation = cms.untracked.string('Summer2012-Z2star sample with PYTHIA6: Monopole MC (Cowden) TuneZ2star, 8 TeV COM')
)


# Path and EndPath definitions
process.generation_step = cms.Path(process.pgen)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step,process.RAWSIMoutput_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 

