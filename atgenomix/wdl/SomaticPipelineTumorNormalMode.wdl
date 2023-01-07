version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/TNpairedMapping.wdl" as mapper
import "subworkflows/TNpairedVariantsCalling.wdl" as caller
import "subworkflows/PickAndAnnotate.wdl" as annotate

# WORKFLOW DEFINITION

# NYCU Dentistry somatic pipeline in Tumor-Normal paired mode
workflow SomaticPipelineTumorNormalMode {
    input {
        Array[File] inFileTumorFastqs
        Array[File] inFileNormalFastqs
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileIntervalBed
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File inFilePON
        File inFilePONindex
        File inDirPCGRref
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String normalSampleName
        String finalOutputName
    }

    call general.FastQC as fastqcTumorFastq {
        input:
            inFileFastqR1 = inFileTumorFastqs[0],
            inFileFastqR2 = inFileTumorFastqs[1]
    }

    call general.FastQC as fastqcNormalFastq {
        input:
            inFileFastqR1 = inFileNormalFastqs[0],
            inFileFastqR2 = inFileNormalFastqs[1]
    }

    call mapper.TNpairedMapping as TNmapping {
        input:
            inFileTumorFastqs = inFileTumorFastqs,
            inFileNormalFastqs = inFileNormalFastqs,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName
    }

    call caller.TNpairedVariantsCalling as variantCalling {
        input:
            inFileTumorBam = TNmapping.outFileTumorBam,
            inFileTumorBamIndex = TNmapping.outFileTumorBamIndex,
            inFileNormalBam = TNmapping.outFileNormalBam,
            inFileNormalBamIndex = TNmapping.outFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inFilePON = inFilePON,
            inFilePONindex = inFilePONindex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = finalOutputName,
            vardictMinimumAF = 0.01
    }

    call annotate.PickAndAnnotate as vcfAnnotate {
        input:
            inFileVcfSS = variantCalling.outFileBamsomaticsniperPyVcfGz,
            inFileVcfMU = variantCalling.outFileMusePyVcfGz,
            inFileVcfM2 = variantCalling.outFileMutect2PyVcfGz,
            inFileVcfLF = variantCalling.outFileLofreqPyVcfGz,
            inFileVcfVD = variantCalling.outFileVardictPyVcfGz,
            infileVcfVS = variantCalling.outFileVarscanPyVcfGz,
            inDirPCGRref = inDirPCGRref,
            refFa = refFa,
            sampleName = finalOutputName
    }
    
    output {
        Array[File] outFileTumorFastqs = TNmapping.outFileTumorFastqs
        Array[File] outFileNormalFastqs = TNmapping.outFileNormalFastqs
        Array[File] outFileTumorFastqcHtmls = fastqcTumorFastq.outFileHtmls
        Array[File] outFileNormalFastqcHtmls = fastqcNormalFastq.outFileHtmls
        Array[File] outFileTumorFastqcZips = fastqcTumorFastq.outFileZips
        Array[File] outFileNormalFastqcZips = fastqcNormalFastq.outFileZips
        File outFileTumorBam = TNmapping.outFileTumorBam
        File outFileNormalBam = TNmapping.outFileNormalBam
        File outFileTumorBamIndex = TNmapping.outFileTumorBamIndex
        File outFileNormalBamIndex = TNmapping.outFileNormalBamIndex
        File outFileTumorRawBam = TNmapping.outFileTumorRawBam
        File outFileNormalRawBam = TNmapping.outFileNormalRawBam
        File outFileStatsTumorBam = TNmapping.outFileStatsTumorBam
        File outFileStatsNormalBam = TNmapping.outFileStatsNormalBam
        File outFileBamsomaticsniperPyVcfGz = variantCalling.outFileBamsomaticsniperPyVcfGz
        File outFileBamsomaticsniperPyVcfIndex = variantCalling.outFileBamsomaticsniperPyVcfIndex
        File outFileLofreqPyVcfGz = variantCalling.outFileLofreqPyVcfGz
        File outFileLofreqPyVcfIndex = variantCalling.outFileLofreqPyVcfIndex
        File outFileMusePyVcfGz = variantCalling.outFileMusePyVcfGz
        File outFileMusePyVcfIndex = variantCalling.outFileMusePyVcfIndex
        File outFileMutect2PyVcfGz = variantCalling.outFileMutect2PyVcfGz
        File outFileMutect2PyVcfIndex = variantCalling.outFileMutect2PyVcfIndex
        File outFileVardictPyVcfGz = variantCalling.outFileVardictPyVcfGz
        File outFileVardictPyVcfIndex = variantCalling.outFileVardictPyVcfIndex
        File outFileVarscanPyVcfGz = variantCalling.outFileVarscanPyVcfGz
        File outFileVarscanPyVcfIndex = variantCalling.outFileVarscanPyVcfIndex
        File outFileBamsomaticsniperVcfGz = variantCalling.outFileBamsomaticsniperVcfGz
        File outFileBamsomaticsniperVcfIndex = variantCalling.outFileBamsomaticsniperVcfIndex
        File outFileLofreqVcfGz = variantCalling.outFileLofreqVcfGz
        File outFileLofreqVcfIndex = variantCalling.outFileLofreqVcfIndex
        File outFileMuseVcfGz = variantCalling.outFileMuseVcfGz
        File outFileMuseVcfIndex = variantCalling.outFileMuseVcfIndex
        File outFileMutect2VcfGz = variantCalling.outFileMutect2VcfGz
        File outFileMutect2VcfIndex = variantCalling.outFileMutect2VcfIndex
        File outFileVardictVcfGz = variantCalling.outFileVardictVcfGz
        File outFileVardictVcfIndex = variantCalling.outFileVardictVcfIndex
        File outFileVarscanVcfGz = variantCalling.outFileVarscanVcfGz
        File outFileVarscanVcfIndex = variantCalling.outFileVarscanVcfIndex
        File outFilePCGRannotatedVcf = vcfAnnotate.outFilePCGRannotatedVcf
        File outFilePCGRannotatedVcfIndex = vcfAnnotate.outFilePCGRannotatedVcfIndex
        File outFileMaf = vcfAnnotate.outFileMaf
        File outFilePCGRflexdbHtml = vcfAnnotate.outFilePCGRflexdbHtml
        File outFilePCGRhtml = vcfAnnotate.outFilePCGRhtml
    }
}