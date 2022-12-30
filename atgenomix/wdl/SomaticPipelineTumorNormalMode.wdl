version 1.0

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
        Float vardictMinimumAF = 0.01
        String libraryKit
        String tumorSampleName
        String normalSampleName
        String finalOutputName
    }

    call mapper.TNpairedMapping as generateBam {
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
            libraryKit = libraryKit,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName
    }

    call caller.TNpairedVariantsCalling as variantCalling {
        input:
            inFileTumorBam = generateBam.outFileTumorBam,
            inFileTumorBamIndex = generateBam.outFileTumorBamIndex,
            inFileNormalBam = generateBam.outFileNormalBam,
            inFileNormalBamIndex = generateBam.outFileNormalBamIndex,
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
            vardictMinimumAF = vardictMinimumAF
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
        Array[File] outFileTumorFastqs = generateBam.outFileTumorFastqs
        Array[File] outFileNormalFastqs = generateBam.outFileNormalFastqs
        Array[File] outFileTumorFastqcHtmls = generateBam.outFileTumorFastqcHtmls
        Array[File] outFileNormalFastqcHtmls = generateBam.outFileNormalFastqcHtmls
        Array[File] outFileTumorFastqcZips = generateBam.outFileTumorFastqcZips
        Array[File] outFileNormalFastqcZips = generateBam.outFileNormalFastqcZips
        File outFileTumorBam = generateBam.outFileTumorBam
        File outFileNormalBam = generateBam.outFileNormalBam
        File outFileTumorBamIndex = generateBam.outFileTumorBamIndex
        File outFileNormalBamIndex = generateBam.outFileNormalBamIndex
        File outFileTumorRawBam = generateBam.outFileTumorRawBam
        File outFileNormalRawBam = generateBam.outFileNormalRawBam
        File outFileStatsTumorBam = generateBam.outFileStatsTumorBam
        File outFileStatsNormalBam = generateBam.outFileStatsNormalBam
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
        File outFileAnnotatedVcf = vcfAnnotate.outFileVcf
        File outFileAnnotatedVcfIndex = vcfAnnotate.outFileVcfIndex
        File outFileMaf = vcfAnnotate.outFileMaf
        File outFilePCGRflexdbHtml = vcfAnnotate.outFileFlexdbHtml
        File outFilePCGRhtml = vcfAnnotate.outFileHtml
    }
}