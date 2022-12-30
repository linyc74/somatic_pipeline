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
            inFileVcfSS = variantCalling.outFileBamsomaticsniperVcf,
            inFileVcfMU = variantCalling.outFileMuseVcf,
            inFileVcfM2 = variantCalling.outFileMutect2Vcf,
            inFileVcfLF = variantCalling.outFileLofreqVcf,
            inFileVcfVD = variantCalling.outFileVardictVcf,
            infileVcfVS = variantCalling.outFileVarscanVcf,
            inDirPCGRref = inDirPCGRref,
            refFa = refFa,
            sampleName = finalOutputName
    }
    
    output {
        File outFileTumorBam = generateBam.outFileTumorBam
        File outFileNormalBam = generateBam.outFileNormalBam
        File outFileTumorBamIndex = generateBam.outFileTumorBamIndex
        File outFileNormalBamIndex = generateBam.outFileNormalBamIndex
        File outFileTumorRawBam = generateBam.outFileTumorRawBam
        File outFileNormalRawBam = generateBam.outFileNormalRawBam
        File outFileVcf = vcfAnnotate.outFileVcf
        File outFileVcfIndex = vcfAnnotate.outFileVcfIndex
        File outFileMaf = vcfAnnotate.outFileMaf
        File outFileFlexdbHtml = vcfAnnotate.outFileFlexdbHtml
        File outFileHtml = vcfAnnotate.outFileHtml
    }
}

# TASK DEFINITIONS