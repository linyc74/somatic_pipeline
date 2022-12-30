version 1.0

import "subworkflows/TNpairedMapping.wdl" as mapper
import "subworkflows/TNpairedVariantsCalling.wdl" as caller
import "subworkflows/PickAndAnnotate.wdl" as annotate
import "subworkflows/CreateMutect2PoN.wdl" as createPoN

# WORKFLOW DEFINITION

# 'Workflow description'
workflow ScatterInProcessPoNSomaticpipelineTumorNormalMode {
    input {
        Array[Array[File]] inFileTumorFastqs
        Array[Array[File]] inFileNormalFastqs
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileIntervalBed
        File inFileGermlineResource
        File inFileGermlineResourceIndex
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
        String ponName
        String extraArgs = "--max-mnp-distance 0"
        Array[String] tumorSampleName
        Array[String] normalSampleName
        Array[String] finalOutputName
    }
 
    scatter (i in range(length(finalOutputName))) {
        Array[File] iFTFs = inFileTumorFastqs[i]
        Array[File] iFNFs = inFileNormalFastqs[i]
        String tSN = tumorSampleName[i]
        String nSN = normalSampleName[i]

        call mapper.TNpairedMapping as generateBam {
            input:
                inFileTumorFastqs = iFTFs,
                inFileNormalFastqs = iFNFs,
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
                tumorSampleName = tSN,
                normalSampleName = nSN
        }
    }

    call createPoN.CreateMutect2PoN {
        input:
            inFileNormalBams = generateBam.outFileNormalBam,
            inFileNormalBamIndexs = generateBam.outFileNormalBamIndex,
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            normalSampleName = normalSampleName,
            extraArgs = extraArgs,
            ponName = ponName
    }

    scatter (i in range(length(finalOutputName))) {
        File iFTB = generateBam.outFileTumorBam[i]
        File iFTBI = generateBam.outFileTumorBamIndex[i]
        File iFNB = generateBam.outFileNormalBam[i]
        File iFNBI = generateBam.outFileNormalBamIndex[i]
        String tSN2 = tumorSampleName[i]
        String nSN2 = normalSampleName[i]
        String fON = finalOutputName[i]

        call caller.TNpairedVariantsCalling as variantCalling {
            input:
                inFileTumorBam = iFTB,
                inFileTumorBamIndex = iFTBI,
                inFileNormalBam = iFNB,
                inFileNormalBamIndex = iFNBI,
                inFileIntervalBed = inFileIntervalBed,
                inFileGermlineResource = inFileGermlineResource,
                inFileGermlineResourceIndex = inFileGermlineResourceIndex,
                inFilePON = CreateMutect2PoN.outFilePoNvcf,
                inFilePONindex = CreateMutect2PoN.outFilePoNvcfIndex,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                tumorSampleName = tSN2,
                normalSampleName = nSN2,
                sampleName = fON,
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
                sampleName = fON
    }
    }

    output {
        Array[Array[File]] outFileTumorFastqs = generateBam.outFileTumorFastqs
        Array[Array[File]] outFileNormalFastqs = generateBam.outFileNormalFastqs
        Array[Array[File]] outFileTumorFastqcHtmls = generateBam.outFileTumorFastqcHtmls
        Array[Array[File]] outFileNormalFastqcHtmls = generateBam.outFileNormalFastqcHtmls
        Array[Array[File]] outFileTumorFastqcZips = generateBam.outFileTumorFastqcZips
        Array[Array[File]] outFileNormalFastqcZips = generateBam.outFileNormalFastqcZips
        Array[File] outFileTumorBam = generateBam.outFileTumorBam
        Array[File] outFileNormalBam = generateBam.outFileNormalBam
        Array[File] outFileTumorBamIndex = generateBam.outFileTumorBamIndex
        Array[File] outFileNormalBamIndex = generateBam.outFileNormalBamIndex
        Array[File] outFileTumorRawBam = generateBam.outFileTumorRawBam
        Array[File] outFileNormalRawBam = generateBam.outFileNormalRawBam
        Array[File] outFileStatsTumorBam = generateBam.outFileStatsTumorBam
        Array[File] outFileStatsNormalBam = generateBam.outFileStatsNormalBam
        Array[File] outFileBamsomaticsniperPyVcfGz = variantCalling.outFileBamsomaticsniperPyVcfGz
        Array[File] outFileBamsomaticsniperPyVcfIndex = variantCalling.outFileBamsomaticsniperPyVcfIndex
        Array[File] outFileLofreqPyVcfGz = variantCalling.outFileLofreqPyVcfGz
        Array[File] outFileLofreqPyVcfIndex = variantCalling.outFileLofreqPyVcfIndex
        Array[File] outFileMusePyVcfGz = variantCalling.outFileMusePyVcfGz
        Array[File] outFileMusePyVcfIndex = variantCalling.outFileMusePyVcfIndex
        Array[File] outFileMutect2PyVcfGz = variantCalling.outFileMutect2PyVcfGz
        Array[File] outFileMutect2PyVcfIndex = variantCalling.outFileMutect2PyVcfIndex
        Array[File] outFileVardictPyVcfGz = variantCalling.outFileVardictPyVcfGz
        Array[File] outFileVarscanPyVcfGz = variantCalling.outFileVarscanPyVcfGz
        Array[File] outFileVarscanPyVcfIndex = variantCalling.outFileVarscanPyVcfIndex
        Array[File] outFileBamsomaticsniperVcfGz = variantCalling.outFileBamsomaticsniperVcfGz
        Array[File] outFileBamsomaticsniperVcfIndex = variantCalling.outFileBamsomaticsniperVcfIndex
        Array[File] outFileLofreqVcfGz = variantCalling.outFileLofreqVcfGz
        Array[File] outFileLofreqVcfIndex = variantCalling.outFileLofreqVcfIndex
        Array[File] outFileMuseVcfGz = variantCalling.outFileMuseVcfGz
        Array[File] outFileMuseVcfIndex = variantCalling.outFileMuseVcfIndex
        Array[File] outFileMutect2VcfGz = variantCalling.outFileMutect2VcfGz
        Array[File] outFileMutect2VcfIndex = variantCalling.outFileMutect2VcfIndex
        Array[File] outFileVardictVcfGz = variantCalling.outFileVardictVcfGz
        Array[File] outFileVardictVcfIndex = variantCalling.outFileVardictVcfIndex
        Array[File] outFileVarscanVcfGz = variantCalling.outFileVarscanVcfGz
        Array[File] outFileVarscanVcfIndex = variantCalling.outFileVarscanVcfIndex
        Array[File] outFilePCGRannotatedVcf = vcfAnnotate.outFilePCGRannotatedVcf
        Array[File] outFilePCGRannotatedVcfIndex = vcfAnnotate.outFilePCGRannotatedVcfIndex
        Array[File] outFileMaf = vcfAnnotate.outFileMaf
        Array[File] outFilePCGRflexdbHtml = vcfAnnotate.outFilePCGRflexdbHtml
        Array[File] outFilePCGRhtml = vcfAnnotate.outFilePCGRhtml
    }
}