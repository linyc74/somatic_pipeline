version 1.0

import "GeneralTask.wdl" as general
import "GenerateReadyBam.wdl" as generateBam
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
        Array[Array[File]] inFilePoNfastqs
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
        String libraryKit
        String ponName
        Array[String] ponSampleName
        Array[String] tumorSampleName
        Array[String] normalSampleName
        Array[String] finalOutputName
    }
 
    scatter (i in range(length(ponSampleName))) {
        Array[File] iFPFs = inFilePoNfastqs[i]
        String pSN = ponSampleName[i]

        call general.TrimGalore as trimPoNfastq {
            input:
                sampleName = pSN,
                inFileFastqR1 = iFPFs[0],
                inFileFastqR2 = iFPFs[1]
        }

        call generateBam.GenerateReadyBam as ponBam {
            input:
                inFileFastqs = trimPoNfastq.outFileFastqs,
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
                sampleName = pSN   
        }
    }

    call createPoN.CreateMutect2PoN {
        input:
            inFileNormalBams = ponBam.outFileBam,
            inFileNormalBamIndexs = ponBam.outFileBamIndex,
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            normalSampleName = normalSampleName,
            extraArgs = "--max-mnp-distance 0",
            ponName = ponName
    }

    scatter (i in range(length(finalOutputName))) {
        Array[File] iFTFs = inFileTumorFastqs[i]
        Array[File] iFNFs = inFileNormalFastqs[i]
        String tSN = tumorSampleName[i]
        String nSN = normalSampleName[i]
        String fON = finalOutputName[i]

        call mapper.TNpairedMapping as TNmapping {
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

        call caller.TNpairedVariantsCalling as variantCalling {
            input:
                inFileTumorBam = TNmapping.outFileTumorBam,
                inFileTumorBamIndex = TNmapping.outFileTumorBamIndex,
                inFileNormalBam = TNmapping.outFileNormalBam,
                inFileNormalBamIndex = TNmapping.outFileNormalBamIndex,
                inFileIntervalBed = inFileIntervalBed,
                inFileGermlineResource = inFileGermlineResource,
                inFileGermlineResourceIndex = inFileGermlineResourceIndex,
                inFilePON = CreateMutect2PoN.outFilePoNvcf,
                inFilePONindex = CreateMutect2PoN.outFilePoNvcfIndex,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                tumorSampleName = tSN,
                normalSampleName = nSN,
                sampleName = fON,
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
                sampleName = fON
    }
    }

    output {
        Array[Array[File]] outFileTumorFastqs = TNmapping.outFileTumorFastqs
        Array[Array[File]] outFileNormalFastqs = TNmapping.outFileNormalFastqs
        Array[Array[File]] outFileTumorFastqcHtmls = TNmapping.outFileTumorFastqcHtmls
        Array[Array[File]] outFileNormalFastqcHtmls = TNmapping.outFileNormalFastqcHtmls
        Array[Array[File]] outFileTumorFastqcZips = TNmapping.outFileTumorFastqcZips
        Array[Array[File]] outFileNormalFastqcZips = TNmapping.outFileNormalFastqcZips
        Array[File] outFileTumorBam = TNmapping.outFileTumorBam
        Array[File] outFileNormalBam = TNmapping.outFileNormalBam
        Array[File] outFileTumorBamIndex = TNmapping.outFileTumorBamIndex
        Array[File] outFileNormalBamIndex = TNmapping.outFileNormalBamIndex
        Array[File] outFileTumorRawBam = TNmapping.outFileTumorRawBam
        Array[File] outFileNormalRawBam = TNmapping.outFileNormalRawBam
        Array[File] outFileStatsTumorBam = TNmapping.outFileStatsTumorBam
        Array[File] outFileStatsNormalBam = TNmapping.outFileStatsNormalBam
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
        File outFilePON = CreateMutect2PoN.outFilePoNvcf
        File outFilePONindex = CreateMutect2PoN.outFilePoNvcfIndex
    }
}