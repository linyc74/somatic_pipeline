version 1.0

import "TrimAndMapping.wdl" as mapping
import "PostMappingProcess.wdl" as postMapping

# WORKFLOW DEFINITION

# Take both tumor and normal paired-end fastq files, using TrimGalore trim then using bwa-mem align to reference genome
workflow TNpairedMapping {
    input {
        Array[File] inFileTumorFastqs
        Array[File] inFileNormalFastqs
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
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
    }

    call mapping.TrimAndMapping as trimAndMapTumorFastq {
        input:
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName,
            inFileFastqs = inFileTumorFastqs
    }

    call mapping.TrimAndMapping as trimAndMapNormalFastq {
        input:
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = normalSampleName,
            inFileFastqs = inFileNormalFastqs
    }

    call postMapping.PostMappingProcess as postMappingTumorBam {
        input:
            inFileUnSortRawBam = trimAndMapTumorFastq.outFileUnSortRawBam,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            sampleName = tumorSampleName
    }

    call postMapping.PostMappingProcess as postMappingNormalBam {
        input:
            inFileUnSortRawBam = trimAndMapNormalFastq.outFileUnSortRawBam,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            sampleName = normalSampleName   
    }
 
    output {
        Array[File] outFileTumorTrimmedFastqs = trimAndMapTumorFastq.outFileTrimmedFastq
        Array[File] outFileNormalTrimmedFastqs = trimAndMapNormalFastq.outFileTrimmedFastq
        File outFileTumorBam = postMappingTumorBam.outFileBam
        File outFileNormalBam = postMappingNormalBam.outFileBam
        File outFileTumorBamIndex = postMappingTumorBam.outFileBamIndex
        File outFileNormalBamIndex = postMappingNormalBam.outFileBamIndex
        File outFileTumorSortedRawBam = postMappingTumorBam.outFileSortedRawBam
        File outFileNormalSortedRawBam = postMappingNormalBam.outFileSortedRawBam
        File outFileStatsTumorBam = postMappingTumorBam.outFileBamStats
        File outFileStatsNormalBam = postMappingNormalBam.outFileBamStats
    }
}