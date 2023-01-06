version 1.0

import "GeneralTask.wdl" as general
import "GenerateReadyBam.wdl" as generateBam

# WORKFLOW DEFINITION

# Take both tumor and normal paired-end fastq files, using TrimGalore trim and get fastqc report then using bwa-mem align to reference genome
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
 
    call general.TrimGalore as trimTumorFastq {
        input:
            sampleName = tumorSampleName,
            inFileFastqR1 = inFileTumorFastqs[0],
            inFileFastqR2 = inFileTumorFastqs[1]
    }

    call general.FastQC as fastqcTumorFastq {
        input:
            inFileFastqR1 = inFileTumorFastqs[0],
            inFileFastqR2 = inFileTumorFastqs[1]
    }

    call general.TrimGalore as trimNormalFastq {
        input:
            sampleName = normalSampleName,
            inFileFastqR1 = inFileNormalFastqs[0],
            inFileFastqR2 = inFileNormalFastqs[1]
    }

    call general.FastQC as fastqcNormalFastq {
        input:
            inFileFastqR1 = inFileNormalFastqs[0],
            inFileFastqR2 = inFileNormalFastqs[1]
    }

    call generateBam.GenerateReadyBam as tumorBam {
        input:
            inFileFastqs = trimTumorFastq.outFileFastqs,
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
            sampleName = tumorSampleName
    }

    call generateBam.GenerateReadyBam as normalBam {
        input:
            inFileFastqs = trimNormalFastq.outFileFastqs,
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
            sampleName = normalSampleName   
    }
 
    output {
        Array[File] outFileTumorFastqs = trimTumorFastq.outFileFastqs
        Array[File] outFileNormalFastqs = trimNormalFastq.outFileFastqs
        Array[File] outFileTumorFastqcHtmls = fastqcTumorFastq.outFileHtmls
        Array[File] outFileNormalFastqcHtmls = fastqcNormalFastq.outFileHtmls
        Array[File] outFileTumorFastqcZips = fastqcTumorFastq.outFileZips
        Array[File] outFileNormalFastqcZips = fastqcNormalFastq.outFileZips
        File outFileTumorBam = tumorBam.outFileBam
        File outFileNormalBam = normalBam.outFileBam
        File outFileTumorBamIndex = tumorBam.outFileBamIndex
        File outFileNormalBamIndex = normalBam.outFileBamIndex
        File outFileTumorRawBam = tumorBam.outFileBam
        File outFileNormalRawBam = normalBam.outFileBam
        File outFileStatsTumorBam = tumorBam.outFileBamStats
        File outFileStatsNormalBam = normalBam.outFileBamStats
    }
}