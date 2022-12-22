version 1.0

import "BamSomaticsniperCallingProcess.wdl" as bamSomaticsniperProcess
import "LoFreqSomaticCallingProcess.wdl" as lofreqProcess
import "MuseCallingProcess.wdl" as museProcess
import "Mutect2CallingProcess.wdl" as mutect2Process
import "VardictPairedCallingProcess.wdl" as vardictProcess
import "VarscanSomaticCallingProcess.wdl" as varscanProcess

# WORKFLOW DEFINITION

# Call variants using multiple caller with tumor normal paired mode
workflow TNpairedVariantsCalling {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File inFilePON
        File inFilePONindex
        File refFa
        File refFai
        File refFaGzi
        File refDict
        Float minimumAF
        Int bgzipIndexThreads
        Int concatThreads
        Int lofreqThreads
        Int m2HmmThreads
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
    
    call bamSomaticsniperProcess.BamSomaticsniperCallingProcess as bamsomaticsniper {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileNormalBam = inFileNormalBam,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call lofreqProcess.LoFreqSomaticCallingProcess as lofreq {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            concatThreads = concatThreads,
            lofreqThreads = lofreqThreads,
            sampleName = sampleName
    }

    call museProcess.MuseCallingProcess as muse {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            sampleName = sampleName
    }

    call mutect2Process.Mutect2CallingProcess as mutect2 {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFilePON = inFilePON,
            inFilePONindex = inFilePONindex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            refDict = refDict,
            m2HmmThreads = m2HmmThreads,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call vardictProcess.VardictPairedCallingProcess as vardict {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            minimumAF = minimumAF,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }
 
    call varscanProcess.VarscanSomaticCallingProcess as varscan {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileNormalBam = inFileNormalBam,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            bgzipIndexThreads = bgzipIndexThreads,
            concatThreads = concatThreads,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }
 
    output {
        File outFileBamsomaticsniperVcf = bamsomaticsniper.outFileVcf
        File outFileLofreqVcf = lofreq.outFileVcf
        File outFileMuseVcf = muse.outFileVcf
        File outFileMutect2Vcf = mutect2.outFileVcf
        File outFileMutect2FilterStats = mutect2.outFileFilterStats
        File outFileVardictVcf = vardict.outFileVcf
        File outFileVarscanVcf = varscan.outFileVcf
    }
}



# TASK DEFINITIONS

#'task start here'