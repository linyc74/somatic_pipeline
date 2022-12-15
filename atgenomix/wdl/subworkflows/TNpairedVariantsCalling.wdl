version 1.0

import "VarscanSomaticCallingProcess.wdl" as varscanProcess

# WORKFLOW DEFINITION

# Call variants using multiple caller with tumor normal paired mode
workflow TNpairedVariantsCalling {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File refFa
        File refFai
        File refFaGzi
        Int bgzipIndexThreads
        Int concatThreads
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call varscanProcess.VarscanSomaticCallingProcess as varscan {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileNormalBam = inFileNormalBam,
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
        File outFileVarscanVcf = varscan.outFileVcf
    }
}



# TASK DEFINITIONS

#'task start here'