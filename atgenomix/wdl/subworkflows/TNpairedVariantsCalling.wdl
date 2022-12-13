version 1.0

import "VarscanSomaticCallingProcess.wdl" as varscan

# WORKFLOW DEFINITION

# Call variants using multiple caller with tumor normal paired mode
workflow TNpairedVariantsCalling {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File refFa
        File refFai
        File refFaGzi
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call varscan {
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
 
#    output {
#        File outFile = TaskName.outFile
#    }
}



# TASK DEFINITIONS

'task start here'