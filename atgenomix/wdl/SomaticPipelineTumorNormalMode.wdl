version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/GenerateReadyBam.wdl" as mapper
import "subworkflows/TNpairedVariantsCalling.wdl" as caller

# WORKFLOW DEFINITION

# NYCU Dentistry somatic pipeline in Tumor-Normal paired mode
workflow SomaticPipelineTumorNormalMode {
    input {
        File inFileTumorFastqR1
        File inFileTumorFastqR2
        File inFileNormalFastqR1
        File inFileNormalFastqR2
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File refFa
        File refFai
        File refDict
        File refFaGzi
        String libraryKit = 'NA'
        String tumorSampleName = 'tumor'
        String normalSampleName = 'normal'
        String sampleName
        BwaIndex bwaIndex
    }
 
    call general.TrimGalore as trimTumorFastq {
        input:
            sampleName = tumorSampleName,
            inFileFastqR1_PAR = inFileTumorFastqR1,
            inFileFastqR2_PAR = inFileTumorFastqR2
    }

    call general.TrimGalore as trimNormalFastq {
        input:
            sampleName = normalSampleName,
            inFileFastqR1_PAR = inFileNormalFastqR1,
            inFileFastqR2_PAR = inFileNormalFastqR2            
    }    

    call mapper.GenerateReadyBam as tumorBam {
        input:
            inFileFastqR1 = trimTumorFastq.outFileFastqR1,
            inFileFastqR2 = trimTumorFastq.outFileFastqR2,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            refFaGzi = refFaGzi,
            libraryKit = libraryKit,
            sampleName = tumorSampleName,
            bwaIndex = bwaIndex
    }

    call mapper.GenerateReadyBam as normalBam {
        input:
            inFileFastqR1 = trimNormalFastq.outFileFastqR1,
            inFileFastqR2 = trimNormalFastq.outFileFastqR2,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            refFaGzi = refFaGzi,
            libraryKit = libraryKit,
            sampleName = normalSampleName,
            bwaIndex = bwaIndex     
    }

    call caller.TNpairedVariantsCalling as callVcf {
        input:
            inFileTumorBam = tumorBam.outputBam,
            inFileNormalBam = normalBam.outputBam,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }
#    output {
#        File outputFile = TaskName.outputFile
#    }
}



# TASK DEFINITIONS