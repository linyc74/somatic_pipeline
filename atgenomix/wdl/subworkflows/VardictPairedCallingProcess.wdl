version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a VarDict paired variant calling processed ready vcf
workflow VardictPairedCallingProcess {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        File refFaGzi
        Float minimumAF = 0.01
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call VardictPaired {
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

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = VardictPaired.outFileVcf,
            sampleName = sampleName
    }    
 
    output {
        File outFileVcf = filter.outFileVcf
    }
}



# TASK DEFINITIONS

# Call variants using VarDict paired variant calling mode
task VardictPaired {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        File refFaGzi
        Float minimumAF = 0.01
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        vardict \
        -G ~{refFa} \
        -f ~{minimumAF} \
        -N ~{tumorSampleName} \
        -b "~{inFileTumorBam} | ~{inFileNormalBam}" \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 \
        ~{inFileIntervalBed} > test.txt \
        | \
        testsomatic.R \
        | \
        var2vcf_paired.pl \
        -N "~{tumorSampleName} | ~{normalSampleName}" \
        -f ~{minimumAF} \
        1 > ~{sampleName}.vcf
    >>>
 
    output {
        File outFileVcf = "~{sampleName}.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}
  