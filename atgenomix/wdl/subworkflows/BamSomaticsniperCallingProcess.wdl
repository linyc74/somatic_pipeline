version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a SomaticSniper processed ready vcf
workflow BamSomaticsniperCallingProcess {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File refFa
        File refFai
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call BamSomaticsniper {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileNormalBam = inFileNormalBam,
            refFa = refFa,
            refFai = refFai,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcfGz = BamSomaticsniper.outFileVcfGz,
            sampleName = sampleName,
            callerName = "somatic-sniper"
    }

    output {
        File outFileVcfGz = BamSomaticsniper.outFileVcfGz
        File outFileVcfIndex = BamSomaticsniper.outFileVcfIndex
        File outFilePythonFilterVcfGz = filter.outFileVcfGz
        File outFilePythonFilterVcfIndex = filter.outFileVcfIndex
    }
}

# TASK DEFINITIONS

# Call variants using SomaticSniper
task BamSomaticsniper {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File refFa
        File refFai
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bam-somaticsniper \
        -t ~{tumorSampleName} \
        -n ~{normalSampleName} \
        -F vcf \
        -f ~{refFa} \
        ~{inFileTumorBam} \
        ~{inFileNormalBam} \
        ~{sampleName}_somatic-sniper.vcf
        bgzip --stdout ~{sampleName}_somatic-sniper.vcf > ~{sampleName}_somatic-sniper.vcf.gz
        tabix --preset vcf ~{sampleName}_somatic-sniper.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}_somatic-sniper.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_somatic-sniper.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}