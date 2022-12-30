version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a LoFreq somatic processed ready vcf
workflow LoFreqSomaticCallingProcess {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        String sampleName
    }
 
    call LoFreqSomatic {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = sampleName
    }

    call general.Concat as concat {
        input:
            inFileSnvVcf = LoFreqSomatic.outFileSnvVcf,
            inFileSnvVcfIndex = LoFreqSomatic.outFileSnvVcfIndex,
            inFileIndelVcf = LoFreqSomatic.outFileIndelVcf,
            infileIndelVcfIndex = LoFreqSomatic.outFileIndelVcfIndex,
            sampleName = sampleName
    }

    call general.BgzipTabix as compressVcf {
        input:
            inFileVcf = concat.outFileVcf,
            sampleName = sampleName
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = concat.outFileVcf,
            sampleName = sampleName
    }

    call general.BgzipTabix as compressPyVcf {
        input:
            inFileVcf = filter.outFileVcf,
            sampleName = sampleName
    }
   
    output {
        File outFileVcfGz = compressVcf.outFileVcfGz
        File outFileVcfIndex = compressVcf.outFileVcfIndex
        File outFilePythonFilterVcfGz = compressPyVcf.outFileVcfGz
        File outFilePythonFilterVcfIndex = compressPyVcf.outFileVcfIndex
    }
}

# TASK DEFINITIONS

# Call somatic variants in matched tumor/normal pairs using Lofreq
task LoFreqSomatic {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        lofreq somatic \
        --normal ~{inFileTumorBam} \
        --tumor ~{inFileNormalBam} \
        --ref ~{refFa} \
        --threads 1 \
        --call-indels \
        -l ~{inFileIntervalBed} \
        -o ~{sampleName}
    >>>
 
    output {
        File outFileSnvVcf = "~{sampleName}somatic_final.snvs.vcf.gz"
        File outFileSnvVcfIndex = "~{sampleName}somatic_final.snvs.vcf.gz.tbi"
        File outFileIndelVcf = "~{sampleName}somatic_final.indels.vcf.gz"
        File outFileIndelVcfIndex = "~{sampleName}somatic_final.indels.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}