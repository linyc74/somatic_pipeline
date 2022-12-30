version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a MuSE processed ready vcf
workflow MuseCallingProcess {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File refFa
        File refFai
        String sampleName
    }
 
    call MuseCall {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            refFa = refFa,
            refFai = refFai,
            sampleName = sampleName
    }

    call MuseSump {
        input:
            inFileMuseResult = MuseCall.outFileMuseResult,
            sampleName = sampleName
    }

    call general.BgzipTabix as compressVcf {
        input:
            inFileVcf = MuseSump.outFileVcf,
            sampleName = sampleName
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = MuseSump.outFileVcf,
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

# Step 1 of MuSE: carries out pre-filtering and calculating position-specific summary statistics using the Markov substitution model.
task MuseCall {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File refFa
        File refFai
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        muse call \
        -f ~{refFa} \
        -O ~{sampleName} \
        ~{inFileTumorBam} \
        ~{inFileNormalBam}
    >>>
 
    output {
        File outFileMuseResult = "~{sampleName}.MuSE.txt"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Step 2 of MuSE: computes tier-based cutoffs from a sample-specific error model.
task MuseSump {
    input {
        File inFileMuseResult
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        muse sump \
        -I ~{inFileMuseResult} \
        -E \
        -O ~{sampleName}.vcf
    >>>
 
    output {
        File outFileVcf = "~{sampleName}.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}