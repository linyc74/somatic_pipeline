version 1.0

# WORKFLOW DEFINITION

# Generate a analysis-ready bam file and a comprehensive statistics report
workflow GenerateReadyBam {
    input {
        Array[File] inFileFastqs
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
        String libraryKit
        String sampleName
    }
    
    call BwaMem {
        input:
            inFileFastqR1 = inFileFastqs[0],
            inFileFastqR2 = inFileFastqs[1],
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            libraryKit = libraryKit,
            sampleName = sampleName
    }

    call Sort { 
        input:
            inFileBam = BwaMem.outFileBam,
            sampleName = sampleName
    }

    call MarkDuplicates {
        input:
            inFileBam = Sort.outFileBam,
            sampleName = sampleName
    }

    call BaseRecalibrator {
        input:
            inFileBam = MarkDuplicates.outFileBam,
            sampleName = sampleName,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict
    }

    call ApplyBqsr {
        input:
            inFileBam = MarkDuplicates.outFileBam,
            inFileRecalibrationTable = BaseRecalibrator.outFileRecalibrationTable,
            sampleName = sampleName,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict
    }

    call BamStats {
        input:
            inFileBam = ApplyBqsr.outFileBam,
            sampleName = sampleName
    }

    output {
        File outFileBam = ApplyBqsr.outFileBam
        File outFileBamIndex = ApplyBqsr.outFileBamIndex
        File outFileRawBam = Sort.outFileBam
        File outFileBamStats = BamStats.outFileBamStats
    }
}


# TASK DEFINITIONS

# Align reads using bwa mem and output a bam file
task BwaMem {
    input {
        File inFileFastqR1
        File inFileFastqR2
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        String libraryKit
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bwa mem \
        -t 4
        -R "@RG\tID:~{sampleName}\tSM:~{sampleName}\tPL:ILLUMINA\tLB:~{libraryKit}" \
        ~{refFa} \
        ~{inFileFastqR1} \
        ~{inFileFastqR2} \
        | \
        samtools view -b - > ~{sampleName}.bam  
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
    }
}

# Sort reads within bam using samtools
task Sort {
    input {
        File inFileBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools sort ~{inFileBam} > ~{sampleName}.bam
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
    }
}

# Mark duplicate reads using GATK
task MarkDuplicates {
    input {
        File inFileBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk MarkDuplicates \
        --INPUT ~{inFileBam} \
        --OUTPUT ~{sampleName}.bam \
        --METRICS_FILE ~{sampleName}_metrics.txt \
        --REMOVE_DUPLICATES false
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
        File outFileMetrics = "~{sampleName}_metrics.txt"
    }
}

# GATK BQSR step1. Generates recalibration table for Base Quality Score Recalibration
task BaseRecalibrator {
    input {
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileBam
        File refFa
        File refFai
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk BaseRecalibrator \
        --input ~{inFileBam} \
        --reference ~{refFa} \
        --known-sites ~{inFileDbsnpVcf} \
        --output ~{sampleName}.recalibration.table
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileRecalibrationTable = '~{sampleName}.recalibration.table'
    }
}

# GATK BQSR step2. Apply base quality score recalibration
task ApplyBqsr {
    input {
        File inFileBam
        File inFileRecalibrationTable
        File refFa
        File refFai
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk ApplyBQSR \
        --input ~{inFileBam} \
        --reference ~{refFa} \
        --bqsr-recal-file ~{inFileRecalibrationTable} \
        --output ~{sampleName}.bam

    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBam = "~{sampleName}.bam"
        File outFileBamIndex = "~{sampleName}.bai"
    }
}

# Generate a comprehensive statistics report from bam file using samtools
task BamStats {
    input {
        File inFileBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools stats ~{inFileBam} > ~{sampleName}_stats.txt
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileBamStats = "~{sampleName}_stats.txt"
    }
}