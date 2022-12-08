version 1.0

# WORKFLOW DEFINITION

# Generate a analysis-ready bam file
workflow GenerateReadyBam {
    input {
        String sampleName
        File refFasta
        File refFastaFai
        File refFastaDict
        File refFastaGzi
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
    }

    call BwaMem {
        input:
            sampleName = sampleName,
            refFasta = refFasta,
            refFastaFai = refFastaFai,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa
    }

    call Sort { 
        input:
            inputBam = BwaMem.outputBam,
            sampleName=sampleName
    }

    call MarkDuplicates {
        input:
            inputBam = Sort.outputBam,
            sampleName = sampleName
    }

    call BaseRecalibrator {
        input:
            inputBam = MarkDuplicates.outputBam,
            sampleName = sampleName,
            refFasta = refFasta,
            refFastaFai = refFastaFai,
            refFastaDict = refFastaDict,
            refFastaGzi = refFastaGzi
    }

    call ApplyBqsr {
        input:
            inputBam = MarkDuplicates.outputBam,
            recalibrationTable = BaseRecalibrator.outputRecalibrationTable,
            sampleName = sampleName,
            refFasta = refFasta,
            refFastaFai = refFastaFai,
            refFastaDict = refFastaDict,
            refFastaGzi = refFastaGzi
    }

    call BamStats {
        input:
            inputBam = ApplyBqsr.outputBam,
            sampleName = sampleName
    }
}


# TASK DEFINITIONS

# Align reads using bwa mem and output a bam file
task BwaMem {
    input {
        File inputFastqR1
        File inputFastqR2
        File refFasta
        File refFastaFai
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        Int threads = 2
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bwa mem \
        -t ~{threads} \
        -R "@RG\tID:~{sampleName}\tSM:~{sampleName}\tPL:ILLUMINA\tLB:~{sampleName}" \
        ~{refFasta} \
        ~{inputFastqR1} \
        ~{inputFastqR2} \
        | \
        samtools view -b - > ~{sampleName}.bam  
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outputBam = "~{sampleName}.bam"
    }
}

# Sort reads within bam using samtools
task Sort {
    input {
        File inputBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools sort ~{inputBam} > ~{sampleName}.bam
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outputBam = "~{sampleName}.bam"
    }
}

# Mark duplicate reads using GATK
task MarkDuplicates {
    input {
        File inputBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk MarkDuplicates \
        --INPUT ~{inputBam} \
        --OUTPUT ~{sampleName}.bam \
        --METRICS_FILE ~{sampleName}_metrics.txt \
        --REMOVE_DUPLICATES false
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outputBam = "~{sampleName}.bam"
        File outputMetrics = "~{sampleName}_metrics.txt"
    }
}

# GATK BQSR step1. Generates recalibration table for Base Quality Score Recalibration
task BaseRecalibrator {
    input {
        File dbsnpVcf
        File dbsnpVcfIndex
        File inputBam
        File refFasta
        File refFastaFai
        File refFastaDict
        File refFastaGzi
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk BaseRecalibrator \
        --input ~{inputBam} \
        --reference ~{refFasta} \
        --known-sites ~{dbsnpVcf} \
        --output ~{sampleName}.recalibration.table
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outputRecalibrationTable = '~{sampleName}.recalibration.table'
    }
}

# GATK BQSR step2. Apply base quality score recalibration
task ApplyBqsr {
    input {
        File inputBam
        File recalibrationTable
        File refFasta
        File refFastaFai
        File refFastaDict
        File refFastaGzi
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk ApplyBQSR \
        --input ~{inputBam} \
        --reference ~{refFasta} \
        --bqsr-recal-file ~{recalibrationTable} \
        --output ~{sampleName}.bam

    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outputBam = "~{sampleName}.bam"
        File outputBamIndex = "~{sampleName}.bai"
    }
}

# Generate a comprehensive statistics report from bam file using samtools
task BamStats {
    input {
        File inputBam
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools stats ~{inputBam} > ~{sampleName}_stats.txt
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File bamStats = "~{sampleName}_stats.txt"
    }
}