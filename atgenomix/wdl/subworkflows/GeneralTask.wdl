version 1.0

# TASK DEFINITIONS

# Compress vcf file using bgzip and index using tabix 
task BgzipTabix {
    input {
        File inFileVcf
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bgzip \
        --stdout ~{inFileVcf} > ~{sampleName}.vcf.gz
        tabix \
        --preset vcf \
        ~{sampleName}.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}.vcf.gz"
        File outFileVcfIndex = "~{sampleName}.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Compress vcf file using bgzip and index using bcftools 
task BgzipBcftoolsIndex {
    input {
        File inFileVcf
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bgzip \
        --force \
        --stdout ~{inFileVcf} > ~{sampleName}.vcf.gz
        bcftools index \
        --force \
        ~{sampleName}.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}.vcf.gz"
        File outFileVcfIndex = "~{sampleName}.vcf.gz.csi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Concatenate SNV.vcf and INDEL.vcf using bcftools
task Concat {
    input {
        File inFileSnvVcf
        File inFileSnvVcfIndex
        File inFileIndelVcf
        File infileIndelVcfIndex
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bcftools concat \
        --allow-overlaps \
        --output-type v \
        --output ~{sampleName}.vcf \
        ~{inFileSnvVcf} \
        ~{inFileIndelVcf}
    >>>
 
    output {
        File outFileVcf = "~{sampleName}.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Generate FastQC report for both fastq_R1 and fastq_R2
task FastQC {
    input {
        File inFileFastqR1
        File inFileFastqR2
    }
 
    command <<<
        set -e -o pipefail
        mkdir fastqc_report
        fastqc ~{inFileFastqR1} ~{inFileFastqR2} -o fastqc_report
        tar cvf fastqc_report.tar fastqc_report
    >>>
 
    output {
        File outFileFastqcReportTar = "fastqc_report.tar"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Generate text pileup output for a bam file using samtools
task Mpileup {
    input {
        File inFileBam
        File inFileIntervalBed
        File refFa
        File refFai
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        samtools mpileup \
        --fasta-ref ~{refFa} \
        --positions ~{inFileIntervalBed} \
        --output ~{sampleName}_pileup.txt \
        ~{inFileBam}
    >>>
 
    output {
        File outFilePileup = "~{sampleName}_pileup.txt"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Generate a filtered vcf using the self maintained python code
task PythonVariantFilter {
    input {
        File inFileVcf
        String flaggingCriteria = "\"LOW_DP: DP<20, HIGH_MQ: MQ>=30\""
        String removalFlags = "panel_of_normal,LOW_DP"
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        python /usr/local/seqslab/variant filtering \
        --input-vcf ~{inFileVcf} \
        --output-vcf ~{sampleName}_Pyfiltered.vcf \
        --variant-flagging-criteria ~{flaggingCriteria}  \
        --variant-removal-flags ~{removalFlags}
    >>>
 
    output {
        File outFileVcf = "~{sampleName}_Pyfiltered.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}