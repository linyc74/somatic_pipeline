version 1.0

# TASK DEFINITIONS

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
        bgzip --stdout ~{sampleName}.vcf > ~{sampleName}.vcf.gz
        tabix --preset vcf ~{sampleName}.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}.vcf.gz"
        File outFileVcfIndex = "~{sampleName}.vcf.gz.tbi"
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
        mkdir -p fastqc_report
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
        File inFileVcfGz
        String flaggingCriteria = "\"LOW_DP: DP<20, HIGH_MQ: MQ>=30\""
        String removalFlags = "panel_of_normal,LOW_DP"
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        zcat ~{inFileVcfGz} > in.vcf
        python /usr/local/seqslab/variant filtering \
        --input-vcf in.vcf \
        --output-vcf ~{sampleName}_Pyfiltered.vcf \
        --variant-flagging-criteria ~{flaggingCriteria}  \
        --variant-removal-flags ~{removalFlags}
        bgzip --stdout ~{sampleName}_Pyfiltered.vcf > ~{sampleName}_Pyfiltered.vcf.gz
        tabix --preset vcf ~{sampleName}_Pyfiltered.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}_Pyfiltered.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_Pyfiltered.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}