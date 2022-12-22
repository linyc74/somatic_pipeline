version 1.0

# TASK DEFINITIONS

# Compress vcf file using bgzip and index using tabix 
task BgzipTabix {
    input {
        File inFileVcf
        Int threads = 16
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bgzip \
        --threads ~{threads} \
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
        Int threads = 16
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bgzip \
        --threads ~{threads} \
        --force \
        --stdout ~{inFileVcf} > ~{sampleName}.vcf.gz
        bcftools index \
        --force \
        --threads ~{threads} \
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
        Int threads = 16
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bcftools concat \
        --allow-overlaps \
        --threads ~{threads} \ 
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

# Generate text pileup output for a bam file using samtools
task Mpileup {
    input {
        File inFileBam
        File inFileIntervalBed
        File refFa
        File refFai
        File refFaGzi
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
        --output-vcf ~{sampleName}_filtered.vcf \
        --variant-flagging-criteria ~{flaggingCriteria}  \
        --variant-removal-flags ~{removalFlags}
    >>>
 
    output {
        File outFileVcf = "~{sampleName}_filtered.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Fastq preprocessing using Trim Galore with paired-end option
task TrimGalore {
    input {
        File inFileFastqR1_PAR
        File inFileFastqR2_PAR
        Int cores = 2
        Int discardReadLength = 20
        Int maxNcount = 0
        Int trimQuality = 20
        String fastqcArg = "\"--threads 16\""
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        trim_galore \
        --paired \
        --quality ~{trimQuality} \
        --phred33 \
        --cores ~{cores} \
        --fastqc_args ~{fastqcArg} \
        --illumina \
        --length ~{discardReadLength} \
        --max_n ~{maxNcount} \
        --trim-n \
        --gzip \
        --output_dir out \
        --basename ~{sampleName} \
        ~{inFileFastqR1_PAR} \
        ~{inFileFastqR2_PAR}
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileFastqR1 = "out/~{sampleName}_val_1.fq.gz"
        File outFileFastqR2 = "out/~{sampleName}_val_2.fq.gz"
        File outFileHtmlR1 = "out/~{sampleName}_val_1_fastqc.html"
        File outFileHtmlR2 = "out/~{sampleName}_val_2_fastqc.html"
        File outFileZipR1 = "out/~{sampleName}_val_1_fastqc.zip"
        File outFileZipR2 = "out/~{sampleName}_val_2_fastqc.zip"
    }
}

# Struct

struct BwaIndex {
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
}