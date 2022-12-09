version 1.0

# TASK DEFINITIONS

# Fastq preprocessing using Trim Galore with paired-end option
task TrimGalore {
    input {
        File inFileFastqR1
        File inFileFastqR2
        Int cores = 2
        Int discardReadLength = 20
        Int maxNcount = 0
        Int trimQuality = 20
        String fastqcArg = "--threads 16"
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
        ~{inFileFastqR1} \
        ~{inFileFastqR2}
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