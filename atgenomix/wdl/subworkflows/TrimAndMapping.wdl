version 1.0

# WORKFLOW DEFINITION

# Take in both R1 and R2 fastq files then generate the trimmed fastq files and bam file.
workflow TrimAndMapping {
    input {
        Array[File] inFileFastqs
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        String sampleName
    }

    call TrimGalore {
        input:
            sampleName = sampleName,
            inFileFastqR1 = inFileFastqs[0],
            inFileFastqR2 = inFileFastqs[1]
    }

    call BwaMem {
        input:
            inFileFastqR1 = TrimGalore.outFileTrimmedFastqs[0],
            inFileFastqR2 = TrimGalore.outFileTrimmedFastqs[1],
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            sampleName = sampleName
    }
 
    output {
        Array[File] outFileTrimmedFastq = TrimGalore.outFileTrimmedFastqs
        File outFileUnSortRawBam = BwaMem.outFileUnSortRawBam
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
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        bwa mem \
        -t 4 \
        -R "@RG\tID:~{sampleName}\tSM:~{sampleName}\tPL:ILLUMINA\tLB:~{sampleName}" \
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
        File outFileUnSortRawBam = "~{sampleName}.bam"
    }
}

# Fastq preprocessing using Trim Galore with paired-end option
task TrimGalore {
    input {
        File inFileFastqR1
        File inFileFastqR2
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        trim_galore \
        --paired \
        --quality 20 \
        --phred33 \
        --illumina \
        --length 20 \
        --max_n 0 \
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
        Array[File] outFileTrimmedFastqs = ["out/~{sampleName}_val_1.fq.gz","out/~{sampleName}_val_2.fq.gz"]
        # File outFileHtmlR1 = "out/~{sampleName}_val_1_fastqc.html"
        # File outFileHtmlR2 = "out/~{sampleName}_val_2_fastqc.html"
        # File outFileZipR1 = "out/~{sampleName}_val_1_fastqc.zip"
        # File outFileZipR2 = "out/~{sampleName}_val_2_fastqc.zip"
    }
}
