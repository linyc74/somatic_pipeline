version 1.0

# WORKFLOW DEFINITION

# The annotation part within NYCU Dentistry somatic pipeline. To pick variants then annotate by PCGR.
workflow PickAndAnnotate {
    input {
        File inFileVcfSS
        File inFileVcfMU
        File inFileVcfM2
        File inFileVcfLF
        File inFileVcfVD
        File infileVcfVS
        File inDirPCGRref
        File refFa
        String sampleName
    }
 
    call PythonVariantPicking {
        input:
            inFileVcfSS = inFileVcfSS,
            inFileVcfMU = inFileVcfMU,
            inFileVcfM2 = inFileVcfM2,
            inFileVcfLF = inFileVcfLF,
            inFileVcfVD = inFileVcfVD,
            infileVcfVS = infileVcfVS,
            refFa = refFa,
            sampleName = sampleName
    }

    call PCGR {
        input:
            inFileVcfGz = PythonVariantPicking.outFileVcfGz,
            inFileVcfIndex = PythonVariantPicking.outFileVcfIndex,
            inDirPCGRref = inDirPCGRref,
            sampleName = sampleName
    }
 
    output {
        File outFilePCGRannotatedVcf = PCGR.outFileVcf
        File outFilePCGRannotatedVcfIndex = PCGR.outFileVcfIndex
        File outFileMaf = PCGR.outFileMaf
        File outFilePCGRflexdbHtml = PCGR.outFileFlexdbHtml
        File outFilePCGRhtml = PCGR.outFileHtml       
    }
}

# TASK DEFINITIONS

# Picking variants from multiple caller's vcf using the self maintained python code
task PythonVariantPicking {
    input {
        File inFileVcfSS
        File inFileVcfMU
        File inFileVcfM2
        File inFileVcfLF
        File inFileVcfVD
        File infileVcfVS
        File refFa
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        python /usr/local/seqslab/variant picking \
        --ref-fa ~{refFa} \
        --somatic-sniper ~{inFileVcfSS} \
        --muse ~{inFileVcfMU} \
        --mutect2 ~{inFileVcfM2} \
        --lofreq ~{inFileVcfLF} \
        --vardict ~{inFileVcfVD} \
        --varscan ~{infileVcfVS} \
        --output-vcf ~{sampleName}_picked.vcf \
        --min-snv-caller 4 \
        --min-indel-callers 2
        bgzip \
        --stdout ~{sampleName}_picked.vcf > ~{sampleName}_picked.vcf.gz
        tabix \
        --preset vcf \
        ~{sampleName}_picked.vcf.gz
    >>>
 
    output {
        File outFileVcfGz = "~{sampleName}_picked.vcf.gz"
        File outFileVcfIndex = "~{sampleName}_picked.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Annotate vcf using PCGR
task PCGR {
    input {
        File inFileVcfGz
        File inFileVcfIndex
        File inDirPCGRref
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        pcgr \
        --input_vcf ~{inFileVcfGz} \
        --pcgr_dir ~{inDirPCGRref} \
        --output_dir pcgr_output \
        --genome_assembly grch38 \
        --vep_buffer_size 30000 \
        --sample_id ~{sampleName}
        # --vcf2maf

        zcat pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz > file_for_vcf2maf.vcf
        vcf2maf.pl \
        --input-vcf file_for_vcf2maf.vcf \
        --output-maf ~{sampleName}.maf \
        --tumor-id ~{sampleName} \
        --inhibit-vep \
        --ncbi-build GRCh38 \
        --ref-fasta ~{inDirPCGRref}/data/grch38/.vep/homo_sapiens/105_GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
    >>>
 
    output {
        File outFileVcf = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz"
        File outFileVcfIndex = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz.tbi"
        File outFileMaf = "~{sampleName}.maf"
        File outFileFlexdbHtml = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.flexdb.html"
        File outFileHtml = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.html"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

