version 1.0
 
import "GeneralTask.wdl" as general

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

    call general.BgzipTabix as vcfCompressIndex {
        input:
            inFileVcf = PythonVariantPicking.outFileVcf,
            sampleName = sampleName
    }

    call PCGR {
        input:
            inFileVcfGz = vcfCompressIndex.outFileVcfGz,
            inFileVcfIndex = vcfCompressIndex.outFileVcfIndex,
            inDirPCGRref = inDirPCGRref,
            sampleName = sampleName
    }
 
    output {
        File outFileVcf = PCGR.outFileVcf
        File outFileVcfIndex = PCGR.outFileVcfIndex
        File outFileMaf = PCGR.outFileMaf
        File outFileFlexdbHtml = PCGR.outFileFlexdbHtml
        File outFileHtml = PCGR.outFileHtml       
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
        --min-snv-caller 2 \
        --min-indel-callers 1
    >>>
 
    output {
        File outFileVcf = "~{sampleName}_picked.vcf"
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
        --sample_id ~{sampleName} \
        --vcf2maf
    >>>
 
    output {
        File outFileVcf = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz"
        File outFileVcfIndex = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.vcf.gz.tbi"
        File outFileMaf = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.maf"
        File outFileFlexdbHtml = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.flexdb.html"
        File outFileHtml = "pcgr_output/~{sampleName}.pcgr_acmg.grch38.html"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

