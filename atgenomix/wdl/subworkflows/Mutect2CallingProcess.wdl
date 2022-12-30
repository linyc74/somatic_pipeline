version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a Mutect2 processed ready vcf
workflow Mutect2CallingProcess {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File? inFileNormalBam
        File? inFileNormalBamIndex
        File? inFileGermlineResource
        File? inFileGermlineResourceIndex
        File? inFilePON
        File? inFilePONindex
        File inFileIntervalBed
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String? normalSampleName
        String sampleName
        String? extraArgs
    }
 
    call Mutect2 {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileGermlineResource = inFileGermlineResource,
            inFilePON = inFilePON,
            inFilePONindex = inFilePONindex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName,
            extraArgs = extraArgs
    }

    call LearnReadOrientationModel {
        input:
            inFileF1R2 = Mutect2.outFileF1R2
    }

    call FilterMutectCalls {
        input:
            inFileArtifactPriors = LearnReadOrientationModel.outFileArtifactPriors,
            inFileVcf = Mutect2.outFileVcf,
            inFileVcfStats = Mutect2.outFileVcfStats,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,   
            sampleName = sampleName     
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = FilterMutectCalls.outFileVcf,
            sampleName = sampleName
    }

    output {
        File outFileM2filterVcf = FilterMutectCalls.outFileVcf
        File outFileVcf = filter.outFileVcf
    }
}

# TASK DEFINITIONS

# Call variants using GATK Mutect2
task Mutect2 {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File? inFileNormalBam
        File? inFileNormalBamIndex
        File? inFileGermlineResource
        File? inFileGermlineResourceIndex
        File? inFilePON
        File? inFilePONindex
        File inFileIntervalBed
        File refFa
        File refFai
        File refDict
        String tumorSampleName
        String? normalSampleName
        String sampleName
        String? extraArgs
    }
 
    command <<<
        set -e -o pipefail
        gatk Mutect2 \
        --reference ~{refFa} \
        --intervals ~{inFileIntervalBed} \
        --input ~{inFileTumorBam} \
        ~{"--input " + inFileNormalBam} \
        --tumor-sample ~{tumorSampleName} \
        ~{"--normal-sample " + normalSampleName} \
        --output ~{sampleName}.vcf \
        --f1r2-tar-gz ~{sampleName}_f1r2.tar.gz \
        --max-reads-per-alignment-start 0 \
        ~{"--germline-resource " + inFileGermlineResource} \
        ~{"--panel-of-normals " + inFilePON} \
        ~{extraArgs}
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileVcf = "~{sampleName}.vcf"
        File outFileF1R2 = "~{sampleName}_f1r2.tar.gz"
        File outFileVcfStats =  "~{sampleName}.vcf.stats"
    }
}

# Get the maximum likelihood estimates of artifact prior probabilities in the orientation bias mixture model filter
task LearnReadOrientationModel {
    input {
        File inFileF1R2
    }
 
    command <<<
        set -e -o pipefail
        gatk LearnReadOrientationModel \
        --input ~{inFileF1R2} \
        --output artifact-prior-table.tar.gz
    >>>
 
    output {
        File outFileArtifactPriors = "artifact-prior-table.tar.gz"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Filter variants using GATK FilterMutectCalls
task FilterMutectCalls {
    input {
        File inFileArtifactPriors
        File inFileVcf
        File inFileVcfStats
        File refFa
        File refFai
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk FilterMutectCalls \
        --variant ~{inFileVcf} \
        --reference ~{refFa} \
        --output ~{sampleName}_filtered.vcf \
        --filtering-stats ~{inFileVcfStats} \
        --orientation-bias-artifact-priors ~{inFileArtifactPriors} \
        --create-output-variant-index false
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileVcf = "~{sampleName}_filtered.vcf"
    }
}