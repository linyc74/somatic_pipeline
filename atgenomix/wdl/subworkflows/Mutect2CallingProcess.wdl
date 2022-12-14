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
        File inFilePON
        File refFa
        File refFai
        File refFaGzi
        File refDict
        String tumorSampleName
        String? normalSampleName
        String sampleName
    }
 
    call Mutect2 {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFilePON = inFilePON,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call LearnReadOrientationModel {
        input:
            inFileF1R2 = Mutect2.outFileF1R2
    }

    call FilterMutectCalls {
        input:
            inFileVcf = Mutect2.outFileVcf,
            inFileArtifactPriors = LearnReadOrientationModel.outFileArtifactPriors,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            refDict = refDict,   
            sampleName = sampleName     
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = FilterMutectCalls.outFileVcf,
            sampleName = sampleName
    }

    output {
        File outFileVcf = filter.outFileVcf
        File outFileFilterStats = FilterMutectCalls.outFileFilterStats
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
        File inFilePON
        File inFilePONindex
        File refFa
        File refFai
        File refFaGzi
        File refDict
        Int threads = 16
        String tumorSampleName
        String? normalSampleName
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk Mutect2 \
        --reference ~{refFa} \
        --input ~{inFileTumorBam} \
        ~{"--input " + inFileNormalBam} \
        --tumor-sample ~{tumorSampleName} \
        ~{"--normal-sample " + normalSampleName} \
        --output ~{sampleName}.vcf \
        --native-pair-hmm-threads ~{threads} \
        --f1r2-tar-gz ~{sampleName}_f1r2.tar.gz \
        --panel-of-normals ~{inFilePON}
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileVcf = "~{sampleName}.vcf"
        File outFileF1R2 = "~{sampleName}_f1r2.tar.gz"
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
        File inFileVcf
        File inFileArtifactPriors
        File refFa
        File refFai
        File refFaGzi
        File refDict
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        gatk FilterMutectCalls \
        --variant ~{inFileVcf} \
        --reference ~{refFa} \
        --ouptput ~{sampleName}_filtered.vcf \
        --filtering-stats ~{sampleName}_filter-stats.tsv \
        --orientation-bias-artifact-priors ~{inFileArtifactPriors} \
        --create-output-variant-index false
    >>>
 
    runtime {
        docker: 'nycu:latest'
    }
 
    output {
        File outFileVcf = "~{sampleName}_filtered.vcf"
        File outFileFilterStats = "~{sampleName}_filter-stats.tsv"
    }
}