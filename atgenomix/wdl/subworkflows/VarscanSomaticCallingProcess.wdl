version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a VarScan processed ready vcf
workflow VarscanSomaticCallingProcess {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File inFileIntervalBed
        File refFa
        File refFai
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call general.Mpileup as tumorMpileup {
        input:
            inFileBam = inFileTumorBam,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = tumorSampleName
    }
 
    call general.Mpileup as normalMpileup {
        input:
            inFileBam = inFileNormalBam,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = normalSampleName
    }

    call VarscanSomatic {
        input:
            inFileTumorPileup = tumorMpileup.outFilePileup,
            inFileNormalPileup = normalMpileup.outFilePileup,
            sampleName = sampleName
    }

    call general.Concat as concat {
        input:
            inFileSnvVcf = VarscanSomatic.outFileSnpVcfGz,
            inFileSnvVcfIndex = VarscanSomatic.outFileSnpVcfIndex,
            inFileIndelVcf = VarscanSomatic.outFileIndelVcfGz,
            infileIndelVcfIndex = VarscanSomatic.outFileIndelVcfIndex,
            sampleName = sampleName
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcfGz = concat.outFileVcfGz,
            sampleName = sampleName
    }

    output {
        File outFileVcfGz = concat.outFileVcfGz
        File outFileVcfIndex = concat.outFileVcfIndex
        File outFilePythonFilterVcfGz = filter.outFileVcfGz
        File outFilePythonFilterVcfIndex = filter.outFileVcfIndex
    }
}


# TASK DEFINITIONS

# Calling somatic variants using VarScan
task VarscanSomatic {
    input {
        File inFileTumorPileup
        File inFileNormalPileup
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        java -jar /usr/local/seqslab/VarScan.v2.3.7.jar  somatic \
        ~{inFileNormalPileup} \
        ~{inFileTumorPileup} \
        --output-snp ~{sampleName}_snp.vcf \
        --output-indel ~{sampleName}_indel.vcf \
        --strand-filter 1 \
        --output-vcf 1
        bgzip --stdout ~{sampleName}_snp.vcf > ~{sampleName}_snp.vcf.gz
        bgzip --stdout ~{sampleName}_indel.vcf > ~{sampleName}_indel.vcf.gz        
        tabix --preset vcf ~{sampleName}_snp.vcf.gz
        tabix --preset vcf ~{sampleName}_indel.vcf.gz
    >>>
 
    output {
        File outFileSnpVcfGz = "~{sampleName}_snp.vcf.gz"
        File outFileIndelVcfGz = "~{sampleName}_indel.vcf.gz"
        File outFileSnpVcfIndex = "~{sampleName}_snp.vcf.gz.tbi"
        File outFileIndelVcfIndex = "~{sampleName}_indel.vcf.gz.tbi"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}