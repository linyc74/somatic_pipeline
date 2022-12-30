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

    call general.BgzipBcftoolsIndex as snpCompressIndex {
        input:
            inFileVcf = VarscanSomatic.outFileSnpVcf,
            sampleName = sampleName  
    }

    call general.BgzipBcftoolsIndex as indelCompressIndex {
        input:
            inFileVcf = VarscanSomatic.outFileIndelVcf,
            sampleName = sampleName  
    }

    call general.Concat as concat {
        input:
            inFileSnvVcf = snpCompressIndex.outFileVcfGz,
            inFileSnvVcfIndex = snpCompressIndex.outFileVcfIndex,
            inFileIndelVcf = indelCompressIndex.outFileVcfGz,
            infileIndelVcfIndex = indelCompressIndex.outFileVcfIndex,
            sampleName = sampleName
    }

    call general.BgzipTabix as compressVcf {
        input:
            inFileVcf = concat.outFileVcf,
            sampleName = sampleName
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = concat.outFileVcf,
            sampleName = sampleName
    }

    call general.BgzipTabix as compressPyVcf {
        input:
            inFileVcf = filter.outFileVcf,
            sampleName = sampleName
    }

    output {
        File outFileVcfGz = compressVcf.outFileVcfGz
        File outFileVcfIndex = compressVcf.outFileVcfIndex
        File outFilePythonFilterVcfGz = compressPyVcf.outFileVcfGz
        File outFilePythonFilterVcfIndex = compressPyVcf.outFileVcfIndex
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
    >>>
 
    output {
        File outFileSnpVcf = "~{sampleName}_snp.vcf"
        File outFileIndelVcf = "~{sampleName}_indel.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}