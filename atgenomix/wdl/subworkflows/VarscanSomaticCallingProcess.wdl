version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Call somatic variants using VarScan
workflow VarscanSomaticCallingProcess {
    input {
        File inFileTumorBam
        File inFileNormalBam
        File refFa
        File refFai
        File refFaGzi
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call general.Mpileup as tumorMpileup {
        input:
            inFileBam = inFileTumorBam,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
            sampleName = tumorSampleName
    }
 
    call general.Mpileup as normalMpileup {
        input:
            inFileBam = inFileNormalBam,
            refFa = refFa,
            refFai = refFai,
            refFaGzi = refFaGzi,
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

    output {
        File outFile = concat.outFileVcf
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
        varscan somatic \
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