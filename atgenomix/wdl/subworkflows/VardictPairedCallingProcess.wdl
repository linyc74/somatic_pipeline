version 1.0

import "GeneralTask.wdl" as general

# WORKFLOW DEFINITION

# Generate a VarDict paired variant calling processed ready vcf
workflow VardictPairedCallingProcess {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        Float minimumAF
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    call VardictPaired {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            minimumAF = minimumAF,
            tumorSampleName = tumorSampleName
    }

    call TestsomaticR {
        input:
            inFileVardictStep1 = VardictPaired.outFileResult
    }

    call Var2Vcf {
        input:
            inFileVardictStep2 = TestsomaticR.outFileResult,
            minimumAF = minimumAF,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call general.BgzipTabix as compressVcf {
        input:
            inFileVcf = Var2Vcf.outFileVcf,
            sampleName = sampleName
    }

    call general.PythonVariantFilter as filter {
        input:
            inFileVcf = Var2Vcf.outFileVcf,
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

# Call variants using VarDict paired variant calling mode step1:vardict
task VardictPaired {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File refFa
        File refFai
        Float minimumAF
        String tumorSampleName
    }
 
    command <<<
        set -e -o pipefail
        VarDict \
        -G ~{refFa} \
        -f ~{minimumAF} \
        -N ~{tumorSampleName} \
        -b "~{inFileTumorBam}|~{inFileNormalBam}" \
        -c 1 \
        -S 2 \
        -E 3 \
        -g 4 \
        ~{inFileIntervalBed} \
        1> VardictStep1
    >>>
 
    output {
        File outFileResult = "VardictStep1"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Call variants using VarDict paired variant calling mode step2:testsomatic.R
task TestsomaticR {
    input {
        File inFileVardictStep1
    }
 
    command <<<
        set -e -o pipefail
        cat ~{inFileVardictStep1} \
        | \
        testsomatic.R \
        1> VardictStep2 
    >>>
 
    output {
        File outFileResult = "VardictStep2"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Call variants using VarDict paired variant calling mode step3:var2vcf
task Var2Vcf {
    input {
        File inFileVardictStep2
        Float minimumAF
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
 
    command <<<
        set -e -o pipefail
        cat ~{inFileVardictStep2} \
        | \
        var2vcf_paired.pl \
        -N "~{tumorSampleName} | ~{normalSampleName}" \
        -f ~{minimumAF} \
        1> ~{sampleName}.vcf
    >>>
 
    output {
        File outFileVcf = "~{sampleName}.vcf"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Call variants using VarDict paired variant calling mode
# task VardictPaired {
#     input {
#         File inFileTumorBam
#         File inFileTumorBamIndex
#         File inFileNormalBam
#         File inFileNormalBamIndex
#         File inFileIntervalBed
#         File refFa
#         File refFai
#         Float minimumAF = 0.01
#         String tumorSampleName
#         String normalSampleName
#         String sampleName
#     }
 
#     command <<<
#         set -e -o pipefail
#         vardict \
#         -G ~{refFa} \
#         -f ~{minimumAF} \
#         -N ~{tumorSampleName} \
#         -b "~{inFileTumorBam} | ~{inFileNormalBam}" \
#         -c 1 \
#         -S 2 \
#         -E 3 \
#         -g 4 \
#         ~{inFileIntervalBed} \
#         | \
#         testsomatic.R \
#         | \
#         var2vcf_paired.pl \
#         -N "~{tumorSampleName} | ~{normalSampleName}" \
#         -f ~{minimumAF} \
#         1 > ~{sampleName}.vcf
#     >>>
 
#     output {
#         File outFileVcf = "~{sampleName}.vcf"
#     }
 
#     runtime {
#         docker: 'nycu:latest'
#     }
# }