version 1.0

import "Mutect2CallingProcess.wdl" as mutect2Process

# WORKFLOW DEFINITION

# Create a Panel of Normal for Mutect2
workflow CreateMutect2PoN {
    input {
        Array[File] inFileNormalBams
        Array[File] inFileNormalBamIndexs
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File inFileIntervalBed
        File refFa
        File refFai
        File refDict
        Array[String] normalSampleName
        String extraArgs = "--max-mnp-distance 0"
        String ponName
    }

    scatter (i in range(length(normalSampleName))) {
        File iFNB = inFileNormalBams[i]
        File iFNBI = inFileNormalBamIndexs[i]
        String nSN = normalSampleName[i]

        call mutect2Process.Mutect2CallingProcess as M2ForPoN {
            input:
                inFileTumorBam = iFNB,
                inFileTumorBamIndex = iFNBI,
                inFileIntervalBed = inFileIntervalBed,
                refFa = refFa,
                refFai = refFai,
                refDict = refDict,
                tumorSampleName = nSN,
                sampleName = nSN,
                extraArgs = extraArgs
        }
    }

    call GenomicDBimport {
        input:
            inFileVcfs = M2ForPoN.outFileVcfGz,
            inFileVcfIndexs = M2ForPoN.outFileVcfIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict
    }

    call CreateSomaticPanelOfNormals {
        input:
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inDirPoNdb = GenomicDBimport.outDirPoNdb,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            ponName = ponName
    }
 
    output {
        File outFilePoNvcf = CreateSomaticPanelOfNormals.outFilePoNvcf
        File outFilePoNvcfIndex = CreateSomaticPanelOfNormals.outFilePoNvcfIndex
    }
}

# TASK DEFINITIONS

# Step 2 of GATK create a panel of normals: Create a GenomicsDB from the normal Mutect2 calls
task GenomicDBimport {
    input {
        Array[File] inFileVcfs
        Array[File] inFileVcfIndexs
        File inFileIntervalBed
        File refFa
        File refFai
        File refDict
    }
 
    command <<<
        set -e -o pipefail
        gatk GenomicsDBImport \
        -R ~{refFa} \
        -L ~{inFileIntervalBed} \
        --genomicsdb-workspace-path pon_db \
        -V ~{sep=' -V ' inFileVcfs} \
        --read-index ~{sep=' --read-index ' inFileVcfIndexs}
    >>>
 
    output {
        File outDirPoNdb = "pon_db"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}

# Step 3 of GATK create a panel of normals: Combine the normal calls using CreateSomaticPanelOfNormals
task CreateSomaticPanelOfNormals {
    input {
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File inDirPoNdb
        File refFa
        File refFai
        File refDict
        String ponName
    }
 
    command <<<
        set -e -o pipefail
        gatk CreateSomaticPanelOfNormals \
        -R ~{refFa} \
        -V gendb://~{inDirPoNdb} \
        -O ~{ponName}.vcf \
        --germline-resource ~{inFileGermlineResource}
    >>>
 
    output {
        File outFilePoNvcf = "~{ponName}.vcf"
        File outFilePoNvcfIndex = "~{ponName}.vcf.idx"
    }
 
    runtime {
        docker: 'nycu:latest'
    }
}


