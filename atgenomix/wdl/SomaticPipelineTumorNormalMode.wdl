version 1.0

import "subworkflows/GeneralTask.wdl" as general
import "subworkflows/GenerateReadyBam.wdl" as mapper
import "subworkflows/TNpairedVariantsCalling.wdl" as caller
import "subworkflows/PickAndAnnotate.wdl" as annotate

# WORKFLOW DEFINITION

# NYCU Dentistry somatic pipeline in Tumor-Normal paired mode
workflow SomaticPipelineTumorNormalMode {
    input {
        File inFileTumorFastqR1
        File inFileTumorFastqR2
        File inFileNormalFastqR1
        File inFileNormalFastqR2
        File inFileDbsnpVcf
        File inFileDbsnpVcfIndex
        File inFileIntervalBed
        File inFilePON
        File inFilePONindex
        File inDirPCGRref
        File refAmb
        File refAnn
        File refBwt
        File refPac
        File refSa
        File refFa
        File refFai
        File refDict
        String libraryKit = 'NA'
        String tumorSampleName = 'tumor'
        String normalSampleName = 'normal'
        String finalOutputName
    }
 
    call general.TrimGalore as trimTumorFastq {
        input:
            sampleName = tumorSampleName,
            inFileFastqR1_PAR = inFileTumorFastqR1,
            inFileFastqR2_PAR = inFileTumorFastqR2
    }

    call general.TrimGalore as trimNormalFastq {
        input:
            sampleName = normalSampleName,
            inFileFastqR1_PAR = inFileNormalFastqR1,
            inFileFastqR2_PAR = inFileNormalFastqR2            
    }    

    call mapper.GenerateReadyBam as tumorBam {
        input:
            inFileFastqR1 = trimTumorFastq.outFileFastqR1,
            inFileFastqR2 = trimTumorFastq.outFileFastqR2,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            libraryKit = libraryKit,
            sampleName = tumorSampleName
    }

    call mapper.GenerateReadyBam as normalBam {
        input:
            inFileFastqR1 = trimNormalFastq.outFileFastqR1,
            inFileFastqR2 = trimNormalFastq.outFileFastqR2,
            inFileDbsnpVcf = inFileDbsnpVcf,
            inFileDbsnpVcfIndex = inFileDbsnpVcfIndex,
            refAmb = refAmb,
            refAnn = refAnn,
            refBwt = refBwt,
            refPac = refPac,
            refSa = refSa,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            libraryKit = libraryKit,
            sampleName = normalSampleName   
    }

    call caller.TNpairedVariantsCalling as variantCalling {
        input:
            inFileTumorBam = tumorBam.outFileBam,
            inFileTumorBamIndex = tumorBam.outFileBamIndex,
            inFileNormalBam = normalBam.outFileBam,
            inFileNormalBamIndex = normalBam.outFileBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            inFilePON = inFilePON,
            inFilePONindex = inFilePONindex,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = finalOutputName
    }

    call annotate.PickAndAnnotate as vcfAnnotate {
        input:
            inFileVcfSS = variantCalling.outFileBamsomaticsniperVcf,
            inFileVcfMU = variantCalling.outFileMuseVcf,
            inFileVcfM2 = variantCalling.outFileMutect2Vcf,
            inFileVcfLF = variantCalling.outFileLofreqVcf,
            inFileVcfVD = variantCalling.outFileVardictVcf,
            infileVcfVS = variantCalling.outFileVarscanVcf,
            inDirPCGRref = inDirPCGRref,
            refFa = refFa,
            sampleName = finalOutputName
    }
    
    output {
        File outFileVcf = vcfAnnotate.outFileVcf
        File outFileVcfIndex = vcfAnnotate.outFileVcfIndex
        File outFileMaf = vcfAnnotate.outFileMaf
        File outFileFlexdbHtml = vcfAnnotate.outFileFlexdbHtml
        File outFileHtml = vcfAnnotate.outFileHtml
    }
}

# TASK DEFINITIONS