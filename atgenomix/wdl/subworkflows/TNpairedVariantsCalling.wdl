version 1.0

import "BamSomaticsniperCallingProcess.wdl" as bamSomaticsniperProcess
import "LoFreqSomaticCallingProcess.wdl" as lofreqProcess
import "MuseCallingProcess.wdl" as museProcess
import "Mutect2CallingProcess.wdl" as mutect2Process
import "VardictPairedCallingProcess.wdl" as vardictProcess
import "VarscanSomaticCallingProcess.wdl" as varscanProcess

# WORKFLOW DEFINITION

# Call variants using multiple caller with tumor normal paired mode
workflow TNpairedVariantsCalling {
    input {
        File inFileTumorBam
        File inFileTumorBamIndex
        File inFileNormalBam
        File inFileNormalBamIndex
        File inFileIntervalBed
        File inFileGermlineResource
        File inFileGermlineResourceIndex
        File? inFilePON
        File? inFilePONindex
        File refFa
        File refFai
        File refDict
        Float vardictMinimumAF
        String tumorSampleName
        String normalSampleName
        String sampleName
    }
    
    call bamSomaticsniperProcess.BamSomaticsniperCallingProcess as bamsomaticsniper {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileNormalBam = inFileNormalBam,
            refFa = refFa,
            refFai = refFai,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call lofreqProcess.LoFreqSomaticCallingProcess as lofreq {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            sampleName = sampleName
    }

    call museProcess.MuseCallingProcess as muse {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            refFa = refFa,
            refFai = refFai,
            sampleName = sampleName
    }

    call mutect2Process.Mutect2CallingProcess as mutect2 {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileGermlineResource = inFileGermlineResource,
            inFileGermlineResourceIndex = inFileGermlineResourceIndex,
            inFilePON = inFilePON,
            inFilePONindex = inFilePONindex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            refDict = refDict,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }

    call vardictProcess.VardictPairedCallingProcess as vardict {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileTumorBamIndex = inFileTumorBamIndex,
            inFileNormalBam = inFileNormalBam,
            inFileNormalBamIndex = inFileNormalBamIndex,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            minimumAF = vardictMinimumAF,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }
 
    call varscanProcess.VarscanSomaticCallingProcess as varscan {
        input:
            inFileTumorBam = inFileTumorBam,
            inFileNormalBam = inFileNormalBam,
            inFileIntervalBed = inFileIntervalBed,
            refFa = refFa,
            refFai = refFai,
            tumorSampleName = tumorSampleName,
            normalSampleName = normalSampleName,
            sampleName = sampleName
    }
 
    output {
        File outFileBamsomaticsniperPyVcfGz = bamsomaticsniper.outFilePythonFilterVcfGz
        File outFileBamsomaticsniperPyVcfIndex = bamsomaticsniper.outFilePythonFilterVcfIndex
        File outFileLofreqPyVcfGz = lofreq.outFilePythonFilterVcfGz
        File outFileLofreqPyVcfIndex = lofreq.outFilePythonFilterVcfIndex
        File outFileMusePyVcfGz = muse.outFilePythonFilterVcfGz
        File outFileMusePyVcfIndex = muse.outFilePythonFilterVcfIndex
        File outFileMutect2PyVcfGz = mutect2.outFilePythonFilterVcfGz
        File outFileMutect2PyVcfIndex = mutect2.outFilePythonFilterVcfIndex
        File outFileVardictPyVcfGz = vardict.outFilePythonFilterVcfGz
        File outFileVardictPyVcfIndex = vardict.outFilePythonFilterVcfIndex
        File outFileVarscanPyVcfGz = varscan.outFilePythonFilterVcfGz
        File outFileVarscanPyVcfIndex = varscan.outFilePythonFilterVcfIndex
        File outFileBamsomaticsniperVcfGz = bamsomaticsniper.outFileVcfGz
        File outFileBamsomaticsniperVcfIndex = bamsomaticsniper.outFileVcfIndex
        File outFileLofreqVcfGz = lofreq.outFileVcfGz
        File outFileLofreqVcfIndex = lofreq.outFileVcfIndex
        File outFileMuseVcfGz = muse.outFileVcfGz
        File outFileMuseVcfIndex = muse.outFileVcfIndex
        File outFileMutect2VcfGz = mutect2.outFileVcfGz
        File outFileMutect2VcfIndex = mutect2.outFileVcfIndex
        File outFileVardictVcfGz = vardict.outFileVcfGz
        File outFileVardictVcfIndex = vardict.outFileVcfIndex
        File outFileVarscanVcfGz = varscan.outFileVcfGz
        File outFileVarscanVcfIndex = varscan.outFileVcfIndex
    }
}