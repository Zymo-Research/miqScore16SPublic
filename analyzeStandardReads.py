import os
import logging
import miqScore16SPublicSupport
import defaults.standard as default
import Figaro
import miqScoreNGSReadCountPublic


def getApplicationParameters():
    parameters = miqScore16SPublicSupport.parameters.environmentParameterParser.EnvParameters()
    parameters.addParameter("sampleName", str, required=True, externalValidation=True)
    parameters.addParameter("maxReadCount", int, default=default.maxReadCount, lowerBound=10000, upperBound=20000000)
    parameters.addParameter("forwardReads", str, default = default.forwardReads, expectedFile=True)
    parameters.addParameter("reverseReads", str, default=default.reverseReads, expectedFile=True)
    parameters.addParameter("forwardPrimerLength", int, lowerBound=15, upperBound=40, required=True)
    parameters.addParameter("reversePrimerLength", int, lowerBound=15, upperBound=40, required=True)
    parameters.addParameter("sequenceFolder", str, default.sequenceFolder, expectedDirectory=True)
    parameters.addParameter("rScriptFolder", str, default=default.rScriptFolder, expectedDirectory=True)
    parameters.addParameter("r1MaxEE", int, default=default.r1MaxEE, lowerBound=0)
    parameters.addParameter("r2MaxEE", int, default=default.r2MaxEE, lowerBound=0)
    parameters.addParameter("truncQ", int, default=default.truncQ, lowerBound=0)
    parameters.addParameter("maxMismatch", int, default=default.maxMismatch, lowerBound=0)
    parameters.addParameter("minOverlap", int, default=default.minOverlap, lowerBound=0)
    parameters.addParameter("noCleanup", bool, default=False)
    parameters.addParameter("outputFolder", str, default=default.outputFolder, createdDirectory=True)
    parameters.addParameter("databaseFile", str, default=default.databaseFile, expectedFile=True)
    parameters.addParameter("trimParameterDownsample", int, default=default.trimParameterDownsample, lowerBound=-1)
    parameters.addParameter("ampliconLength", int, lowerBound = 100, upperBound = 300000, required = True)
    parameters.addParameter("noAutoTrimParameters", bool, default = False)
    parameters.addParameter("forwardReadTrim", int, default = 0, lowerBound = 0)
    parameters.addParameter("reverseReadTrim", int, default = 0, lowerBound = 0)
    parameters.addParameter("fileNamingStandard", str, default="zymo", externalValidation=True)
    parameters.addParameter("trimParameterPercentile", int, default = 83, upperBound=90, lowerBound=10)
    parameters.addParameter("trimParameterPickle", str, default="", externalValidation=True)
    parameters.addParameter("dada2OutputFiles", str, default="", externalValidation=True)
    parameters.addParameter("debug", bool, default=False)
    test = dict(os.environ)
    requiredCombinedLength = parameters.ampliconLength.value + parameters.minOverlap.value
    parameters.sideLoadParameter("minCombinedReadLength", requiredCombinedLength)
    if not parameters.fileNamingStandard.value.lower() in Figaro.figaroSupport.fileNamingStandards.aliasList.keys():
        raise ValueError("%s is not a valid naming standard alias" %parameters.fileNamingStandard.value)
    if not validSampleName(parameters.sampleName.value):
        logger.error("Invalid sample name given: %s" %parameters.sampleName.value)
        raise ValueError("Invalid sample name given: %s" %parameters.sampleName.value)
    parameters.checkCreatedFileStructures()
    if parameters.dada2OutputFiles.value or parameters.trimParameterPickle.value:
        if not parameters.debug.value:
            raise ValueError("Debugging options used without debug being set to true. Please check into this.")
    return parameters


def validateFastqPair(forwardPath:str, reversePath:str):
    readCount = miqScore16SPublicSupport.formatReaders.fastq.fastqHandler.validFastqPair(forwardPath, reversePath)
    if not readCount:
        if readCount == False:
            errorMessage = "Fastq file failed validation checks"
        elif readCount == 0:
            errorMessage = "Fastq files appear to be empty of reads"
        else:
            raise RuntimeError("Fastq validation should only return either False or a non-negative integer. This code should be unreachable and this is a bug.")
        logger.error(errorMessage)
        raise miqScore16SPublicSupport.formatReaders.fastq.fastqHandler.FastqFormatError(errorMessage)
    if readCount > parameters.maxReadCount.value:
        raise RuntimeError("Max fastq read count exceeded for this sample. Max: %s. Read count: %s" %(parameters.maxReadCount.value, readCount))
    return True


def validSampleName(name:str):
    validNameCharacters = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890.-_ "
    invalidStartingCharacters = ".-_ "
    if name[0] in invalidStartingCharacters:
        return False
    for character in name:
        if not character in validNameCharacters:
            return False
    return True


def getLoggingParameters():
    loggingParameters = miqScore16SPublicSupport.parameters.environmentParameterParser.EnvParameters()
    loggingParameters.addParameter("logFile", str, default=default.logFile, createdFile=True)
    loggingParameters.addParameter("logLevel", str, default=default.loggingLevel, logLevel=True)
    loggingParameters.addParameter("streamOff", bool, default=False)
    loggingParameters.addParameter("streamLoglevel", str, default=default.loggingLevel, logLevel=True)
    loggingParameters.addParameter("fileLogLevel", str, default=default.loggingLevel, logLevel=True)
    logFilePath = os.path.split(loggingParameters.logFile.value)[0]
    if not os.path.isdir(logFilePath):
        os.makedirs(logFilePath)
    loggingParameters.checkCreatedFileStructures()
    return loggingParameters


def loadDefaultPackage():
    defaultParameters = miqScore16SPublicSupport.parameters.environmentParameterParser.EnvParameters()
    defaultParameters.addParameter("defaultPackageName", str, default="standard", externalValidation=True)
    return miqScore16SPublicSupport.parameters.defaultParser.loadDefaultModule(defaultParameters.defaultPackageName.value)


def setLogging():
    loggingParameters = getLoggingParameters()
    formatter = logging.Formatter(loggingFormat)
    logStreamHandle = logging.StreamHandler()
    logStreamHandle.setFormatter(formatter)
    if not loggingParameters.streamLogLevel.usingDefaultValue:
        logStreamHandle.setLevel(loggingParameters.streamLogLevel.value)
    else:
        logStreamHandle.setLevel(loggingParameters.logLevel.value)
    logFileHandle = logging.FileHandler(loggingParameters.logFile.value)
    logFileHandle.setFormatter(formatter)
    if not loggingParameters.fileLogLevel.usingDefaultValue:
        logFileHandle.setLevel(loggingParameters.fileLogLevel.value)
    else:
        logFileHandle.setLevel(loggingParameters.logLevel.value)
    logger.addHandler(logFileHandle)
    if not loggingParameters.streamOff:
        logger.addHandler(logStreamHandle)


class Dada2OutputFiles(object):

    __slots__ = ["rawReads", "trimmedReads", "errorModels", "amplicons", "chimeraFreeAmplicons", "taxa", "rdsFile"]

    def __init__(self):
        self.rawReads = None
        self.trimmedReads = None
        self.errorModels = None
        self.amplicons = None
        self.chimeraFreeAmplicons = None
        self.taxa = None
        self.rdsFile = None


def getReadLengthsFromFastq(path:str):
    return miqScore16SPublicSupport.formatReaders.fastq.fastqHandler.estimateReadLength(path)

 
def buildRscriptArguments(arguments:dict):
    rscriptArgumentsList = []
    for flag in arguments:
        rscriptArgumentsList.append("-%s %s" % (flag, arguments[flag]))
    rscriptArguments = " ".join(rscriptArgumentsList)
    return rscriptArguments


def getTrimmingParametersWithFigaro(sequenceFolder):
    import Figaro
    if parameters.trimParameterPickle.value:
        import pickle
        resultTablePickle = open(parameters.trimParameterPickle.value, 'rb')
        resultTable = pickle.load(resultTablePickle)
        resultTablePickle.close()
    else:
        resultTable, forwardCurve, reverseCurve = Figaro.figaro.runAnalysis(sequenceFolder, parameters.minCombinedReadLength.value, parameters.forwardPrimerLength.value, parameters.reversePrimerLength.value, parameters.fileNamingStandard.value, parameters.trimParameterDownsample.value, parameters.trimParameterPercentile.value)
        # import pickle
        # file = open("/data/output/trimParameters.pkl", 'wb')
        # pickle.dump(resultTable, file)
        # file.close()
    return resultTable[0]


def dada2Trim(forwardReads:str, reverseReads:str, trimParameters:Figaro.figaroSupport.trimParameterPrediction.TrimParameterSet):
    r1TrimmedFile = os.path.join(parameters.outputFolder.value, parameters.sampleName.value + ".forward.trimmed.fastq.gz")
    r2TrimmedFile = os.path.join(parameters.outputFolder.value, parameters.sampleName.value + ".reverse.trimmed.fastq.gz")
    arguments = {
        "f" : forwardReads,
        "r" : reverseReads,
        "x" : r1TrimmedFile,
        "y" : r2TrimmedFile,
        "a" : parameters.forwardPrimerLength.value,
        "b" : parameters.reversePrimerLength.value,
        "c" : trimParameters.forwardTrimPosition,
        "d" : trimParameters.reverseTrimPosition,
        "q" : parameters.truncQ,
        "m" : trimParameters.forwardMaxExpectedError,
        "n" : trimParameters.reverseMaxExpectedError
    }
    rscriptCommand = "%s %s" %(default.rScriptExecutable, os.path.join(parameters.rScriptFolder.value, "dada2.trimreads.R"))
    rscriptArguments = buildRscriptArguments(arguments)
    command = "%s %s" %(rscriptCommand, rscriptArguments)
    logger.info("Running trim command: %s" %command)
    exitCode = os.system(command)
    if exitCode == 0:
        logger.info("Trim command returned exit status 0")
    else:
        logger.error("Trim command %s returned exit status %s" %(command, exitCode))
        raise RuntimeError("Error in trimming operation")
    return r1TrimmedFile, r2TrimmedFile


def dada2BuildErrorModelRunner(inputFileList:str, outputFileName:str):
    arguments = {
        "i" : inputFileList,
        "o" : outputFileName
    }
    rscriptCommand = "%s %s" %(default.rScriptExecutable, os.path.join(parameters.rScriptFolder.value, "dada2.builderrormodels.R"))
    rscriptArguments = buildRscriptArguments(arguments)
    command = "%s %s" %(rscriptCommand, rscriptArguments)
    logger.info("Running error model command: %s" %command)
    exitCode = os.system(command)
    if exitCode == 0:
        logger.info("Error model command returned exit status 0")
    else:
        logger.error("Error model command returned exit status %s" %exitCode)
        raise RuntimeError("Error in error modeling process.")
    return exitCode


def dada2BuildErrorModels(forwardTrimmedReads:str, reverseTrimmedReads:str):
    pe1FileList = os.path.join(parameters.outputFolder.value, "read1.fileList.txt")
    pe2FileList = os.path.join(parameters.outputFolder.value, "read2.fileList.txt")
    pe1ErrorModelFile = os.path.join(parameters.outputFolder.value, parameters.sampleName.value + "_1.Rda")
    pe2ErrorModelFile = os.path.join(parameters.outputFolder.value, parameters.sampleName.value + "_2.Rda")
    pe1List = [forwardTrimmedReads]
    pe2List = [reverseTrimmedReads]
    pe1File = open(pe1FileList, 'w')
    for path in pe1List:
        print(path, file=pe1File)
    pe1File.close()
    pe2File = open(pe2FileList, 'w')
    for path in pe2List:
        print(path, file=pe2File)
    pe2File.close()
    dada2BuildErrorModelRunner(pe1FileList, pe1ErrorModelFile)
    dada2BuildErrorModelRunner(pe2FileList, pe2ErrorModelFile)
    if not parameters.noCleanup:
        os.remove(pe1FileList)
        os.remove(pe2FileList)
    return pe1ErrorModelFile, pe2ErrorModelFile


def dada2GetAmplicons(forwardTrimmedReads:str, reverseTrimmedReads:str, forwardErrorModel:str, reverseErrorModel:str):
    outputFolder = parameters.outputFolder.value
    seqTableCSV = os.path.join(outputFolder, parameters.sampleName.value + ".SV.csv")
    outputTaxa = os.path.join(outputFolder, parameters.sampleName.value + ".SV.taxa.csv")
    outputTaxaChimeraFree = os.path.join(outputFolder, parameters.sampleName.value + ".SV.nochimera.csv")
    outputSequenceTable = os.path.join(outputFolder, parameters.sampleName.value + ".seqtable.rds")
    rdpDataBase = parameters.databaseFile.value
    arguments = {
        "f" : forwardTrimmedReads,
        "r" : reverseTrimmedReads,
        "x" : forwardErrorModel,
        "y" : reverseErrorModel,
        "o" : seqTableCSV,
        "t" : outputTaxa,
        "k" : outputTaxaChimeraFree,
        "s" : outputSequenceTable,
        "d" : rdpDataBase,
        "l" : parameters.minOverlap.value,
        "m" : parameters.maxMismatch.value
    }
    rscriptCommand = "%s %s" %(default.rScriptExecutable, os.path.join(parameters.rScriptFolder.value, "dada2.getamplicons.R"))
    rscriptArguments = buildRscriptArguments(arguments)
    command = "%s %s" %(rscriptCommand, rscriptArguments)
    logger.info("Running amplicon calling command: %s" %command)
    exitCode = os.system(command)
    if exitCode == 0:
        logger.info("Amplicon calling command returned exit status 0")
    else:
        logger.error("Amplicon calling command returned exit status %s" %exitCode)
        raise RuntimeError("Error in amplicon calling step")
    return seqTableCSV, outputTaxaChimeraFree, outputTaxa, outputSequenceTable


def addErrorModelsToSampleTreeAndMasterFile(errorModelTable:dict, samples:miqScore16SPublicSupport.projectData.microbiome.sixteenS.metadata.masterTable.MasterTable, sampleFileTree:dict):
    for sample in samples:
        errorModelFilePe1 = errorModelTable[sample.errorModelGroup][1]
        errorModelFilePe2 = errorModelTable[sample.errorModelGroup][2]
        sample.addCreatedFile(errorModelFilePe1, "errorModel", 1)
        sample.addCreatedFile(errorModelFilePe2, "errorModel", 2)
        sampleFileTree[sample.baseName]["errorModel"] = {}
        sampleFileTree[sample.baseName]["errorModel"][1] = os.path.split(errorModelFilePe1)[1]
        sampleFileTree[sample.baseName]["errorModel"][2] = os.path.split(errorModelFilePe2)[1]


def removeFileRedundantEntries(masterTable:miqScore16SPublicSupport.projectData.microbiome.sixteenS.metadata.masterTable.MasterTable):
    deduped = []
    alreadyEntered = set()
    for sample in masterTable:
        if sample.baseName in alreadyEntered:
            continue
        else:
            deduped.append(sample)
            alreadyEntered.add(sample.baseName)
    return deduped

def runDada2Functions(forwardReads:str, reverseReads:str):
    if parameters.dada2OutputFiles.value:
        import pickle
        if parameters.trimParameterPickle.value:
            readFolder = os.path.split(os.path.abspath(forwardReads))[0]
            trimParameters = getTrimmingParametersWithFigaro(readFolder)
        file = open(parameters.dada2OutputFiles.value, 'rb')
        outputFiles = pickle.load(file)
        outputFiles.rawReads = (forwardReads, reverseReads)
        file.close()
    else:
        import miqScore16SPublicSupport.projectData.microbiome.dada2Outputs
        outputFiles = Dada2OutputFiles()
        outputFiles.rawReads = (forwardReads, reverseReads)
        readFolder = os.path.split(os.path.abspath(forwardReads))[0]
        trimParameters = getTrimmingParametersWithFigaro(readFolder)
        outputFiles.trimmedReads = dada2Trim(forwardReads, reverseReads, trimParameters)
        outputFiles.errorModels = dada2BuildErrorModels(*outputFiles.trimmedReads)
        outputFiles.amplicons, outputFiles.chimeraFreeAmplicons, outputFiles.taxa, outputFiles.rdsFile = dada2GetAmplicons(*outputFiles.trimmedReads, *outputFiles.errorModels)
        # import pickle
        # file = open("/data/output/dada2OutputFiles.pkl", 'wb')
        # pickle.dump(outputFiles, file)
        # file.close()
    return outputFiles


def dada2GenusCountTable(chimeraFreeReadCounts:miqScore16SPublicSupport.projectData.microbiome.dada2Outputs.Dada2AmpliconCount, taxaMapping:miqScore16SPublicSupport.projectData.microbiome.dada2Outputs.Dada2KingdomGenusCalls):
    ampliconToGenusTable = {}
    for amplicon in taxaMapping.callTable:
        ampliconToGenusTable[amplicon] = taxaMapping.callTable[amplicon][5]
    genusReadCounts = {}
    for amplicon in chimeraFreeReadCounts.ampliconTable:
        if amplicon in ampliconToGenusTable:
            genus = ampliconToGenusTable[amplicon]
        else:
            import hashlib
            genus = hashlib.md5(amplicon.encode("utf-8")).hexdigest()
        if not genus in genusReadCounts:
            genusReadCounts[genus] = 0
        genusReadCounts[genus] += chimeraFreeReadCounts.ampliconTable[amplicon]
    return genusReadCounts

def getDada2Results(dada2OutputFiles:Dada2OutputFiles):
    from . import miqScore16SPublicSupport
    totalReadInput = miqScore16SPublicSupport.formatReaders.fastq.fastqHandler.countReads(dada2OutputFiles.rawReads[0])
    trimmedReadInput = miqScore16SPublicSupport.formatReaders.fastq.fastqHandler.countReads(dada2OutputFiles.trimmedReads[0])
    ampliconsWithChimera = miqScore16SPublicSupport.projectData.microbiome.dada2Outputs.Dada2AmpliconCount(dada2OutputFiles.amplicons)
    ampliconsWithoutChimera = miqScore16SPublicSupport.projectData.microbiome.dada2Outputs.Dada2AmpliconCount(dada2OutputFiles.chimeraFreeAmplicons)
    totalReads = totalReadInput
    filteredReads = totalReads - trimmedReadInput
    unmergedReads = trimmedReadInput - ampliconsWithChimera.readCount
    chimericReads = ampliconsWithChimera.readCount - ampliconsWithoutChimera.readCount
    readAssignments = miqScore16SPublicSupport.projectData.microbiome.dada2Outputs.Dada2KingdomGenusCalls(dada2OutputFiles.taxa)
    genusReadCounts = dada2GenusCountTable(ampliconsWithoutChimera, readAssignments)
    readFateTable = {"totalReads" : totalReads,
                     "filteredReads" : filteredReads,
                     "unmergedReads" : unmergedReads,
                     "chimericReads" : chimericReads,
                     "calledReads" : ampliconsWithoutChimera.readCount,
                     "genusCalls" : genusReadCounts}
    return readFateTable


readFatePrintNames = {"filteredReads": "Failed Quality Filter",
                      "unmergedReads": "Failed To Merge",
                      "chimericReads": "Chimeric",
                      "Reference": "Aligned To Reference"}


def analyzeStandardResult(dada2ResultTable:dict):
    cleanedTable = dada2ResultTable.copy()
    del cleanedTable["totalReads"]
    del cleanedTable["calledReads"]
    for genus in cleanedTable["genusCalls"]:
        cleanedTable[genus] = cleanedTable["genusCalls"][genus]
    del cleanedTable["genusCalls"]
    referenceDataFile = os.path.join(os.path.split(__file__)[0], "reference", "zrCommunityStandard.json")
    referenceData = miqScoreNGSReadCountPublic.referenceHandler.StandardReference(referenceDataFile)
    calculator = miqScoreNGSReadCountPublic.MiqScoreCalculator(referenceData, analysisMethod="16s", percentToleranceInStandard=15, floor=0)
    miqScoreResult = calculator.calculateMiq(cleanedTable, parameters.sampleName.value)
    miqScoreResult.makeReadFateChart(readFatePrintNames=readFatePrintNames)
    miqScoreResult.makeRadarPlots()
    goodMiqPath = os.path.join(os.path.split(__file__)[0], "reference", "goodMiq.json")
    badMiqPath = os.path.join(os.path.split(__file__)[0], "reference", "badMiq.json")
    goodComposition, badComposition = miqScoreNGSReadCountPublic.loadReferenceCompositionFromExampleMiq(goodMiqPath, badMiqPath)
    miqScoreResult.makeCompositionBarPlot(goodComposition, badComposition)
    return miqScoreResult


def saveResult(result:miqScoreNGSReadCountPublic.MiqScoreData):
    outputFilePath = os.path.join(parameters.outputFolder.value, "%s.json" %parameters.sampleName.value)
    print("Output results to %s" %outputFilePath)
    outputFile = open(outputFilePath, 'w')
    outputFile.write(result.jsonOutput())
    outputFile.close()
    return outputFilePath


def generateReport(result:miqScoreNGSReadCountPublic.MiqScoreData):
    referenceDataFile = os.path.join(os.path.split(__file__)[0], "reference", "zrCommunityStandard.json")
    referenceData = miqScoreNGSReadCountPublic.referenceHandler.StandardReference(referenceDataFile)
    templateFilePath = os.path.join(os.path.split(os.path.abspath(__file__))[0], "reference", "16SReportTemplate.html")
    templateFile = open(templateFilePath, 'r')
    template = templateFile.read()
    templateFile.close()
    goodMiqPath = os.path.join(os.path.split(__file__)[0], "reference", "goodMiq.json")
    badMiqPath = os.path.join(os.path.split(__file__)[0], "reference", "badMiq.json")
    goodMiq, badMiq = miqScoreNGSReadCountPublic.loadExampleData(goodMiqPath, badMiqPath, referenceData)
    replacementTable = miqScore16SPublicSupport.reporting.generateReplacementTable(result, goodMiq, badMiq, readFatePrintNames=readFatePrintNames)
    report = miqScoreNGSReadCountPublic.reportGeneration.generateReport(template, replacementTable)
    reportFilePath = os.path.join(parameters.outputFolder.value, "%s.html" % parameters.sampleName.value)
    print("Output report to %s" % reportFilePath)
    outputFile = open(reportFilePath, 'w')
    outputFile.write(report)
    outputFile.close()
    return reportFilePath


if __name__ == "__main__":
    default = loadDefaultPackage()
    loggingFormat = "%(levelname)s:%(name)s:%(message)s"
    logger = logging.getLogger(__name__)
    logger.setLevel(
        logging.DEBUG)  # Do not change this line unless you know exactly what you are doing any why you are doing it. This will mess up logging in a way that can be hard to trace back.
    setLogging()
    parameters = getApplicationParameters()
    logger.debug("Starting analysis")
    dada2Outputs = runDada2Functions(parameters.forwardReads.value, parameters.reverseReads.value)
    dada2Results = getDada2Results(dada2Outputs)
    standardAnalysisResults = analyzeStandardResult(dada2Results)
    saveResult(standardAnalysisResults)
    generateReport(standardAnalysisResults)
    exit(0)