from .... import generics
from ....utilities import validations

parameterKeys = ("fwdPrimerLength",
                 "revPrimerLength",
                 "fwdReadLengthWithPrimer",
                 "revReadLengthWithPrimer",
                 "minSeqSize",
                 "refSeq",
                 "refTaxa",
                 "refTree",
                 "refAlign",
                 "rdpDatabase")

class PipelineParameters(generics.InputFile):

    def parseInputFile(self):
        if self.filePath.lower().endswith(".json"):
            self.loadJSON()
        elif self.filePath.lower().endswith(".csv"):
            self.loadCSV()
        else:
            file = open(self.filePath, 'r')
            firstCharacter = file.read(1)
            file.close()
            if firstCharacter == "{":
                self.loadJSON()
            else:
                self.loadCSV()
        self.lookupTable = {}
        self.makeLookupTable()

    def loadJSON(self):
        import json
        self.parameterTable = {}
        file = open(self.filePath, 'r')
        self.rawParameters = json.load(file)
        file.close()
        for pipeline in self.rawParameters:
            for key in parameterKeys:
                if not key in self.parameterTable[pipeline]:
                    raise KeyError("Unable to find the %s parameter for the %s pipeline." %(key, pipeline))
        for pipelineName in self.rawParameters:
            pipeline = self.rawParameters[pipelineName]
            self.parameterTable[pipelineName] = PipelineParameterSet(pipelineName, pipeline["fwdPrimerLength"], pipeline["revPrimerLength"], pipeline["fwdReadlengthWithPrimer"], pipeline["revReadLengthWithPrimer"], pipeline["minSeqSize"], pipeline["refSeq"], pipeline["refTaxa"], pipeline["refTree"], pipeline["refAlign"], pipeline["rdpDatabase"])

    def loadCSV(self):
        import csv
        self.parameterTable = {}
        file = open(self.filePath, 'r')
        csvHandle = csv.reader(file)
        for line in csvHandle:
            if line[0].startswith("#"):
                continue
            else:
                if not len(line) == 11:
                    raise ValueError("Pipeline parameter dataline length not equal to 11. %s items found.\n%s" %(len(line), line))
                seqType, fwdPrimerLength, revPrimerLength, fwdReadLengthWithPrimer, revReadLengthWithPrimer, minSeqSize, refSeq, refTaxa, refTree, refAlign, rdpDatabase = line
                fwdPrimerLength = int(fwdPrimerLength)
                revPrimerLength = int(revPrimerLength)
                revReadLengthWithPrimer = int(revReadLengthWithPrimer)
                fwdReadLengthWithPrimer = int(fwdReadLengthWithPrimer)
                minSeqSize = int(minSeqSize)
                if seqType in self.parameterTable:
                    raise KeyError("Duplicate sequence type found: %s" %seqType)
                parameterSet = PipelineParameterSet(seqType, fwdPrimerLength, revPrimerLength, fwdReadLengthWithPrimer, revReadLengthWithPrimer, minSeqSize, refSeq, refTaxa, refTree, refAlign, rdpDatabase)
                self.parameterTable[seqType] = parameterSet
        file.close()

    def saveJSON(self, outputFilePath:str):
        import json
        outputData = {}
        for pipelineName in self.parameterTable:
            outputData[pipelineName] = self.parameterTable[pipelineName].makeDict()
        outputFile = open(outputFilePath, 'w')
        json.dump(outputData, outputFile, indent = 4)
        outputFile.close()

    def makeLookupTable(self):
        for key in self.parameterTable:
            if key in self.lookupTable:
                raise KeyError("Error: Duplicate pipeline keys found: %s" %key)
            if key.lower() in self.lookupTable:
                raise KeyError("Error: Casing collision found for key %s" %key)
            self.lookupTable[key] = key
            self.lookupTable[key.lower()] = key

    def __getitem__(self, item):
        if not item.lower() in self.lookupTable:
            raise KeyError("Unable to find a pipeline parameter set called " %item)
        return self.parameterTable[self.lookupTable[item.lower()]]

    def __str__(self):
        return str(self.parameterTable)

    def __contains__(self, item):
        return item in self.lookupTable


class PipelineParameterSet(object):

    def __init__(self, name:str, fwdPrimerLength:int, revPrimerLength:int, fwdReadLengthWithPrimer:int, revReadLengthWithPrimer:int, minSeqSize:int, refSeq:str, refTaxa:str, refTree:str, refAlign:str, rdpDatabase:str):
        self.name = name
        self.fwdPrimerLength = fwdPrimerLength
        self.revPrimerLength = revPrimerLength
        self.fwdReadLengthWithPrimer = fwdReadLengthWithPrimer
        self.revReadLengthWithPrimer = revReadLengthWithPrimer
        self.minSeqSize = minSeqSize
        self.refSeq = refSeq
        self.refTaxa = refTaxa
        self.refTree = refTree
        self.refAlign = refAlign
        self.rdpDatabase = rdpDatabase
        self.performValidations()

    def makeDict(self):
        dictionary = {"fwdPrimerLength": self.fwdPrimerLength,
                      "revPrimerLength": self.revPrimerLength,
                     "fwdReadLengthWithPrimer": self.fwdReadLengthWithPrimer,
                     "revReadLengthWithPrimer": self.revReadLengthWithPrimer,
                     "minSeqSize": self.minSeqSize,
                     "refSeq": self.refSeq,
                     "refTaxa": self.refTaxa,
                     "refTree": self.refTree,
                     "refAlign": self.refAlign,
                     "rdpDatabase": self.rdpDatabase}
        return dictionary

    def performValidations(self):
        assert validations.numerical.isPositiveInteger(self.fwdPrimerLength), "Error. Forward primer length must be a positive integer. %s was given." %self.fwdPrimerLength
        self.fwdPrimerLength = int(self.fwdPrimerLength)
        assert validations.numerical.isPositiveInteger(self.revPrimerLength), "Error. Reverse primer length must be a positive integer. %s was given." %self.revPrimerLength
        self.revPrimerLength = int(self.revPrimerLength)
        assert self.fwdReadLengthWithPrimer > self.fwdPrimerLength, "Error. Forward read length with primer should always be longer than the primer itself. Primer length: %s, forward read length with primer: %s" %(self.fwdPrimerLength, self.fwdReadLengthWithPrimer)
        self.fwdReadLengthWithPrimer = int(self.fwdReadLengthWithPrimer)
        assert self.revReadLengthWithPrimer > self.revPrimerLength, "Error. Reverse read length with primer should always be longer than the primer itself. Primer length: %s, reverse read length with primer: %s" %(self.revPrimerLength, self.revReadLengthWithPrimer)
        self.revReadLengthWithPrimer = int(self.revReadLengthWithPrimer)
        assert validations.numerical.isPositiveInteger(self.minSeqSize), "Error. Minimum sequence size must be a positive integer. %s was given." %self.minSeqSize
        self.minSeqSize = int(self.minSeqSize)
        # KILLING THESE VALIDATIONS, THIS NEEDS TO HAPPEN OUTSIDE OF HERE, AS THE FOLDER WILL BE UNKNOWN BY THIS FUNCTION
        # if self.refSeq:
        #     refSeq = validations.system.fileExistsAndAbsolutePath(self.refSeq)
        #     if not refSeq:
        #         pass #raise FileNotFoundError("Unable to find reference sequence file at %s" %self.refSeq)
        #     self.refSeq = refSeq
        # else:
        #     self.refSeq = None
        # if self.refTaxa:
        #     refTaxa = validations.system.fileExistsAndAbsolutePath(self.refTaxa)
        #     if not refTaxa:
        #         pass #raise FileNotFoundError("Unable to find reference taxa file at %s" %self.refTaxa)
        #     self.refTaxa = refTaxa
        # else:
        #     self.refTaxa = None
        # if self.refTree:
        #     refTree = validations.system.fileExistsAndAbsolutePath(self.refTree)
        #     if not refTree:
        #         pass #raise FileNotFoundError("Unable to find reference tree file at %s" %self.refTree)
        #     self.refTree = refTree
        # else:
        #     self.refTree = None
        # if self.refAlign:
        #     refAlign = validations.system.fileExistsAndAbsolutePath(self.refAlign)
        #     if not refAlign:
        #         pass #raise FileNotFoundError("Unable to find reference alignment file at %s" %self.refAlign)
        #     self.refAlign = refAlign
        # else:
        #     self.refAlign = None
        # if self.rdpDatabase:
        #     rdpDatabase = validations.system.fileExistsAndAbsolutePath(self.rdpDatabase)
        #     if not rdpDatabase:
        #         pass #raise FileNotFoundError("Unable to find RDP database file at %s" %self.rdpDatabase)
        #     self.rdpDatabase = rdpDatabase
        # else:
        #     self.rdpDatabase = None

    def __str__(self):
        return str(self.makeDict())

