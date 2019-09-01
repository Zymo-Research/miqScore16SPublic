from .... import generics
from ....utilities import validations

class GroupInfo(object):

    def __init__(self, seqType:str):
        self.seqType = seqType
        self.members = []

    def append(self, member:str):
        self.members.append(member)

    def __iter__(self):
        for member in self.members:
            yield member

    def __str__(self):
        return "%s: %s" %(self.seqType, self.members)


class MasterTable(generics.InputFile):

    def parseInputFile(self):
        import csv
        inputList = []
        masterTableFile = open(self.filePath, 'r')
        masterTable = csv.reader(masterTableFile)
        for line in masterTable:
            if line[0].startswith("#"):
                self.categoryTitles = line[6:]
                continue
            inputList.append(MasterTableLine(line))
        masterTableFile.close()
        if not self.validTableData(inputList):
            raise ValueError("Got invalid sample data in %s" %inputList)
        self.samples = inputList.copy()
        self.errorModelGroups = self.getErrorModelGroups()
        self.groupInfo = self.makeGroupInfoTable()
        self.lookupTable = self.makeLookupTable()
        self.categoryCount, self.categorySpace = self.getGroupInfo()

    def makeLookupTable(self):
        lookupTable = {}
        for index, sample in enumerate(self.samples):
            lookupName = "%s-%s" %(sample.sampleID, sample.sampleLabel)
            lookupTable[lookupName] = index
        return lookupTable

    def makeGroupInfoTable(self):
        groupInfo = {}
        for sample in self.samples:
            if not sample.groupID in groupInfo:
                groupInfo[sample.groupID] = GroupInfo(sample.seqType)
            groupInfo[sample.groupID].append(sample.sampleLabel)
        return groupInfo

    def getErrorModelGroups(self):
        errorModelGroups = set()
        for sample in self.samples:
            errorModelGroups.add(sample.errorModelGroup)
        return errorModelGroups

    def validTableData(self, table):
        uniqueIdentifiers = [line.identifier for line in table]
        testSet = set()
        for identifier in uniqueIdentifiers:
            if not identifier in testSet:
                testSet.add(identifier)
            else:
                raise ValueError("Found a duplicate Group ID and Unique Label for %s" %identifier)
        return True

    def getGroupInfo(self):
        categorySpace = []
        categoryCount = len(self.samples[0].categories)
        for i in range(categoryCount):
            categorySpace.append(set())
        for sample in self.samples:
            if not len(sample.categories) == categoryCount:
                raise ValueError("Found a line with a different number of categories than others. Please check formatting of the input csv.\nCSV: %s\nLine: %s" %(self.filePath, sample))
            for index, category in enumerate(sample.categories):
                categorySpace[index].add(category)
        return categoryCount, categorySpace

    def __getitem__(self, item):
        if type(item) == int:
            return self.samples[item]
        else:
            return self.samples[self.lookupTable[item]]

    def __iter__(self):
        for line in self.samples:
            yield line

    def __str__(self):
        output = ""
        for line in self.samples:
            output += str(line) + "\n"
        return output


class MasterTableLine(object):

    def __init__(self, lineArray:list):
        if not len(lineArray) >= 6:
            raise ValueError("Each line in the master table should have at least 6 elements.  A line (below) had %s elements.\n%s" %(len(lineArray), lineArray))
        number, projectID, runID, groupID, seqType, uniqueLabel = lineArray[:6]
        categories = lineArray[6:]
        self.setNumber(number)
        self.setProjectID(projectID)
        self.setRunID(runID)
        self.setGroupID(groupID)
        self.setSeqType(seqType)
        self.setUniqueLabel(uniqueLabel)
        self.setCategories(categories)
        self.identifier = (self.groupID.lower(), self.uniqueLabel.lower())
        self.sampleID = "%s_%s" %(projectID, number)
        self.groupID = "%s_%s" %(groupID, seqType)
        self.read1 = "%s_R1.fastq.gz" %self.sampleID
        self.read2 = "%s_R2.fastq.gz" %self.sampleID
        self.errorModelGroup = self.getErrorModelGroup()
        self.createdFiles = {}
        self.baseName, self.read1Base, self.read2Base = self.getBaseNames()
        if not self.uniqueLabel:
            self.sampleLabel = self.sampleID
        else:
            self.sampleLabel = self.uniqueLabel

    def getBaseNames(self):
        read1Base = self.read1.split(".")[0]
        read2Base = self.read2.split(".")[0]
        base = "_".join(read1Base.split("_")[:2])
        return base, read1Base, read2Base

    def setNumber(self, sampleName):
        self.sampleName = validations.naming.alphaNumericString(sampleName.strip(), replacement="")

    def setProjectID(self, projectID:str):
        self.projectID = validations.naming.alphaNumericString(projectID.strip(), replacement="")

    def setRunID(self, runID:str):
        import datetime
        runID = runID.strip()
        self.runID = validations.naming.alphaNumericString(runID, replacement="")
        rawDate = runID[-6:]
        if not rawDate.isdigit():
            self.runDate = None
        day = int(rawDate[0:2])
        month = int(rawDate[2:4])
        year = int(rawDate[4:6])
        self.runDate = datetime.date(year, month, day)

    def setGroupID(self, groupID:str):
        self.groupID = validations.naming.alphaNumericString(groupID, replacement="")

    def setSeqType(self, seqType:str):
        self.seqType = seqType.strip()

    def setUniqueLabel(self, uniqueLabel:str):
        self.uniqueLabel = validations.naming.alphaNumericString(uniqueLabel.strip(), replacement="")

    def setCategories(self, categories:list):
        self.categories = []
        for category in categories:
            if category:
                self.categories.append(validations.naming.alphaNumericString(category.strip()))
            else:
                self.categories.append(None)

    def getErrorModelGroup(self):
        return "%s_%s" %(self.runID, self.seqType)

    def addCreatedFile(self, path:str, identifier:str, direction:int = 1):
        try:
            hashval = hash(identifier)
        except TypeError:
            raise ValueError("Identifier value must be hashable. %s of type %s was given, which is not." %(identifier, type(identifier)))
        if not direction in [1, 2]:
            print("WARNING: Direction was given as something other than 1 or 2.")
        if not identifier in self.createdFiles:
            self.createdFiles[identifier] = {}
        self.createdFiles[identifier][direction] = path

    def __str__(self):
        identifier = "%s:%s" %(self.sampleID, self.seqType)
        return identifier