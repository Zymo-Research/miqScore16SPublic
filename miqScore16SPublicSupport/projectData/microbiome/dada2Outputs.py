import os
import logging
logger = logging.getLogger(__name__)

class Dada2AmpliconCount(object):

    def __init__(self, path:str):
        if not os.path.isfile(path):
            logger.critical("Tried to load dada2 amplicon count table, but could not find file at %s" %path)
            raise FileNotFoundError("Unable to find file %s" %path)
        self.path = path
        self.ampliconTable = self.readTable(path)
        self.readCount = self.getReadCount()
        self.ampliconList = tuple(self.ampliconTable.keys())

    def readTable(self, path):
        import csv
        rawTable = []
        file = open(path, 'r')
        csvHandle = csv.reader(file)
        for line in csvHandle:
            rawTable.append(line[1:])
        file.close()
        zipped = zip(rawTable[0], rawTable[1])
        ampliconTable = {}
        for line in zipped:
            amplicon, count = line
            ampliconTable[amplicon] = int(count)
        return ampliconTable

    def getReadCount(self):
        readCount = 0
        for amplicon in self.ampliconTable:
            readCount += self.ampliconTable[amplicon]
        return readCount

    def __iter__(self):
        for amplicon in self.ampliconList:
            yield amplicon

    def __getitem__(self, item):
        return self.ampliconTable[item]

    def __str__(self):
        return "Dada2 amplicon table. %s amplicons. %s reads. From %s" %(len(self.ampliconList), self.readCount, self.path)


class Dada2GenusSpeciesCalls(object):

    def __init__(self, path:str):
        if not os.path.isfile(path):
            logger.critical("Tried to load dada2 amplicon count table, but could not find file at %s" %path)
            raise FileNotFoundError("Unable to find file %s" %path)
        self.path = path
        self.callTable = self.makeSequenceCallTable()
        self.ampliconList = list(self.callTable.keys())

    def makeSequenceCallTable(self):
        file = open(self.path, "r")
        callTable = {}
        import csv
        csvHandle = csv.reader(file)
        for line in csvHandle:
            sequence, genus, species = line
            if not sequence:
                continue
            callTable[sequence] = (genus, species)
        return callTable

    def __iter__(self):
        for amplicon in self.ampliconList:
            yield amplicon

    def __getitem__(self, item):
        return self.callTable[item]

    def __str__(self):
        return "Dada2 genus/species call table. %s amplicons. From %s" %(len(self.ampliconList), self.path)


class Dada2KingdomGenusCalls(object):

    def __init__(self, path:str):
        if not os.path.isfile(path):
            logger.critical("Tried to load dada2 amplicon count table, but could not find file at %s" %path)
            raise FileNotFoundError("Unable to find file %s" %path)
        self.path = path
        self.callTable = self.makeSequenceCallTable()
        self.ampliconList = list(self.callTable.keys())

    def makeSequenceCallTable(self):
        file = open(self.path, "r")
        callTable = {}
        import csv
        csvHandle = csv.reader(file)
        for line in csvHandle:
            sequence, kingdom, phylum, className, order, family, genus = line
            if not sequence:
                continue
            callTable[sequence] = (kingdom, phylum, className, order, family, genus)
        return callTable

    def __iter__(self):
        for amplicon in self.ampliconList:
            yield amplicon

    def __getitem__(self, item):
        return self.callTable[item]

    def __str__(self):
        return "Dada2 genus/species call table. %s amplicons. From %s" %(len(self.ampliconList), self.path)


def getHitCountByTaxa(ampliconCountTable:[str, Dada2AmpliconCount], taxaCallTable:[str, Dada2GenusSpeciesCalls, Dada2KingdomGenusCalls], genusSpeciesTable:bool = False):
    hitCountTable = {}
    if type(ampliconCountTable) == str:
        ampliconCountTable = Dada2AmpliconCount(ampliconCountTable)
    if type(taxaCallTable) == str:
        if not genusSpeciesTable:
            taxaCallTable = Dada2KingdomGenusCalls(taxaCallTable)
        else:
            taxaCallTable = Dada2GenusSpeciesCalls(taxaCallTable)
    for sequence in ampliconCountTable:
        if not sequence in taxaCallTable:
            continue
        taxa = taxaCallTable[sequence]
        if not taxa in hitCountTable:
            hitCountTable[taxa] = 0
        hitCountTable[taxa] += ampliconCountTable[sequence]
    return hitCountTable


if __name__ == "__main__":
    chimeras = Dada2AmpliconCount("in1055_1.SV.csv")
    chimeraFree = Dada2AmpliconCount("in1055_1.SV.nochimera.csv")
    taxa = Dada2KingdomGenusCalls("in1055_1.SV.taxa.csv")
    hitCounts = getHitCountByTaxa(chimeraFree.ampliconTable, taxa.callTable)
    print("something")