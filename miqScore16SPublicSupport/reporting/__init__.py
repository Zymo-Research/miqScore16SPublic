def sumDictionary(dictionary:dict):
    sum = 0
    for key in dictionary:
        if key:
            sum += dictionary[key]
    return sum


def generateReadFateChartBody(readFateTable:dict, readFatePrintNames:dict=None):
    printReadFateTable = {}
    if readFatePrintNames:
        for readFate in readFateTable:
            if readFate in readFatePrintNames:
                printReadFateTable[readFatePrintNames[readFate]] = readFateTable[readFate]
            else:
                printReadFateTable[readFate] = readFateTable[readFate]
    else:
        printReadFateTable = readFateTable.copy()
    outputTable = ""
    for fate in printReadFateTable:
        outputTable += '\
            <tr style="height: 21px;">\n\
            <td style="width: 50%%; height: 21px;">%s</td>\n\
            <td style="width: 50%%; height: 21px;">%s</td>\n\
            </tr>\
    ' %(fate, round(printReadFateTable[fate], 2))
    return outputTable


def generateAbsoluteReadFateCounts(miqScore):
    referenceCounts = sumDictionary(miqScore.referenceReadCounts)
    absoluteReadFates = miqScore.nonreferenceReadCounts.copy()
    absoluteReadFates["Reference"] = referenceCounts
    return absoluteReadFates


def generateReplacementTable(sampleMiq, goodExampleMiq, badExampleMiq, readFatePrintNames:dict=None):
    readFateTable = generateReadFateChartBody(generateAbsoluteReadFateCounts(sampleMiq), readFatePrintNames)
    replacementTable = {"SAMPLENAME": sampleMiq.sampleID,
                        "MIQSCORE": str(round(sampleMiq.miqScore)),
                        "READFATETABLE": readFateTable,
                        "READFATECHART": sampleMiq.plots["readFates"],
                        "COMPOSITIONBARPLOT": sampleMiq.plots["compositionPlot"],
                        "GOODRADARPLOT": goodExampleMiq.plots["radarPlots"]["Lysis Difficulty"],
                        "SAMPLERADARPLOT": sampleMiq.plots["radarPlots"]["Lysis Difficulty"],
                        "BADRADARPLOT": badExampleMiq.plots["radarPlots"]["Lysis Difficulty"]}
    return replacementTable