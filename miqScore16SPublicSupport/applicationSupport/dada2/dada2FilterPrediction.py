from ... import formatReaders


def predictReadsKeptAfterExpectedErrorFilteringPaired(forwardReads:str, reverseReads:str, trimPositions:tuple, maxExpectedError:tuple, percentage:bool=False, forwardSkip:int = 0, reverseSkip:int = 0):
    '''
    Simulates the dada2 trimming process to see how many reads pass without as much work
    :param forwardReads: Path to the forward reads
    :param reverseReads:  Path to the reverse reads
    :param trimPositions: Either a tuple with a forward trim position and reverse trim position (both integers) -OR- a tuple of these tuples to try multiples
    :param maxExpectedError: Either a tuple with a forward max expected error and reverse max expected error (both integers) -OR- a tuple of these tuples to try multiples
    :param percentage: If true, results will be given in percentages instead of read counts
    :param forwardSkip Use if the first n bases of the sequence are being left off the analysis (removing bases from before the first potential trim point)
    :param reverseSkip Just like the forward skip, but for the other direction. Useful to keep matrices smaller.
    :return: If a single set of values are given for both expected error and trim positions, a single number will be returned.  Otherwise, a dictionary will be where the first key will be trimming and second will be expected error.
    '''
    def validPairParameters(parameter):
        if not type(parameter) == tuple:
            return False
        if not type(parameter[0]) == tuple:
            if not len(parameter) == 2:
                return False
            for value in parameter:
                if not type(value) == int:
                    return False
        else:
            for pair in parameter:
                if not len(pair) == 2:
                    return False
                for value in pair:
                    if not type(value) == int:
                        return False
        return True
    if not validPairParameters(trimPositions):
        raise TypeError("Trim positions should either be a single list of two integers with a forward and reverse trim position to test or a list of these value lists")
    if not validPairParameters(maxExpectedError):
        raise TypeError("Max expected error should be either a single list of two integers with a forward and reverse max expected errors to test or a list of these value lists.")
    if type(maxExpectedError[0]) == int:
        maxExpectedErrorTuple = (maxExpectedError)
    else:
        maxExpectedErrorTuple = maxExpectedError
    if type(trimPositions[0]) == int:
        trimPositionsTuple = (trimPositions)
    else:
        trimPositionsTuple = trimPositions
    maxExpectedErrorSet = set(maxExpectedErrorTuple)
    trimPositions = set(trimPositionsTuple)
    forwardMatrix, reverseMatrix = formatReaders.fastq.fastqAnalysis.buildExpectedErrorMatrixPaired(forwardReads, reverseReads)
    if not len(forwardMatrix) == len(reverseMatrix):
        raise formatReaders.fastq.fastqHandler.FastqValidationError(
            "Forward and reverse read fastq files %s and %s appear to have different numbers of reads: %s and %s (respectively)." % (
            forwardReads, reverseReads, len(forwardMatrix), len(reverseMatrix)))
    forwardMatrix = forwardMatrix.transpose()
    reverseMatrix = reverseMatrix.transpose()  #making read positions rows for faster searching by position across all reads
    resultTable = {}
    for trimPositionSet in trimPositions:
        trimPositionSet = tuple(trimPositionSet)
        resultTable[trimPositionSet] = {}
        forwardTrimPosition, reverseTrimPosition = trimPositionSet
        forwardExpectedErrors = forwardMatrix[forwardTrimPosition - forwardSkip]
        reverseExpectedErrors = reverseMatrix[reverseTrimPosition - reverseSkip]
        observedForwardAndReverseExpectedErrorsAtTrimPosition = [item for item in zip(forwardExpectedErrors, reverseExpectedErrors)]
        for maxExpectedErrorValuePair in maxExpectedErrorSet:
            maxExpectedErrorValuePair = tuple(maxExpectedErrorValuePair)
            forwardMaxExpectedError, reverseMaxExpectedError = maxExpectedErrorValuePair
            totalReads = 0
            keptReads = 0
            for observedForwardExpectedError, observedReverseExpectedError in observedForwardAndReverseExpectedErrorsAtTrimPosition:
                if observedForwardExpectedError <= forwardMaxExpectedError and observedReverseExpectedError <= reverseMaxExpectedError:
                    keptReads += 1
                totalReads += 1
            if percentage:
                storeValue = keptReads/totalReads
            else:
                storeValue = keptReads
            resultTable[trimPositionSet][maxExpectedErrorValuePair] = storeValue
    if len(resultTable) == 1:
        trimKey = list(resultTable.keys())[0]
        if len(resultTable[trimKey]) == 1:
            errorKey = list(resultTable.keys())[0]
            return resultTable[trimKey][errorKey]
    else:
        return resultTable


def predictReadsKeptAfterExpectedErrorFilteringSingle(path:str, trimPositions:[int, tuple], maxExpectedError:[int, tuple], percentage:bool=False):
    '''
    Simulates the dada2 trimming process to see how many reads pass without as much work
    :param path: Path to the fastq file
    :param trimPositions: Either a single positive integer value or a tuple of positive integer values to use as the trim position
    :param maxExpectedError: Either a single positive integer value or a tuple of positive integer values to use as the max expected error
    :param percentage: If true, results will be given in percentages instead of read counts
    :return: If a single value is given for both expected error and trim positions, a single number will be returned.  Otherwise, a dictionary will be where the first key will be trimming and second will be expected error.
    '''
    import numpy
    def validSingleReadParameters(parameter):
        if not type(parameter) in [int, tuple] :
            return False
        if type(parameter) == int:
            if parameter <= 0:
                return False
        else:
            for value in parameter:
                if not type(value) == int:
                    return False
                if value <= 0:
                    return False
        return True
    if not validSingleReadParameters(trimPositions):
        raise TypeError("Trim positions should be either a single positive integer value or a tuple of positive integer values")
    if not validSingleReadParameters(maxExpectedError):
        raise TypeError("Max expected error should be either a single positive integer value or a tuple of positive integer values")
    if type(maxExpectedError) == int:
        maxExpectedErrorTuple = (maxExpectedError)
    else:
        maxExpectedErrorTuple = maxExpectedError
    if type(trimPositions) == int:
        trimPositionsList = (trimPositions)
    else:
        trimPositionsList = trimPositions
    maxExpectedErrorSet = set(maxExpectedErrorTuple)
    trimPositionSet = set(trimPositionsList)
    expectedErrorMatrix = formatReaders.fastq.fastqAnalysis.buildExpectedErrorMatrix(path)
    expectedErrorMatrix = expectedErrorMatrix.transpose() #making read positions rows for faster searching by position across all reads. If you are curious why, talk to Mike W.; this is a good CS learning opportunity for you.
    resultTable = {}
    for trimPosition in trimPositionSet:
        resultTable[trimPosition] = {}
        observedExpectedErrorsAtTrimPosition = expectedErrorMatrix[trimPosition]
        for maxExpectedErrorValue in maxExpectedErrorSet:
            totalReads = 0
            keptReads = 0
            for observedExpectedError in numpy.nditer(observedExpectedErrorsAtTrimPosition):
                if observedExpectedError <= maxExpectedErrorValue:
                    keptReads += 1
                totalReads += 1
            if percentage:
                storeValue = keptReads/totalReads
            else:
                storeValue = keptReads
            resultTable[trimPosition][maxExpectedErrorValue] = storeValue
    if len(resultTable) == 1:
        trimKey = list(resultTable.keys())[0]
        if len(resultTable[trimKey]) == 1:
            errorKey = list(resultTable.keys())[0]
            return resultTable[trimKey][errorKey]
    else:
        return resultTable


def plotFilteringPredictionsSingleEndScatter(filterSimulationResults:dict):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    trimPositions = tuple(filterSimulationResults.keys())
    maxExpectedErrors = tuple(filterSimulationResults[trimPositions[0]].keys())
    x = []
    y = []
    z = []
    usingPercentages = False
    for trimPosition in trimPositions:
        for maxExpectedError in maxExpectedErrors:
            x.append(trimPosition)
            y.append(maxExpectedError)
            zValue = filterSimulationResults[trimPosition][maxExpectedError]
            if zValue > 0 and zValue < 1:
                usingPercentages = True
            z.append(zValue)
    xmax = max(x)
    ymax = max(y)
    xmin = min(x)
    ymin = min(y)

    #making plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.scatter(x, y, z, alpha=0.8)
    ax.set_xlabel("Trim Position")
    ax.set_ylabel("Max Expected Error")
    if usingPercentages:
        ax.set_zlabel("Percent remaining reads")
    else:
        ax.set_zlabel("Remaining reads")

    plt.show()

''' Non-working/incomplete plotting functions
def plotFilteringPredictionsSingleEndChecker(filterSimulationResults:dict):
    import numpy
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    trimPositions = tuple(filterSimulationResults.keys())
    maxExpectedErrors = tuple(filterSimulationResults[trimPositions[0]].keys())
    x = []
    y = []
    usingPercentages = False
    for trimPosition in trimPositions:
        for maxExpectedError in maxExpectedErrors:
            x.append(trimPosition)
            y.append(maxExpectedError)
            zValue = filterSimulationResults[trimPosition][maxExpectedError]
            if zValue > 0 and zValue < 1:
                usingPercentages = True
            z.append(zValue)
    xmax = max(x)
    ymax = max(y)
    xmin = min(x)
    ymin = min(y)
    divisions = 50
    xdiv = (xmax - xmin) / divisions
    ydiv = (ymax - ymin) / divisions
    xMesh = numpy.arange(xmin, xmax, xdiv)
    yMesh = numpy.arange(ymin, ymax, ydiv)
    xMeshLen = len(xMesh)
    yMeshLen = len(yMesh)
    xGrid, yGrid = numpy.meshgrid(xMesh, yMesh)

    #mapping color squares
    colorTuple = ('g', 'y')
    xGridShape = xGrid.shape
    yGridShape = yGrid.shape
    colorGrid = numpy.empty(xGridShape, dtype=str)
    for yVal in range(yMeshLen):
        for xVal in range(xMeshLen):
            colorGrid[xVal, yVal] = colorTuple[(xVal + yVal) % len(colorTuple)]

    #making a figure
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surface = ax.plot_surface(xGrid, yGrid, z, facecolors=colorGrid, linewidth=0)

    #setting the Z-axis
    zmin = min(z)
    zmax = max(z)
    ax.set_zlim(zmin, zmax)
    ax.w_zadis.set_major_locator(matplotlib.ticker.LinearLocator(6))

    plt.show()

def plotFilteringPredictionsSingleEndContour(filterSimulationResults:dict):
    import numpy
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    trimPositions = tuple(filterSimulationResults.keys())
    maxExpectedErrors = tuple(filterSimulationResults[trimPositions[0]].keys())
    x = []
    y = []
    z = []
    usingPercentages = False
    for trimPosition in trimPositions:
        for maxExpectedError in maxExpectedErrors:
            x.append(trimPosition)
            y.append(maxExpectedError)
            zValue = filterSimulationResults[trimPosition][maxExpectedError]
            if zValue > 0 and zValue < 1:
                usingPercentages = True
            z.append(zValue)
    xmax = max(x)
    ymax = max(y)
    xmin = min(x)
    ymin = min(y)

    #making contour plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')

def plotFilteringPredictionsSingleEndBar(filterSimulationResults:dict):
    import numpy
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    import matplotlib.ticker
    trimPositions = tuple(filterSimulationResults.keys())
    maxExpectedErrors = tuple(filterSimulationResults[trimPositions[0]].keys())
    x = []
    y = []
    z = []
    usingPercentages = False
    for trimPosition in trimPositions:
        for maxExpectedError in maxExpectedErrors:
            x.append(trimPosition)
            y.append(maxExpectedError)
            zValue = filterSimulationResults[trimPosition][maxExpectedError]
            if zValue > 0 and zValue < 1:
                usingPercentages = True
            z.append(zValue)
    xmax = max(x)
    ymax = max(y)
    xmin = min(x)
    ymin = min(y)

    #making plot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection = '3d')
    ax.bar(x, y, z, alpha=0.8)
    ax.set_xlabel("Trim Position")
    ax.set_ylabel("Max Expected Error")
    ax.set_zlabel("Remaining reads")

    plt.show()
'''

