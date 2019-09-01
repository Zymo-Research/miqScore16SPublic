from . import applicationSupport
from . import formatReaders


class FastqSamples(object):

    def __init__(self):
        import os
        thisFolder = os.path.split(__file__)[0]
        self.forward = os.path.join(thisFolder, "fileSamples/forward.good.fastq")
        self.reverse = os.path.join(thisFolder, "fileSamples/reverse.good.fastq")
        self.forwardBadHeader = os.path.join(thisFolder, "fileSamples/forward.badHeaderLine.fastq")
        self.forwardMismatchScoreAndSeqLength = os.path.join(thisFolder, "fileSamples/forward.seqQScoreLengthMismatch.fastq")
        self.reverseReadCountMismatch = os.path.join(thisFolder, "fileSamples/reverse.missingRead.fastq")
        self.reverseHeaderMismatch = os.path.join(thisFolder, "fileSamples/reverse.headerMismatch.fastq")
        self.directorySample = os.path.join(thisFolder, "fileSamples/sampleDirectory")

class Dada2_tests(object):

    def test_ExpectedErrorCurveGeneration(self):
        forwardCurve, reverseCurve = applicationSupport.dada2.expectedErrorCurve.calculateExpectedErrorCurvesForFastqFolder(FastqSamples().directorySample, formatReaders.fastq.fileNamingStandards.ZymoServicesNamingStandard, makePNG=True)
        assert abs(forwardCurve.a - 1.82514740818) < 0.00000005
        assert abs(forwardCurve.b - 0.00658991293689) < 0.00000005
        assert abs(forwardCurve.c + 2.42451396627) < 0.00000005
        #large = applicationSupport.dada2.expectedErrorCurve.calculateExpectedErrorCurvesForFastqFolder("c:/Users/mweinstein/dada2big/input/sequence", formatReaders.fastq.fileNamingStandards.ZymoServicesNamingStandard, subsample=200, makeSVG=True)

    def test_ExpectedErrorReadLossPaired(self):
        expectedErrorTries = (
            (1, 2),
            (3, 5),
            (6, 9)
        )
        trimPositionTries = (
            (167, 150),
            (200, 175)
        )
        result = applicationSupport.dada2.dada2FilterPrediction.predictReadsKeptAfterExpectedErrorFilteringPaired(FastqSamples().forward, FastqSamples().reverse, trimPositionTries, expectedErrorTries)
        assert result == {(167, 150): {(1, 2): 10, (6, 9): 20, (3, 5): 18}, (200, 175): {(1, 2): 9, (6, 9): 19, (3, 5): 11}}

    def test_ExpectedErrorReadLossSingle(self):
        expectedErrorTries = (1, 3, 6)
        trimPositionTries = (167, 200)
        result = applicationSupport.dada2.dada2FilterPrediction.predictReadsKeptAfterExpectedErrorFilteringSingle(FastqSamples().forward, trimPositionTries, expectedErrorTries)
        assert result == {200: {1: 9, 3: 12, 6: 19}, 167: {1: 10, 3: 18, 6: 20}}

    def test_PlotExpectedErrorReadLossSingle(self):
        expectedErrorTries = (1, 2, 3, 4, 5, 6, 7, 8)
        trimPositionTries = (150, 175, 200, 225, 250, 275, 300)
        result = applicationSupport.dada2.dada2FilterPrediction.predictReadsKeptAfterExpectedErrorFilteringSingle(FastqSamples().forward, trimPositionTries, expectedErrorTries, percentage=True)
        assert result == {225: {1: 0.35, 2: 0.5, 3: 0.55, 4: 0.55, 5: 0.75, 6: 0.85, 7: 0.9, 8: 1.0}, 200: {1: 0.45, 2: 0.55, 3: 0.6, 4: 0.75, 5: 0.9, 6: 0.95, 7: 1.0, 8: 1.0}, 300: {1: 0.1, 2: 0.25, 3: 0.35, 4: 0.5, 5: 0.5, 6: 0.55, 7: 0.55, 8: 0.6}, 175: {1: 0.5, 2: 0.55, 3: 0.85, 4: 0.95, 5: 1.0, 6: 1.0, 7: 1.0, 8: 1.0}, 275: {1: 0.1, 2: 0.35, 3: 0.5, 4: 0.55, 5: 0.55, 6: 0.55, 7: 0.7, 8: 0.75}, 150: {1: 0.55, 2: 0.7, 3: 0.95, 4: 1.0, 5: 1.0, 6: 1.0, 7: 1.0, 8: 1.0}, 250: {1: 0.3, 2: 0.45, 3: 0.55, 4: 0.55, 5: 0.55, 6: 0.75, 7: 0.75, 8: 0.85}}
        #formatReaders.fastq.filtering.dada2FilterPrediction.plotFilteringPredictionsSingleEndScatter(result) #example code, not for running here


if __name__ == "__main__":
    Dada2_tests().test_ExpectedErrorCurveGeneration()