from . import formatReaders
from . import applicationSupport


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


class FastqReader_tests():

    def test_fileLoad(self):
        fastq = formatReaders.fastq.fastqHandler.FastqFile(FastqSamples().forward)
        fastq.close()

    def test_iterationSingleEnd(self):
        fastq = formatReaders.fastq.fastqHandler.FastqFile(FastqSamples().forward)
        for read in fastq:
            assert type(read) == formatReaders.fastq.fastqHandler.FastqLineSet

    def test_iterationPairdEnd(self):
        fastqPair = formatReaders.fastq.fastqHandler.FastqFilePair(FastqSamples().forward, FastqSamples().reverse)
        for pair in fastqPair:
            assert type(pair) == tuple
            assert len(pair) == 2
            for read in pair:
                assert type(read) == formatReaders.fastq.fastqHandler.FastqLineSet

    def test_ValidationOnValids(self):
        assert formatReaders.fastq.fastqHandler.validFastqFile(FastqSamples().forward)
        assert formatReaders.fastq.fastqHandler.validFastqPair(FastqSamples().forward, FastqSamples().reverse)

    def test_ValidationOnBadHeader(self):
        assert not formatReaders.fastq.fastqHandler.validFastqFile(FastqSamples().forwardBadHeader)

    def test_ValidationOnMismatchedScoreAndSeqLength(self):
        assert not formatReaders.fastq.fastqHandler.validFastqFile(FastqSamples().forwardMismatchScoreAndSeqLength)

    def test_ValidationOnMismatchedHeaderPair(self):
        assert not formatReaders.fastq.fastqHandler.validFastqPair(FastqSamples().forward, FastqSamples().reverseHeaderMismatch)

    def test_ValidationOnReadCountMismatchedPair(self):
        assert not formatReaders.fastq.fastqHandler.validFastqPair(FastqSamples().forward, FastqSamples().reverseReadCountMismatch)

    def test_LengthMeasure(self):
        assert formatReaders.fastq.fastqHandler.estimateReadLength(FastqSamples().forward, getVariance=True) == (321, 0)

    def test_EncodingDetection(self):
        assert formatReaders.fastq.fastqHandler.findQualityScoreEncoding(FastqSamples().forward) == formatReaders.qualityScore.qualityScoreHandler.encodingSchemes.sanger

    def test_FastqPairFinder(self):
        fastqTable = formatReaders.fastq.fastqHandler.getSamplePairTableFromFolder(FastqSamples().directorySample)
        assert len(fastqTable) == 3
        assert len(fastqTable["unpaired"]) == 2
        assert len(fastqTable["sample", "1"]) == 2
        assert len(fastqTable["sample", "2"]) == 2


class FastqAnalysis_tests(object):

    def test_QualityMatrix(self):
        result = formatReaders.fastq.fastqAnalysis.buildQualityMatrix(FastqSamples().forward)
        assert result.sum() == 180673

    def test_ExpectedErrorMatrix(self):
        result = formatReaders.fastq.fastqAnalysis.buildExpectedErrorMatrix(FastqSamples().forward)
        assert int(result.sum()) == 12992

    def test_PercentileCutoff(self):
        result = formatReaders.fastq.fastqAnalysis.findCutoffByPercentile(FastqSamples().forward, 10, 10)
        assert result == 167

    def test_QualityCountMatrix(self):
        result = formatReaders.fastq.fastqAnalysis.makeQualityMatrix(FastqSamples().forward)
        assert result.sum() == 6420



    def test_DirectoryPlotter(self):
        import hashlib
        result = formatReaders.fastq.fastqAnalysis.plotFastqFilesInFolder(FastqSamples().directorySample, formatReaders.fastq.fileNamingStandards.ZymoServicesNamingStandard, base64Output=True, outputFormat="png")
        resultHashes = set()
        expectedSet = {'ff17b3895d9b58cb39bbf1e163be7f09', '91e6a896f1249bb0d743c891117b95c1', 'a5ba0cc981315c26c29b78d7ffb84bf6', '814cb11021efb23a36f976abb5ca52f3'}
        for key in result:
            resultHashes.add(hashlib.md5(result[key]).hexdigest())
        assert resultHashes == expectedSet


if __name__=="__main__":
    fastqTests = FastqReader_tests()