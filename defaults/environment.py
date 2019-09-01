import os
import datetime
import __main__
timestamp = str(datetime.datetime.now().timestamp()).replace(".", "")
rScriptExecutable = "/usr/local/bin/Rscript"
dataFolder = "/data"
projectFolder = os.path.split(os.path.split(os.path.abspath(__file__))[0])[0]
inputFolder = os.path.join(dataFolder, "input")
outputFolder = os.path.join(dataFolder, "output")
sequenceFolder = os.path.join(inputFolder, "sequence")
forwardReads = os.path.join(sequenceFolder, "standard_submitted_R1.fastq")
reverseReads = os.path.join(sequenceFolder, "standard_submitted_R2.fastq")
rScriptFolder = os.path.join(projectFolder, "rscripts")
databaseFile = os.path.join(os.path.split(__main__.__file__)[0], "reference", "rdp_train_set_16.fa.gz")
logFile = os.path.join(outputFolder, "dada2.%s.log" %timestamp)
