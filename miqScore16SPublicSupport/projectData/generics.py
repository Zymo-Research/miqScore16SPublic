class InputFile(object):

    def __init__(self, filePath:str):
        import os
        filePath = os.path.abspath(filePath)
        if not os.path.isfile(filePath):
            raise FileNotFoundError("Unable to find input file at %s" %filePath)
        self.filePath = filePath
        self.parseInputFile()

    def parseInputFile(self):
        raise RuntimeError("This was always meant to be overridden and should not be getting hit during the program.  This is a bug.")