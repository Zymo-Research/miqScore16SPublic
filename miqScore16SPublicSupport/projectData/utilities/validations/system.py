def fileExistsAndAbsolutePath(fileName):
    import os
    if not os.path.isfile(fileName):
        return False
    else:
        return os.path.abspath(fileName)

def directoryExistsAndAbsolutePath(directory):
    import os
    if not os.path.isdir(directory):
        return False
    else:
        return os.path.abspath(directory)