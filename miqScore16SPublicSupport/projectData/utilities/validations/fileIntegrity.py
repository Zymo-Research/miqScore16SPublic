import os
import logging
logger = logging.getLogger(__name__)

def md5File(path:str, dontWarnOnEmptyFile = False):
    import hashlib
    path = os.path.abspath(path)
    if not os.path.isfile(path):
        logger.critical("Unable to find file for hashing at %s" %path)
        raise FileNotFoundError("File not found: %s" %path)
    md5 = hashlib.md5()
    with open(path, 'rb') as file:
        for chunk in iter(lambda: file.read(4096), b""):
            md5.update(chunk)
    digest = md5.hexdigest()
    if digest == "d41d8cd98f00b204e9800998ecf8427e" and not dontWarnOnEmptyFile:
        logger.warning("MD5 value d41d8cd98f00b204e9800998ecf8427e is the hash of an empty file. Someone should probably review this. File path: %s" %path)
    logger.info("MD5 checksum for %s was %s" %(path, digest))
    return digest

def getFileSize(path:str):
    if not os.path.isfile(path):
        return None
    size = os.path.getsize(path)
    return size

def logFileInfo(path:str):
    md5 = md5File(path)
    size = getFileSize(path)
    logger.info("File integrity info for %s: MD5=%s SIZE=%s" %(path, md5, size))
    return (md5, size)