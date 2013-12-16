#!/usr/bin/env python

import sys
from glob import glob
import os.path
import md5

#===============================================================================
# Functions
#===============================================================================
def getChecksums(fileNames):
    return dict([ (os.path.basename(f), checksum(f)) for f in fileNames ])

def checksum(fileName):
    return md5.new(file(fileName).read()).hexdigest()

def readChecksums(fileName):
    checksums = dict()
    if os.path.exists(fileName):
        for line in file(fileName):
            checkFile, checksumValue = line.strip().split("\t")
            checksums[checkFile] = checksumValue
    return checksums    

def writeChecksums(checksums, fileName):
    outFile = file(fileName, 'w')
    for fileName, value in checksums.items():
        print >>outFile, "%s\t%s" % (fileName, value)
    outFile.close()

#===============================================================================
# Main
#===============================================================================
path = sys.argv[1]
checksumFileName = "%s/checksums" % path
dataFiles = set(glob("%s/data/*.txt" % path))
dataFiles.discard("%s/data/all.txt" % path)
checksums = getChecksums(dataFiles)
previousChecksums = readChecksums(checksumFileName)
if checksums != previousChecksums:
    writeChecksums(checksums, checksumFileName)
