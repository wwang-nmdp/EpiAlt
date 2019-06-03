import sys
import gzip
import os
from pathlib import Path
import statistics

from os import listdir
from os.path import isfile, join

class Range:
    def __init__(self, chrome, start, end, symbol):
        self.chrome = chrome
        self.start = start
        self.end = end
        self.symbol = symbol

class methData:
    def __init__(self, chrome, location, meth, umeth):
        self.chrome = chrome
        self.location = location
        self.meth = meth
        self.umeth = umeth


def processFolder(folderPath):
    #Get folder name
    path = Path(folderPath)
    folderName = path.name
    print("process folder "+ folderName)
    #Process files in the folder
    onlyfiles = [f for f in listdir(folderPath) if isfile(join(folderPath, f))]
    for file in onlyfiles:
        readfile(folderPath + "//" + file)

def readfile(input):
    my_map = {}
    inputfile = open(input, "r")
    for line in inputfile:
        #print(line)
        words = line.split()
        chrome = words[0]
        location = words[1]
        meth = words[4]
        umeth = words[5]
        meth_data = methData(chrome, location, meth, umeth)
        if chrome in my_map:
            list = my_map.get(chrome)
            list.append(meth_data)
        else:
            list = []
            list.append(meth_data)
            my_map[chrome] = list
    print("Finish read "+input)
    print("Find symbols in "+input)

    input_path = Path(input)
    output = input_path.name
    for symbol in rangelist:
        if symbol.chrome in my_map:
            count = 0
            meth = 0
            umeth = 0
            genelist = my_map.get(symbol.chrome)
            for gene in genelist:
                if int(symbol.start) <= int(gene.location) <= int(symbol.end):
                    count = count + 1
                    meth = meth + 1
                    umeth = umeth + 1
            if count > 0:
                output = output + "\t" + str(count)
                output = output + "\t" + str(meth/(meth + umeth))
            else:
                output = output + "\t" + str(count)
                output = output + "\t" + "NA"
        else:
            output = output + "\t" + "0"
            output = output + "\t" + "NA"
    outputlins.append(output + "\n")


def readRange(input):
    inputFile = open(input, "r")
    #Skip the column names
    inputFile.readline()
    title = "fileName"
    for line in inputFile:
        #print(line)
        words = line.split()
        range = Range(words[0], words[1], words[2], words[3])
        rangelist.append(range)
        title = title + "\t" + words[3]+"_count" + "\t"+words[3]+"_ratio"

    print("finish read range")
    outputlins.append(title + "\n")



# First param is a file defines the chromosome coordinates to check methylation level
# Second param is folder path with cov
# Third param is output file
rangeFile = sys.argv[1]
mypath = sys.argv[2]
output = sys.argv[3]
rangelist = []
outputlins = []
readRange(rangeFile)
processFolder(mypath)
outputFile = open(output, "w")
for line in outputlins:
    outputFile.write(line)
outputFile.close()