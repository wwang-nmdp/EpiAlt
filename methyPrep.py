import sys
import gzip
import os
from pathlib import Path
import statistics

from os import listdir
from os.path import isfile, join
counter = 0


def readfile(inputfile, folder):
    with gzip.open(inputfile, 'rt') as f:
        file = f.read()
        print("read "+ inputfile)
    tempFile = open("temp.txt", "w")
    tempFile.write(file)
    tempFile.close()

    print("temp file is generated")
    inputFile = open("temp.txt", "r")
    for line in inputFile:
        #print(line)
        global counter
        counter= counter + 1
        if counter % 1000000 == 0:
            print("processed "+ str(counter) + "lines")

        word = line.split()
        if len(word) != 6:
            continue
        else:
            chrome_data = ChromeDataHolder(word[0], word[1], float(word[3]), int(word[4]), int(word[5]), folder)
            me = int(word[4])
            un = int(word[5])
            if(me + un < checkSum):
                # print("skip because less than "+ str(checkSum))
                continue
            else:
                if chrome_data.getKey() in my_map:
                    # print("merge into existing data")
                    existData = my_map.get(chrome_data.getKey())
                    existData.merge(chrome_data)
                else:
                    # print("add new data")
                    my_map[chrome_data.getKey()] = chrome_data
    inputFile.close()
    print("input file closed")
    if os.path.exists("temp.txt"):
        os.remove("temp.txt")
        print("remove temp file")
    print("read file finished")

class ChromeDataHolder:

    def __init__(self, chrome, location, ratio, meth, unmeth, folderName):
        self.chrome = chrome
        self.location = location
        self.ratio = meth/(meth + unmeth)
        self.meth = meth
        self.unmeth = unmeth
        self.counter = 1
        self.folderName = folderName
        self.folderMap = {folderName: Data(meth, unmeth)}
        self.ratioList = []
        self.ratioList.append(ratio)

    def getKey(self):
        return self.chrome + self.location

    def merge(self,  another_chrome):
        self.meth = self.meth + another_chrome.meth
        self.unmeth = self.unmeth + another_chrome.unmeth
        self.ratio = self.meth / (self.meth + self.unmeth)
        self.counter = self.counter + 1
        self.ratioList.append(another_chrome.ratio)
        if another_chrome.folderName in self.folderMap.keys():
            anotherData = Data(another_chrome.meth, another_chrome.unmeth)
            self.folderMap[another_chrome.folderName].merge(anotherData)
        else:
            self.folderMap[another_chrome.folderName] = another_chrome.folderMap[another_chrome.folderName]
        # print("meth is" + str(self.meth))
        # print("unmeth is" + str(self.unmeth))
        # print("ratio is" + str(self.ratio))


    def shouldPrint(self):
        for key in self.folderMap:
            if(self.folderMap[key].getCount() <= minimumCount):
                return False
        return True

    def toString(self):
        folderCount = ""
        for key in sorted(self.folderMap):
            folderCount = folderCount + str(self.folderMap[key].getRatio()) + "\t" + str(self.folderMap[key].getCount()) + "\t"

        return self.chrome + "\t" + self.location + "\t"  + self.location+ "\t" + str(self.ratio) + "\t"  + \
               str(int(self.meth/self.counter)) + "\t" + str(int(self.unmeth/self.counter)) + "\t" + str(self.counter) + "\t"+ folderCount + "\t"+ str(self.getStdv())+ "\n"

    def getStdv(self):
        return statistics.stdev(self.ratioList)

class Data:
    def __init__(self, meth, unmeth):
        self.meth = meth
        self.unmeth = unmeth
        self.count = 1

    def merge(self, anotherData):
        self.meth = self.meth + anotherData.meth
        self.unmeth = self.unmeth + anotherData.unmeth
        self.count = self.count + 1

    def getRatio(self):
        return self.meth/(self.meth + self.unmeth)

    def getCount(self):
        return self.count

def processFolder(folderPath):
    #Get folder name
    path = Path(folderPath)
    folderName = path.name
    print("process folder "+ folderName)
    #Process files in the folder
    onlyfiles = [f for f in listdir(folderPath) if isfile(join(folderPath, f))]
    for file in onlyfiles:
        print(file)
        readfile(folderPath + "//" + file, folderName)


# First param is folder paths, which is separated with ;
# Second param is output file
# Third param is check sum for meth
# Fourth param is minimum count for each folder
mypaths = sys.argv[1]
folders = mypaths.split("+")
output = sys.argv[2]
checkSum = int(sys.argv[3])
minimumCount = int(sys.argv[4])
my_map = {}
outputFile = open(output, "w")

#Process files in each folder
for folder in folders:
    processFolder(folder)

#generate output
print("Generating output file")
hasTitle = False
for key in sorted(my_map):
    if my_map.get(key).shouldPrint():
        if not hasTitle:
            folders = ""
            for folder in sorted(my_map.get(key).folderMap):
                folders = folders + folder + "_ratio\t" + folder+"_count\t"
            outputFile.write("chrome\tstart\tend\tratio\tmeth\tunmeth\tsum\t" + folders + "\t STDV"+"\n")
            hasTitle = True

        outputFile.write(my_map.get(key).toString())

outputFile.close()


