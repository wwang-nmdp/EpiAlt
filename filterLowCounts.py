import sys
import gzip
import os
input = sys.argv[1]
output = sys.argv[2]
checkSum = int(sys.argv[3])

with gzip.open(input, 'rt') as f:
    file = f.read()
tempFile = open("temp.txt", "w")
outputFile = open(output, "w")
tempFile.write(file)
tempFile.close()

inputFile = open("temp.txt", "r")
for line in inputFile:
    print(line)
    word = line.split()
    me = int(word[4])
    un = int(word[5])
    if(me + un >= checkSum):
        outputFile.write(line)
inputFile.close()
outputFile.close()
if os.path.exists("temp.txt"):
    os.remove("temp.txt")