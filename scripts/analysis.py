#Read in file
import sys

mixedFileName = 'results/fasta/mixed_5.txt'
haplogrepFileName = 'results/haplogrep_5.txt'

mixedFile = open(mixedFileName, 'r')
haplogrepFile = open(haplogrepFileName, 'r')
mixedFileLines = mixedFile.read().split("\n")
mixedFileContents = mixedFileLines[2].split(", ")
haplogrepFileContents = haplogrepFile.read().replace(" ", ' ').split("\n")
haplogrepFileContents[0] = haplogrepFileContents[0].replace(' ', '').replace('""', ',').replace('"', '').split(',')
haplogrepFileContents[1] = haplogrepFileContents[1].replace('  ', '').replace('"""', '"" "').replace('" "', ',').replace('"', '').split('\t')
types = haplogrepFileContents[1][6].split(" ")
count = 0
for i in types:
    for j in mixedFileContents:
        if(i == j):
            count += 1

print(count, " out of ", len(mixedFileContents), " found from ", mixedFileLines[1])

mixedFile.close()
haplogrepFile.close()

mixedFile = open(mixedFileName, 'a')
mixedFile.write(str(count))
mixedFile.close()
