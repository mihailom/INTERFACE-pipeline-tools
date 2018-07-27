import argparse
import csv
import os
'''
This script will take in a bed file and create a matrix based on the values
of the first read.
Input will be a .BED file
'''


def parseArguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="Input the path of the file you want to create a matrix for")
    parser.add_argument("-s", "--start", type=int,
                        help="Starting coulumn for the bed file. Starting from 1.")
    parser.add_argument("-e", "--end", type=int,
                        help="Ending coulumn for the bed file. Starting from 1.")
    parser.add_argument("-n", "--name", type=int,
                        help="Name coulumn for the bed file. Starting from 1.")
    parser.add_argument("-o", "--output", type=str,
                        help="Output name of the file. Default value will be the name of the input file.")
    parser.add_argument("-ends", "--endsOnly", type=bool,
                        help="Tells the script to only count end frequencies.")
    args = parser.parse_args()
    if (args.endsOnly == None):
        args.endsOnly = False
    return args


def readFile(infile):
    # return the row as a generator object
    try:
        with open(infile, "r") as bedfile:
            for row in bedfile:
                yield row.split()

    except Exception as e:
        print("error {} occured".format(e))


class Matrix(object):

    def __init__(self, nameRow=1, startRow=2, endRow=3, ends=False):
        self.nameRow = nameRow - 1
        self.startRow = startRow - 1
        self.endRow = endRow - 1
        self.ends = ends
        self.matrix = {}
        self.max_len = 1
        self.gene_names = []  # used to call in order

    # vector addition of two arrays
    def extend_vector(self, name, start, end):
        current_vector = self.matrix[name]  # grab the current vector
        if len(current_vector) < self.max_len:
            # add zeros to the vector
            extended_vector = [0 for i in range(self.max_len - len(current_vector))]
            current_vector += extended_vector  # extend the size of the current vector

        # add only the end to the vector
        # EX:
        # Probe-1 ID 1 4 etc... ==> [0,0,0,1]
        # Probe-1 ID 1 3 etc... ==> [0,0,1,1]
        # Probe-1 ID 1 4 etc... ==> [0,0,1,2]
        if self.ends == True:
            current_vector[end - 1] += 1
        elif self.ends == False:   # I know this is ugly but Im just trying to finish fast OK
            # use this vector to check the values in between
            new_vector = list(range(start - 1, end))
            for count, value in enumerate(current_vector):
                if count in new_vector:
                    current_vector[count] += 1

    def convertBadRead(self, value):
        # There was a read that had a value of 82 have 1000 bytes
        toReturn = ""
        try:
            toReturn = int(value)
        except ValueError as e:
            for val in value:
                if ord(val) == 0:
                    break
                toReturn += val
        return int(toReturn)

    def add(self, row):
        # add a probe to the set of matricies
        # currently place holders for test, change if needed
        # name is a string and start and end are ints
        name = row[self.nameRow].strip()
        start = self.convertBadRead(row[self.startRow])
        end = self.convertBadRead(row[self.endRow])
        if end > self.max_len:
            self.max_len = end
        if name in self.matrix:
            # add to the vector
            self.extend_vector(name, start, end)
        else:
            if self.ends:
                sRNA_vector = [0 for i in range(self.max_len)]
                sRNA_vector[end - 1] += 1
            else:
                # create an entry for the sRNA where the list is the max size
                # new vector to add to the matrix dictionary
                new_vector = list(range(start - 1, end))
                sRNA_vector = [1 if i in new_vector else 0 for i in range(self.max_len)]
            self.gene_names.append(name)
            self.matrix[name] = sRNA_vector

    # add zeros to the row
    def normalizeMatrixRow(self, row):
        if len(row) < self.max_len:
            # add zeros to the vector
            extended_vector = [0 for i in range(self.max_len - len(row))]
            row += extended_vector  # extend the size of the
        return row

    def matrixToCSV(self, filename):
        with open("./{}.csv".format(filename), "w") as CSVfile:
            header = "gene,"
            positions = ",".join([str(i) for i in range(1, self.max_len + 1)])
            header += positions + "\n"
            CSVfile.write(header)
            for gene in self.gene_names:
                line = gene + ","
                row = self.matrix[gene]
                row = self.normalizeMatrixRow(row)
                linePositions = ",".join([str(i) for i in row])
                line += linePositions + "\n"
                CSVfile.write(line)

    def __str__(self):
        returnString = ""
        for value in self.gene_names:
            returnString += "{}:{}\n".format(value, self.matrix[value])
        return returnString


def setOutput(outputName, filePath):
    if outputName == None:
        outputName = filePath.split("/")
        outputName = outputName[-1]
        outputName = os.path.splitext(outputName)[0]
    return outputName


def main():
    args = parseArguments()
    filePath, start, end, name, outputName, endsOnly = args.path, args.start, args.end, args.name, args.output, args.endsOnly
    # print(filePath, name, start, end, outputName)
    outputName = setOutput(outputName, filePath)

    print("creating matrix using parameters:\nfilepath: {}\nOutput Name: {}.csv".format(filePath, outputName))
    bedfile = readFile(filePath)
    matrix = Matrix(name, start, end, endsOnly)
    for count, row in enumerate(bedfile):
        # print(count)
        matrix.add(row)

    matrix.matrixToCSV(outputName)
    # print(matrix)


if __name__ == '__main__':
    main()
