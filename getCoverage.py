# Script to summarize coverage across genome from Mosdepth output
# Usage python getCoverage.py <mosdepth_coverage.bed.gz>
# Function to calculate the coverage of genomic regions
def calsum(source, totallen, chromos=[], chromolen=[]):
    # If chromosomes are provided, initialize a list to track coverage for each chromosome
    if len(chromos) > 0:
        covnums = []  # List to store coverage for each chromosome
        for i in tuple(chromos):
            covnums.append(0)  # Initialize each chromosome's coverage to 0
    coveredlen = 0  # Initialize the total covered length to 0
    chromos = tuple(chromos)  # Convert chromos to a tuple for faster lookup

    import gzip  # Import the gzip module for reading compressed files

    # Case when no specific chromosomes are given
    if len(chromos) == 0:
        with gzip.open(source, 'r') as b:  # Open the gzipped source file for reading
            for line in b:
                line = line.decode("utf-8")  # Decode each line from bytes to string
                line = tuple(line.split("\t"))  # Split the line by tab and convert it to a tuple
                if line[-1] == "0\n":  # Skip lines where the last column is '0' (no data)
                    continue
                # Calculate the covered length for the region (difference between start and end positions)
                coveredlen += (int(line[2]) - int(line[1]))

        return str(round(coveredlen / totallen, 2))  # Return the coverage as a percentage of the total length

    # Case when specific chromosomes are provided
    else:
        with gzip.open(source, 'r') as b:  # Open the gzipped source file for reading
            for line in b:
                line = line.decode("utf-8")  # Decode each line from bytes to string
                line = tuple(line.split("\t"))  # Split the line by tab and convert it to a tuple
                if line[-1] == "0\n":  # Skip lines where the last column is '0' (no data)
                    continue
                if line[0] in set(chromos):  # If the chromosome is in the specified list
                    # Calculate the covered length for this region and add it to the total
                    coveredlen += (int(line[2]) - int(line[1]))
                    # Add the covered length to the corresponding chromosome's coverage
                    covnums[chromos.index(line[0])] += (int(line[2]) - int(line[1]))

        # Calculate the coverage percentage for each chromosome and store it in a list
        rezo = []
        for i in range(len(covnums)):
            rezo.append(str(round(covnums[i] * 100 / chromolen[i], 2)) + "%")  # Percentage of coverage for each chromosome
        
        # Return the total coverage percentage and the per-chromosome coverage percentages
        return str(round(coveredlen * 100 / totallen, 2)) + "%", tuple(rezo)

# Main function to execute the script
def main():
    import sys  # Import the sys module to access command line arguments
    # Call calsum with the necessary arguments: file path, total length, list of chromosomes, and list of chromosome lengths
    a = calsum(sys.argv[1], 3088269832, ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"],
               [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415])
    
    # Print the total coverage percentage
    print(a[0])

# Entry point of the script
if __name__ == '__main__':
    main()
