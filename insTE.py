    """
    Usage:
        python script.py <vcf_file> <sample_index1> <sample_index2>
        
    Arguments:
        - vcf_file: Path to the VCF file to analyze.
        - sample_index1: Index of the SC sample in the VCF file.
        - sample_index2: Index of the bulk sample in the VCF file.
        
    Input:
        Repeatmasker .out file into stdin
          
    Output:
        Prints categorized counts of insertions.
    """

def countSCb(vcffile):
    """
    Reads a VCF file and identifies insertions annotated as "Sniffles2.INS" for specific samples.
    
    Args:
        vcffile (str): Path to the VCF file to be analyzed.

    Returns:
        tuple: 
            - set(SConly): A set of insertion names present only in the SC sample.
            - set(SCbulk): A set of insertion names present in both SC and bulk samples.
            - total (int): The total count of "Sniffles2.INS" insertions across all analyzed samples.

    Note:
        The sample indices to be analyzed are provided as command-line arguments.
    """
    import sys

    SConly = []
    SCbulk = []
    with open(vcffile, 'r') as vcff:
        total = 0
        for line in vcff:
            # Skip header lines in the VCF file
            if line.startswith("#"):
                continue
            else:
                line = tuple(line.split("\t"))
                name = line[2]  # Extract the name of the variant (3rd column)
                test = tuple(line[9:])  # Extract genotype information from the 10th column onward
                
                # Check for "Sniffles2.INS" in specified samples
                if "Sniffles2.INS" in test[int(sys.argv[2])] and "Sniffles2.INS" in test[int(sys.argv[3])]:
                    SCbulk.append(name)
                    total += 1
                elif "Sniffles2.INS" in test[int(sys.argv[3])]:
                    SConly.append(name)
                    total += 1

    return set(SConly), set(SCbulk), total


def getloc(line):
    """
    Splits a line into a tuple of fields.

    Args:
        line (str): A single line from input.

    Returns:
        tuple: A tuple containing the split fields of the line.
    """
    return tuple(line.split())


def main():
    """
    Main function to analyze insertions in a VCF file and categorize them 
    by presence in "SC" and "SC + bulk" groups.
    
    Input:
        - A VCF file path provided as the first command-line argument.
        - Two sample indices provided as the second and third command-line arguments.
        - Additional input is read from stdin in tabular format.

    Output:
        - The counts of Alu and LINE1 insertions in SC-only and SC+bulk groups.
        - The total number of insertions analyzed.
    """
    import sys

    resA = []  # Stores Alu insertions from stdin input
    resL = []  # Stores LINE1 insertions from stdin input

    for i in sys.stdin:
        i = getloc(i)
        if len(i) > 0:
            # Check if the line corresponds to an Alu or LINE1 insertion
            if i[0][0].isdigit() == 1:  # Ensures the first character is a digit (chromosome number)
                if i[10] == "SINE/Alu" and int(i[5]) <= 1000 and 1000 <= int(i[6]):
                    resA.append(i[4].split("_")[0])
                elif i[10] == "LINE/L1" and int(i[5]) <= 1000 and 1000 <= int(i[6]):
                    resL.append(i[4].split("_")[0])

    # Remove duplicates
    resA = list(set(resA))
    resL = list(set(resL))

    # Initialize counts
    alusc = 0
    alubulksc = 0
    linesc = 0
    linebulksc = 0

    # Analyze SC and SC+bulk insertions
    isSC = countSCb(sys.argv[1])
    for i in resA:
        if i in isSC[0]:
            alusc += 1
        if i in isSC[1]:
            alubulksc += 1
    for j in resL:
        if j in isSC[0]:
            linesc += 1
        if j in isSC[1]:
            linebulksc += 1

    # Print results
    print("SC + bulk Alu: " + str(alubulksc))
    print("SC Alu: " + str(alusc))
    print("SC + bulk LINE: " + str(linebulksc))
    print("SC LINE: " + str(linesc))
    print("Total number of INS: " + str(isSC[2]))


if __name__ == "__main__":

    main()
