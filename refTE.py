    """
    Usage:
       cat <repeatmasker_insertions.out> python refTE.py <repeatmasker_insertion_reference_loci.out>  <sample_index1> <sample_index2> <vcf_file>
    
    Arguments:
        repeatmasker_insertions.out: repeatmasker .out file
        repeatmasker_insertion_reference_loci.out: repeatmasker .out file for +-1kb sequences of insertion breakpoints from the reference.
        
        sample_index1: Index of SC sample in the VCF file.
        sample_index2: Index of bulk sample in the VCF file.
        vcf_file: VCF containing insertions
    """

def countSCb(vcffile):
    """
    Reads a VCF file and identifies "Sniffles2.INS" insertions in specific samples.
    
    Args:
        vcffile (str): Path to the VCF file.
    
    Returns:
        tuple:
            - set(SConly): A set of insertion names present only in the SC sample.
            - set(SCbulk): A set of insertion names present in both SC and bulk samples.
            - total (int): Total count of "Sniffles2.INS" insertions.
    """
    import sys

    SConly = []
    SCbulk = []
    with open(vcffile, 'r') as vcff:
        total = 0
        for line in vcff:
            if line.startswith("#"):  # Skip header lines
                continue
            else:
                line = tuple(line.split("\t"))
                name = line[2]  # Variant name
                test = tuple(line[9:])  # Sample genotype fields

                # Check if "Sniffles2.INS" is in the specified sample indices
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
        line (str): A line from input.
    
    Returns:
        tuple: A tuple of split fields.
    """
    return tuple(line.split())


def main():
    """
    Main function to analyze insertions based on SC and SC+bulk conditions,
    and compare them with reference data.
    
    Input:
        - stdin: Tabular input with element data.
        - sys.argv:
            1. Reference data file path.
            2. VCF file path.
            3. Sample index 1.
            4. Sample index 2.
    
    Output:
        Prints counts and comparisons for "same element" and "same type of element".
    """
    import sys

    records = []  # Stores processed records from stdin

    # Process stdin input
    for i in sys.stdin:
        i = getloc(i)
        if len(i) > 0:
            if i[0][0].isdigit() == 1:  # Ensure line starts with a digit (chromosome)
                if i[10] == "SINE/Alu":
                    if len(records) == 0:
                        records.append([i[4], [i[9]], []])
                    else:
                        found = 0
                        for j in range(len(records)):
                            if records[j][0] == i[4]:
                                found = 1
                                records[j][1] = list(set(records[j][1] + [i[9]]))
                                break
                        if not found:
                            records.append([i[4], [i[9]], []])
                elif i[10] == "LINE/L1":
                    if len(records) == 0:
                        records.append([i[4], [], [i[9]]])
                    else:
                        found = 0
                        for j in range(len(records)):
                            if records[j][0] == i[4]:
                                found = 1
                                records[j][2] = list(set(records[j][2] + [i[9]]))
                                break
                        if not found:
                            records.append([i[4], [], [i[9]]])

    # Read reference data
    with open(sys.argv[1], 'r') as f:
        refel = []
        for line in f:
            if line != "\n":
                if tuple(line.split())[0][0].isdigit() == 1:
                    line = tuple(line.split())
                    name = line[4]
                    loc = set(range(int(line[5]), int(line[6]) + 1))
                    if 1000 in loc and 1001 in loc:
                        refel.append(tuple([name, line[9], line[10]]))
        refel = tuple(refel)

    # Analyze against SC and SC+bulk data
    isSC = countSCb(sys.argv[4])
    typesame = 0
    namesame = 0

    # SC only analysis
    print("SC only:")
    for j in refel:
        if j[0].split("_")[0] not in isSC[0]:
            continue
        for i in records:
            if j[0].split("_")[0] == i[0]:
                if "LINE/L1" in j[2] and len(i[2]) > 0:
                    typesame += 1
                    if j[1] in set(i[2]):
                        namesame += 1
                elif "Alu" in j[2] and len(i[1]) > 0:
                    typesame += 1
                    if j[1] in set(i[1]):
                        namesame += 1
    print("Same element: " + str(namesame))
    print("Same type of element: " + str(typesame))
    print("Total number of INS: " + str(isSC[2]))

    # Reset counts for SC+bulk analysis
    typesame = 0
    namesame = 0

    # SC+bulk analysis
    print("SC + bulk:")
    for j in refel:
        if j[0].split("_")[0] not in isSC[1]:
            continue
        for i in records:
            if j[0].split("_")[0] == i[0]:
                if "LINE/L1" in j[2] and len(i[2]) > 0:
                    typesame += 1
                    if j[1] in set(i[2]):
                        namesame += 1
                elif "Alu" in j[2] and len(i[1]) > 0:
                    typesame += 1
                    if j[1] in set(i[1]):
                        namesame += 1
    print("Same element: " + str(namesame))
    print("Same type of element: " + str(typesame))
    print("Total number of INS: " + str(isSC[2]))


if __name__ == "__main__":

    main()
