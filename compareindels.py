## Script to extract small ins numbers from indels
## Modify lines 104-115 with <input.vcf.gz> <bulk_index> <SC_index> <Indexes_of_SCs> <name>

# Function to process VCF files and count specific occurrences
def countsup(fil, bulkl, scl, scloc, namesc):
    import gzip  # For reading gzipped files
    
    # Initialize counters
    bulk = 0  # Indicates if 'bulk' condition is met
    sc = 0  # Indicates if 'sc' condition is met
    bulkonly = 0  # Count for entries exclusive to 'bulk'
    sconly = 0  # Count for entries exclusive to 'sc'
    bulksc = 0  # Count for entries in both 'bulk' and 'sc'
    multi = [0, 0, 0, 0, 0, 0]  # Counts for multi-level 'sc' events (1sc to 6sc)
    inbulk = 0  # Insertion count for 'bulk'
    insc = 0  # Insertion count for 'sc'
    inbulksc = 0  # Insertion count for overlapping 'bulk' and 'sc'
    
    # Open the gzipped VCF file
    with gzip.open(fil, 'r') as f:
        for line in f:
            # Decode the line from binary to string
            line = line.decode("UTF-8")
            
            # Skip header lines
            if line.startswith("#"):
                continue
            
            # Initialize flags for the current line
            sc = 0
            bulk = 0
            
            # Parse the line into columns
            liner = tuple(line.rstrip("\n").split("\t"))
            
            # Extract relevant sample columns (9th to 17th)
            line = liner[9:18]
            
            countl = 0  # Track column index
            ml = -1  # Multi-level tracker
            
            # Process each sample column
            for i in line:
                loco = tuple(i.split(":"))  # Split by ":" for genotype details
                
                # Check if the genotype is valid
                if loco[0] != "./.":
                    # Check for depth > 2 in alternate allele depths
                    for l in loco[2].split(",")[1:]:
                        if int(l) > 2:
                            # Set 'bulk' or 'sc' flag based on the column
                            if countl == bulkl:
                                bulk = 1
                            elif countl == scl:
                                sc = 1
                countl += 1  # Increment column index
            
            # Analyze scenarios based on 'bulk' and 'sc' flags
            if sc == 1 and bulk == 0:
                # Specific to 'sc' only
                for i in range(len(line)):
                    if i in set(scloc):  # Check if column index matches given locations
                        loco = tuple(line[i].split(":"))
                        # Ensure valid depth
                        if loco[2] == "." or loco[0] == "./.":
                            continue
                        if int(loco[2].split(",")[1]) > 0:
                            ml += 1
                multi[ml] += 1  # Increment appropriate multi-level counter
            
            if bulk == 1 and sc == 1:
                # Both 'bulk' and 'sc'
                bulksc += 1
                # Check for insertions in this case
                if "," in liner[4]:
                    if any(len(j) < len(liner[3]) for j in tuple(liner[4].split(","))):
                        inbulksc += 1
                elif len(liner[4]) > 1:
                    inbulksc += 1
            elif bulk == 1 and sc == 0:
                # Specific to 'bulk' only
                bulkonly += 1
                if "," in liner[4]:
                    if any(len(j) < len(liner[3]) for j in tuple(liner[4].split(","))):
                        inbulk += 1
                elif len(liner[4]) > 1:
                    inbulk += 1
            elif sc == 1 and bulk == 0:
                # Specific to 'sc' only
                sconly += 1
                if "," in liner[4]:
                    if any(len(j) < len(liner[3]) for j in tuple(liner[4].split(","))):
                        insc += 1
                elif len(liner[4]) > 1:
                    insc += 1

    # Return results for the current file
    return namesc, bulkonly, sconly, bulksc, multi, inbulksc, inbulk, insc

# Main function to process multiple VCF files
def main():
    print("name\tFN\tFP\tTP\t1sc\t2sc\t3sc\t4sc\t5sc\t6sc\tinsBulkSC\tinsBulk\tinsSc")
    # Call countsup for each file and print results
    a, b, c, d, e, f, g, h = countsup("PromMinMerges.SNV.minion2.indel.vcf.gz", 4, 0, [0, 1, 2, 3, 7, 8], "minion2")
    print(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + " " + str(g) + " " + str(h))
    a, b, c, d, e, f, g, h = countsup("PromMinMerges.SNV.minion45.indel.vcf.gz", 4, 1, [0, 1, 2, 3, 7, 8], "minion45")
    print(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + " " + str(g) + " " + str(h))
    a, b, c, d, e, f, g, h = countsup("PromMinMerges.SNV.minion67c.indel.vcf.gz", 5, 2, [0, 1, 2, 3, 7, 8], "minion67C")
    print(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + " " + str(g) + " " + str(h))
    a, b, c, d, e, f, g, h = countsup("PromMinMerges.SNV.minion67m.indel.vcf.gz", 6, 3, [0, 1, 2, 3, 7, 8], "minion67M")
    print(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + " " + str(g) + " " + str(h))
    a, b, c, d, e, f, g, h = countsup("PromMinMerges.SNV.promethionC.indel.vcf.gz", 5, 8, [0, 1, 2, 3, 7, 8], "promethionC")
    print(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + " " + str(g) + " " + str(h))
    a, b, c, d, e, f, g, h = countsup("PromMinMerges.SNV.promethionM.indel.vcf.gz", 6, 7, [0, 1, 2, 3, 7, 8], "promethionM")
    print(str(a) + " " + str(b) + " " + str(c) + " " + str(d) + " " + str(e) + " " + str(f) + " " + str(g) + " " + str(h))

# Run the main function
if __name__ == "__main__":
    main()