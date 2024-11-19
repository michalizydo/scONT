# Script to summarize statistics on SINE/Alu and LINE/L1 elements
# Indexes start at 0 for the first sample in vcf.
# 
# Usage python alulineins.py <bulk_vcf_index> <sc_vcf_index> <vcf_calls.vcf>

def countSCb(vcffile):
    # Initialize lists to hold different categories of names and a total counter.
    bulkOnly = []
    SConly = []
    SCbulk = []
    
    # Open the VCF file for reading.
    with open(vcffile, 'r') as vcff:
        total = 0
        # Iterate through each line in the VCF file.
        for line in vcff:
            # Skip header lines indicated by '#'.
            if line.startswith("#"):
                continue
            else:
                # Split the line into a tuple and extract relevant fields.
                line = tuple(line.split("\t"))
                name = line[2]  # Extract the name (usually the ID field).
                test = tuple(line[9:])  # Extract test results (starting from 9th field).

                # Check for specific conditions to classify the current entry.
                if "Sniffles2.INS" in test[int(sys.argv[1])] and "Sniffles2.INS" in test[int(sys.argv[2])]:
                    SCbulk.append(name)  # Both conditions matched: add to SCbulk.
                    total += 1
                elif "Sniffles2.INS" in test[int(sys.argv[1])] and "Sniffles2.INS" not in test[int(sys.argv[2])]:
                    bulkOnly.append(name)  # Only first condition matched: add to bulkOnly.
                    total += 1
                elif "Sniffles2.INS" in test[int(sys.argv[2])]:
                    SConly.append(name)  # Only second condition matched: add to SConly.
                    total += 1

    # Return unique sets and total count.
    return set(SConly), set(SCbulk), set(bulkOnly), total


def getloc(line):
    # Split the line into a tuple and return it.
    line = tuple(line.split())
    return line


def main():
    records = []  # Initialize a list to hold records.
    import sys  # Import the sys module to handle command-line arguments.

    # Iterate through standard input (assumed to be genomic data).
    for i in sys.stdin:
        i = getloc(i)  # Process the line to extract fields.
        if len(i) > 0:
            # Check if the first field is a digit or a header line.
            if i[0][0].isdigit() == 1 or i[0][0] == "#":
                # Handle 'SINE/Alu' type records.
                if i[10] == "SINE/Alu":
                    if len(records) == 0:
                        # Add a new record if the list is empty.
                        records.append([i[4], [i[9] + "/" + str(i[5] + ":" + i[6] + ":" + i[7].lstrip("(").rstrip(")") + "_" + i[14] + "_" + i[6]], [], [0, 0]])
                    else:
                        found = 0
                        # Search for an existing record by name.
                        for j in range(len(records)):
                            if records[j][0] == i[4]:
                                found = 1
                                # Update the existing record with a new entry.
                                records[j][1] = list(set(records[j][1] + [i[9] + "/" + str(i[5] + ":" + i[6] + ":" + i[7].lstrip("(").rstrip(")") + "_" + i[14] + "_" + i[6]]))
                                break
                            if found == 1:
                                break
                        if found == 0:
                            # Create a new record if none was found.
                            records.append([i[4], [i[9] + "/" + str(i[5] + ":" + i[6] + ":" + i[7].lstrip("(").rstrip(")") + "_" + i[14] + "_" + i[6]], [], [0, 0]])

                # Handle 'LINE/L1' type records similarly.
                elif i[10] == "LINE/L1":
                    if len(records) == 0:
                        records.append([i[4], [], [i[9] + "/" + str(i[5] + ":" + i[6] + ":" + i[7].lstrip("(").rstrip(")") + "_" + i[14] + "_" + i[6]], [0, 0]])
                    else:
                        found = 0
                        for j in range(len(records)):
                            if records[j][0] == i[4]:
                                found = 1
                                records[j][2] = list(set(records[j][2] + [i[9] + "/" + str(i[5] + ":" + i[6] + ":" + i[7].lstrip("(").rstrip(")") + "_" + i[14] + "_" + i[6]]))
                                break
                            if found == 1:
                                break
                        if found == 0:
                            records.append([i[4], [], [i[9] + "/" + str(i[5] + ":" + i[6] + ":" + i[7].lstrip("(").rstrip(")") + "_" + i[14] + "_" + i[6]], [0, 0]])

    # Call the countSCb function to gather statistics from the provided VCF file.
    a = countSCb(sys.argv[3])
    
    # Process each record and print the results.
    for r in range(len(records)):
        print(r)
        linefl = 0  # Counter for LINE records.
        alufl = 0   # Counter for Alu records.

        print(records[r][0])  # Print the name of the record.
        print(records[r][1])  # Print Alu records.
        cat = []  # Category list for unique Alu types.
        
        # Categorize Alu records by type
        for i in records[r][1]:
            cat.append(i.split("_")[1])
        cat = list(set(cat))  # Get unique types.
        
        # Create a structured list based on Alu types
        for j in range(len(cat)):
            cat[j] = [cat[j]]  # Create a list for each category.
            for i in records[r][1]:
                if i.split("_")[1] == cat[j][0]:
                    cat[j].append(i)  # Append relevant records.
            cat[j] = cat[j][1:]  # Remove the category name.

        print(cat)  # Print the categorized Alu records.

        aluranges = 0  # Initialize a range counter for Alu.
        
        # Process Alu records if present.
        if len(cat) > 0:
            inslen = int(cat[0][0].split("/")[1].split("_")[0].split(":")[1]) + int(cat[0][0].split("/")[1].split("_")[0].split(":")[2])
            print("INS LEN: " + str(inslen))  # Print the total length of insertions.
            
            # Find the effective ranges covered by Alu records.
            for i in cat:
                ranges = set([])  # Set to store unique ranges.
                for j in i:
                    ranges.update(set(range(int(j.split("/")[1].split(":")[0]), int(j.split("/")[1].split(":")[1]) + 1)))
                print(len(ranges))  # Print the number of unique ranges.
                alufl += 1  # Increment Alu record count.
                aluranges += len(ranges)  # Add to total range count.
            print("Alu combined: " + str(aluranges))
            
            # If the Alu records cover more than 80% of insertion length, increment the Alu flag.
            if float(aluranges) / float(inslen) > 0.8:
                print("ALU! " + str(alufl))
                records[r][3][0] += alufl  # Update the count of Alu records.

        print(records[r][2])  # Print LINE records.

        cat = []  # Reset category list for LINE records.
        
        # Categorize LINE records by type.
        for i in records[r][2]:
            cat.append(i.split("_")[1])
        cat = list(set(cat))  # Get unique types.
        for j in range(len(cat)):
            cat[j] = [cat[j]]  # Initialize LINE categories.
            for i in records[r][2]:
                if i.split("_")[1] == cat[j][0]:
                    cat[j].append(i)  # Append relevant records.
            cat[j] = cat[j][1:]  # Remove the category name.

        print(cat)  # Print LINE record categories.

        lineranges = 0  # Initialize a range counter for LINE.
        
        # Process LINE records if present.
        if len(cat) > 0:
            inslen = int(cat[0][0].split("/")[1].split("_")[0].split(":")[1]) + int(cat[0][0].split("/")[1].split("_")[0].split(":")[2])
            print("INS LEN: " + str(inslen))  # Print the total length of insertions.
            
            for i in cat:
                ranges = set([])  # Set for unique ranges.
                for j in i:
                    ranges.update(set(range(int(j.split("/")[1].split(":")[0]), int(j.split("/")[1].split(":")[1]) + 1)))
                    
                print(len(ranges))  # Print the number of unique ranges.
                linefl += 1  # Increment LINE record count.
                lineranges += len(ranges)  # Add to total range count.
            print("LINE combined: " + str(lineranges))

            # If the LINE records cover more than 80% of insertion length, increment the LINE flag.
            if float(lineranges) / float(inslen) > 0.8:
                print("LINE! " + str(linefl))
                records[r][3][1] += linefl  # Update the count of LINE records.

    # Initialize various counters for analysis.
    sc_alu = 0
    sc_alue = 0
    sc_alub_fl = 0
    scl = 0
    scle = 0
    sc_lb_fl = 0
    scbl = 0
    scble = 0
    scb_lb_fl = 0
    scb_alu = 0
    scb_alue = 0
    scb_alub_fl = 0
    b_alu = 0
    b_alue = 0
    b_alub_fl = 0
    bl = 0
    ble = 0
    blb_fl = 0

    # Initialize arrays for classification.
    aluSCf = [0, 0, 0, 0]  # Alu classification for SConly (AluJ, AluS, AluY, other).
    alubSCf = [0, 0, 0, 0]  # Alu classification for SCbulk.
    alubf = [0, 0, 0, 0]    # Alu classification for bulk only.

    # Analyze records based on categories obtained from countSCb.
    for i in tuple(records):
        print(i)  # Print the current record for debugging.
        
        # Check if record name is in SConly set.
        if i[0].split("_")[2] in a[0]:
            # If Alu records are present.
            if len(i[1]) > 0:
                sc_alu += 1  # Increment Alu count for SConly.
                sc_alub_fl += i[3][0]  # Update Alu baseflow count.
                sc_alue += len(i[1])  # Update total Alu entries.
                
                # Count types of Alu entries.
                for al in i[1]:
                    if "AluJ" in al:
                        aluSCf[0] += 1
                    elif "AluS" in al:
                        aluSCf[1] += 1
                    elif "AluY" in al:
                        aluSCf[2] += 1
                    else:
                        aluSCf[3] += 1
            
            # If LINE records are present.
            if len(i[2]) > 0:
                scl += 1  # Increment LINE count for SConly.
                sc_lb_fl += i[3][1]  # Update LINE baseflow count.
                scle += len(i[2])  # Update total LINE entries.
        
        # Check if record name is in SCbulk set.
        elif i[0].split("_")[2] in a[1]:
            if len(i[1]) > 0:  # If Alu records are present.
                scb_alu += 1  # Increment Alu count for SCbulk.
                scb_alub_fl += i[3][0]  # Update Alu baseflow count.
                scb_alue += len(i[1])  # Update total Alu entries.
                
                # Count types of Alu entries.
                for al in i[1]:
                    if "AluJ" in al:
                        alubSCf[0] += 1
                    elif "AluS" in al:
                        alubSCf[1] += 1
                    elif "AluY" in al:
                        alubSCf[2] += 1
                    else:
                        alubSCf[3] += 1
            
            # If LINE records are present.
            if len(i[2]) > 0:
                scbl += 1  # Increment LINE count for SCbulk.
                scb_lb_fl += i[3][1]  # Update LINE baseflow count.
                scble += len(i[2])  # Update total LINE entries.
        
        # Check if record name is in bulkOnly set.
        elif i[0].split("_")[2] in a[2]:
            if len(i[1]) > 0:  # If Alu records are present.
                b_alu += 1  # Increment Alu count for bulk only.
                b_alub_fl += i[3][0]  # Update Alu baseflow count.
                b_alue += len(i[1])  # Update total Alu entries.
                
                # Count types of Alu entries.
                for al in i[1]:
                    if "AluJ" in al:
                        alubf[0] += 1
                    elif "AluS" in al:
                        alubf[1] += 1
                    elif "AluY" in al:
                        alubf[2] += 1
                    else:
                        alubf[3] += 1
            
            # If LINE records are present.
            if len(i[2]) > 0:
                bl += 1  # Increment LINE count for bulk only.
                blb_fl += i[3][1]  # Update LINE baseflow count.
                ble += len(i[2])  # Update total LINE entries.

    # Print final summary of counts and statistics.
    print("SC_alu INS: " + str(sc_alu), end=' ')
    print(aluSCf)
    print("SC_line INS: " + str(scl) + "\nSCbulk_alu INS: " + str(scb_alu), end=' ')
    print(alubSCf)
    print("SCbulk_line INS: " + str(scbl) + "\nBulk only alu INS:" + str(b_alu), end=' ')
    print(alubf)
    print("Bulk only line INS:" + str(bl) + "\nTotal merged INS: " + str(a[3]))
    print("SC_alub_TE: " + str(sc_alue) + " SC_lineb_TE: " + str(scle) + " SCbulk_alub_TE: " + str(scb_alue) + " SCbulk_lineb_TE: " + str(scble) + " bulk_alub_TE: " + str(b_alue) + " bulk_lineb_TE: " + str(ble))
    print("SC_alub_fl: " + str(sc_alub_fl) + " SC_lineb_fl: " + str(sc_lb_fl) + " SCbulk_alub_fl: " + str(scb_alub_fl) + " SCbulk_lineb_fl: " + str(scb_lb_fl) + " bulk_alub_fl: " + str(b_alub_fl) + " bulk_lineb_fl: " + str(blb_fl))


if __name__ == "__main__":
    main()  # Call the main function to execute the script.
