# Script to calculate basic stats from detected SVs.
#
# Usage stats_SV.py vcf_file.vcf.gz bulk_index SC_index

def main():
    # Initialize result counters for insertions and deletions:
    # rezM2i/d - for Promethion control vs bulk
    # Format: [True Positive (TP), False Positive (FP), False Negative (FN), Single cell+bulk Mosaic, Bulk-only Mosaic]
    rezM2i = [0, 0, 0, 0, 0]
    rezM2d = [0, 0, 0, 0, 0]
    import sys

    # Open the VCF file for reading
    with open(sys.argv[1], 'r') as f:
        for line in f:
            # Skip header lines
            if line.startswith("#"):
                continue

            # Process non-header lines for SVTYPE field (INS or DEL)
            svtype = line.split("\t")[7].split(";")[1].lstrip("SVTYPE=")
            const = []
            constp = []
            sup = 0  

            # Parse the format fields
            for i in line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").rstrip("\t").split("\t"):
                if "NULL" not in i:
                    const.append(1)
                    # Track maximum support score
                    if int(i.split(":")[-2]) > sup:
                        sup = int(i.split(":")[-2])
                    # Identify partial matches based on GT and DV fields
                    constp.append(1 if i.startswith("0/0") and i.split(":")[3] != "0" else 0)
                else:
                    const.append(0)
                    constp.append(0)

            # Ensure sufficient support score
            if sup < 3:
                print(f"{svtype} SUP ERROR!!!")

            # Analyze insertions
            if svtype == "INS":
                if const[int(sys.argv[2])] == 1 and const[int(sys.argv[3])] == 1:
                    rezM2i[0] += 1  # True Positive
                if constp[int(sys.argv[2])] == 1 and (constp[int(sys.argv[3])] == 1 or const[int(sys.argv[3])] == 1):
                    rezM2i[3] += 1  # Single cell+bulk Mosaic
                if const[int(sys.argv[2])] == 0 and const[int(sys.argv[3])] == 1:
                    rezM2i[1] += 1  # False Positive
                if const[int(sys.argv[2])] == 1 and const[int(sys.argv[3])] == 0:
                    rezM2i[2] += 1  # False Negative
                if constp[int(sys.argv[2])] == 1 and constp[int(sys.argv[3])] == 0 and const[int(sys.argv[3])] == 0:
                    rezM2i[4] += 1  # Bulk-only Mosaic

            # Analyze deletions
            elif svtype == "DEL":
                if const[int(sys.argv[2])] == 1 and const[int(sys.argv[3])] == 1:
                    rezM2d[0] += 1  # True Positive
                if constp[int(sys.argv[2])] == 1 and (constp[int(sys.argv[3])] == 1 or const[int(sys.argv[3])] == 1):
                    rezM2d[3] += 1  # Single cell+bulk Mosaic
                if const[int(sys.argv[2])] == 0 and const[int(sys.argv[3])] == 1:
                    rezM2d[1] += 1  # False Positive
                if const[int(sys.argv[2])] == 1 and const[int(sys.argv[3])] == 0:
                    rezM2d[2] += 1  # False Negative
                if constp[int(sys.argv[2])] == 1 and constp[int(sys.argv[3])] == 0 and const[int(sys.argv[3])] == 0:
                    rezM2d[4] += 1  # Bulk-only Mosaic

    # Calculate and display precision, recall, and F1-score for insertions
    print("## INSERTIONS")
    print("## TEST NAME \tTP\tFP\tFN\tSingle cell+bulk Mosaic\tBulk-only Mosaic\tPRECISION\tRECALL\tF1")
    precision_i = rezM2i[0] / (rezM2i[0] + rezM2i[1])
    recall_i = rezM2i[0] / (rezM2i[0] + rezM2i[2])
    f1_i = (2 * precision_i * recall_i) / (precision_i + recall_i)
    print("INSERTIONS ", rezM2i + [precision_i, recall_i, f1_i])

    # Calculate and display precision, recall, and F1-score for deletions
    print("## DELETIONS")
    print("## TEST NAME \tTP\tFP\tFN\tSingle cell+bulk Mosaic\tBulk-only Mosaic\tPRECISION\tRECALL\tF1")
    precision_d = rezM2d[0] / (rezM2d[0] + rezM2d[1])
    recall_d = rezM2d[0] / (rezM2d[0] + rezM2d[2])
    f1_d = (2 * precision_d * recall_d) / (precision_d + recall_d)
    print("DELETIONS", rezM2d + [precision_d, recall_d, f1_d])


if __name__=="__main__":
    main()