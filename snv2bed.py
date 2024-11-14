def calculate_variant_coordinates(ref_allele, alt_alleles, chrom, pos):
    """
    Calculate the BED format coordinates for deletions, insertions (only one), and substitutions.
    
    Parameters:
    - ref_allele: String representing the reference allele.
    - alt_alleles: List of strings representing the alternate alleles.
    - chrom: String representing the chromosome.
    - pos: Integer representing the 1-based position in the VCF file.
    
    Returns:
    - A list of dictionaries with BED format coordinates and the alternate allele for each variant.
    """
    ref_length = len(ref_allele)
    bed_entries = []

    # VCF position is 1-based, but BED format uses 0-based coordinates
    bed_start = pos - 1
    insertion_handled = False  # Track if an insertion has already been processed

    for alt in alt_alleles:
        alt_length = len(alt)
        
        if alt_length < ref_length:  # Deletion
            # The 'start' for the deletion should be +1 to skip the base preceding the deletion
            bed_deletion_start = bed_start + 1
            bed_end = bed_start + ref_length  # The 'end' position where the deletion ends
            bed_entries.append({
                "chrom": chrom,
                "start": bed_deletion_start,
                "end": bed_end,
                "alt": alt,
                "type": "deletion"
            })
        elif alt_length > ref_length and not insertion_handled:  # Insertion (only process the first one)
            # For insertion, the 'end' position should be one base beyond the 'start'
            bed_end = bed_start + 1  # BED end is exclusive, so end = start + 1
            bed_entries.append({
                "chrom": chrom,
                "start": bed_start,
                "end": bed_end,
                "alt": alt,
                "type": "insertion"
            })
            insertion_handled = True  # Mark that we have processed one insertion
        elif alt_length == ref_length:  # Substitution (SNV)
            # For substitution, the 'start' and 'end' positions are the same, covering the exact changed position
            bed_end = bed_start + ref_length  # Substitution spans the same length as reference
            bed_entries.append({
                "chrom": chrom,
                "start": bed_start,
                "end": bed_end,
                "alt": alt,
                "type": "substitution"
            })
    
    return bed_entries

def process_vcf_line_to_bed(vcf_line):
    """
    Process a single VCF line to extract deletions, insertions (only one), and substitutions, 
    and output them in BED format.
    
    Parameters:
    - vcf_line: String representing a single line of the VCF file.
    
    Returns:
    - A list of BED format dictionaries for each alternate allele.
    """
    # Split the VCF line by tab to extract fields
    fields = vcf_line.strip().split("\t")
    
    if len(fields) < 5:
        return []  # Not a valid VCF line
    
    chrom = fields[0]  # Chromosome
    pos = int(fields[1])  # 1-based position
    ref_allele = fields[3]  # Reference allele
    alt_alleles = fields[4].split(",")  # Alternate alleles (separated by commas)
    
    # Calculate the BED format for deletions/insertions/substitutions
    return calculate_variant_coordinates(ref_allele, alt_alleles, chrom, pos)

def format_bed_entry(bed_entry):
    """
    Format a single BED entry as a string in BED format.
    
    Parameters:
    - bed_entry: A dictionary representing a single BED entry.
    
    Returns:
    - A string representing the BED entry in 'chrom start end alt' format.
    """
    return f"{bed_entry['chrom']}\t{bed_entry['start']}\t{bed_entry['end']}"


def main(): # Example usage with a VCF line
    #vcf_line = "chr9\t12479735\t.\tAAAAAT\tAAAAATAAAATAAAAT,A\t40\t.\tAF=0.277778,0.111111;AQ=40,40\tGT:DP:AD:GQ:PL:RNC\t./1:5:0,3,0:1:40,5,0,990,990,990:D.\t0/1:30:16,8,0:4:26,2,0,990,990,990:..\t./.:0:0,0:1:0,0,0,0,0,0:DD\t./.:15:15,0,0:0:0,45,449,45,449,449:II\t0/1:71:57,8,0:18:25,0,77,990,990,990:..\t0/2:10:6,0,4:17:26,990,990,0,990,60:..\t0/2:22:9,0,13:16:40,990,990,0,990,28:..\t0/0:28:28,0,0:50:0,63,869,63,869,869:..\t./.:0:0,0,0:1:0,0,0,0,0,0:DD"
    import sys
    for vcf_line in sys.stdin:
        if vcf_line.startswith("#"):
            continue

        bed_entries = process_vcf_line_to_bed(vcf_line)

# Print all BED entries for the VCF line
        for entry in bed_entries:
            print(format_bed_entry(entry))
if __name__ == '__main__':
	main()