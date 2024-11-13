# Filtering chimeric reads by number of splits and % mismatch to reference
#
# Usage samtools view -H in.bam | python MD_1SA5pp.py <max_number_of_splits> <max_%_mismatch_allowed>

import re


def cigar2(cigar_string):
    # Initialize the read length
    read_length = 0
    # Parse the CIGAR string using regular expressions
    cigar_tokens = re.findall(r'(\d+)([MIDNSHP=X])', cigar_string)
    # Calculate read length
    for length, operation in cigar_tokens:
        if operation in ('M', 'I', 'S', '=', 'X'):
            read_length += int(length)
    
    return read_length
   # Function to count mismatches from MD tag
def count_mismatches(md_tag):
	    mismatch_count = 0
	    i = 0
	    # parse MD tag
	    while i < len(md_tag):
	    	# Count matches
	        if md_tag[i].isdigit():
	            i += 1
	            continue
	        # count mismatches    
	        elif md_tag[i].isalpha():
	            mismatch_count += 1
	            i += 1
	        elif md_tag[i] == '^':
	            # Skip over deletions indicated by '^' followed by the deletion sequence
	            i += 1
	            while i < len(md_tag) and md_tag[i].isalpha():
	                i += 1
	        else:
	        	# error for debugging
	        	raise ValueError(f"Unexpected character in MD tag: {md_tag[i]} \n"+ md_tag)
	    return mismatch_count

def main():
	banned=[]
	import sys
	# parse from stdin
	for line in sys.stdin:
		# pass headers
		if line.startswith("@"):
			print(line,end='')

		else:
			liner=tuple(line.split("\t"))
			# look for splits in line
			if "SA:Z:" in line:
				tomiss=0
				for j in liner:
					if 'SA:Z:' in j:
						if j.count(";")>int(sys.argv[1]):
							tomiss=1
							break
			# discard if number of splits > argument 1
				if tomiss==1:
					continue
			mismatches=-1
			# look for MD tag
			for i in liner[11:]:
				if i.startswith("MD:Z:"):
					mismatches = count_mismatches(i.lstrip("MD:Z:"))
					break
			# discard if mismatch % > argument 2
			if mismatches>-1:
				if mismatches<cigar2(liner[5])*float(sys.argv[2].rstrip("%").rstrip(" "))/100:
					print(line,end='')
			else:
				print(line,end='')	

if __name__=="__main__":
	main()
