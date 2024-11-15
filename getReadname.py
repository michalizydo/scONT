## Script to trace back SVs from merge into individual single-cells
##
## Usage python getReadname.py vcf_merge_file.vcf.gz vcf_SC_calls.vcf.gz vcf_SC_mosaic_calls.vcf.gz SC_bam_alignment.bam bulk_index_in_vcf_merge sc_index_in_vcf_merge

# Extracting SVs from merged file
def getfilteredSV(vcffile):
	def fil(line):
		liner=tuple(line.rstrip("\n").split("\t"))
		
		ids=[]
		blk=liner[int(sys.argv[5])]
		
		sc=liner[int(sys.argv[6])]
		# Check if variant is present in SC		
		if "Sniffles2.INS" in sc or "Sniffles2.DEL" in sc:
				
			sc=tuple(sc.split(":"))[-1]
			for j in sc.split(","):
				ids.append(j)
		# Mark SC+bulk variants
		if "NULL" not in blk:
			line="!"+line			
						
		return tuple([line,set(ids)])
	# Extract data from vcf.gz SV calls using above procedure
	if vcffile.endswith("gz"):
		import gzip
		with gzip.open(vcffile,'r') as f:
			header=[]
			SV=[]	
			for line in f:
				if line.startswith("#") or ("Sniffles2.INS" not in line and "Sniffles2.DEL" not in line):
					continue
				else:
					out=fil(line)
					SV.append(out)
		f.close()
	else:
		# Extract data from vcf SV calls using above procedure
		with open(vcffile,'r') as f:
			header=[]
			SV=[]	
			for line in f:
				if line.startswith("#") or ("Sniffles2.INS" not in line and "Sniffles2.DEL" not in line):
					continue
				else:
					out=fil(line)
					SV.append(out)
		f.close()
	# Return required SVs and metadata.
	return SV
# Function to extract read names from original vcf.gz calls of SC.
def getReadnames(vcfReadName):
	import gzip
	with gzip.open(vcfReadName,'r') as f:
		readnames=[]
		# Parse all SVs
		for line in f:
			line=line.decode("UTF-8")
			if line.startswith("#"):
				continue
			else:
				liner=tuple(line.split("\t"))
				# Get only INS and DEL - parse RNAMES from calls
				if liner[2].startswith("Sniffles2.INS") or liner[2].startswith("Sniffles2.DEL"):
					readnames.append(tuple([liner[2],tuple(tuple(liner[7].split("RNAMES="))[1].split(";"))[0].split(",")]))

	f.close()
	# Return variants with their RNAMEs.
	return readnames
# Function to parse bam alignment files and extract reads that correspond to variants.
def findRG(bamfile,pairedres,allreadnames):
	import pysam
	samfile = pysam.AlignmentFile(bamfile, "rb")
# Fetch all reads
	for read in samfile.fetch():
		# Select reads that were selected as parts of variants.
		if read.qname in allreadnames:
			# If the read is selected, find variant that is represented by that read
			for j in range(len(pairedres)):
				if read.qname in pairedres[j][1]:
					# If no list for SC file names,attach one
					if len(pairedres[j])==2:
						pairedres[j].append([])
					# Get SC file name from RG tag of the bam merge and append it to the variant.	
					for rt in read.tags:
						if rt[0]=="RG" or rt[0]=="rg":
							pairedres[j][-1].append(read.qname+":"+rt[1])
							break
		test=[]
	# Remove redundant SC file names in lists appended to variants
	for pr in range(len(pairedres)):
		pairedres[pr][-1]=tuple(set(pairedres[pr][-1]))

	return pairedres
def main():
	import sys
	# Get SVs from merge
	res=getfilteredSV(sys.argv[1])

	# Get SVs and their read names from SC vcf calls.
	readnames=getReadnames(sys.argv[2])
	readnames=readnames+getReadnames(sys.argv[3]) # mosaic call vcf

	pairedres=[]
	totalres=[]
	# Extract read names that are searched for and pair with varaint names
	for i in range(len(res)):
		pairedreadnames=[]

		for j in readnames:

			if j[0] in res[i][1]:
				pairedreadnames+=j[1]
		if len(pairedreadnames)>0:
			pairedres.append([res[i][0],set(pairedreadnames)])
			totalres+=list(set(pairedreadnames))
	pairedtest=[]
	# Get SC file names from original bam merge,
	pairedfiles=findRG(sys.argv[4],pairedres,set(totalres))

	
	# Calculate how many variants come from how many single cells.
	statsFPi=[]
	statsTPi=[]
	statsFPd=[]
	statsTPd=[]

	for pf in pairedfiles:
		
		if pf[0].startswith("!"):
			fil=[]
			
			for j in pf[-1]:
				fil.append(j.split(":")[-1])
			if "INS" in pf[0]:
				statsTPi.append(len(set(fil)))
			else:
				statsTPd.append(len(set(fil)))

		else:
			fil=[]
			
			for j in pf[-1]:
				fil.append(j.split(":")[-1])

			if "INS" in pf[0]:
				statsFPi.append(len(set(fil)))
			else:
				statsFPd.append(len(set(fil)))
	# Print results
	print("TP insertions")
	for i in range(max(statsTPi)+1):
		print(str(i)+" "+str(statsTPi.count(i)))
	print("FP insertions")
	for i in range(max(statsFPi)+1):
		print(str(i)+" "+str(statsFPi.count(i)))
	print("TP deletions")
	for i in range(max(statsTPd)+1):
		print(str(i)+" "+str(statsTPd.count(i)))
	print("FP deletions")
	for i in range(max(statsFPd)+1):
		print(str(i)+" "+str(statsFPd.count(i)))
	
if __name__=="__main__":
    main()