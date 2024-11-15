# Script filtering the single cell data for PASS and >2 support.
#
# Usage: filterSVmerge.py bulk_vcf_limited_with_SC_coverage_bed.vcf.gz bulk_bam_file_name.bam SC_vcf_limited_with_SC_coverage_bed.vcf.gz SC_bam_file_name.bam merged_vcf_name.vcf.gz
def openMerge(mergename,svd):
	import gzip
	with gzip.open(mergename,'r') as f:
		for line in f:
			line=line.decode('UTF-8')
			if line.startswith("#CHROM"):
				filenames=tuple(line.rstrip("\n").split("FORMAT	")[1].split("\t"))
			if line.startswith("#"):
				print(line,end='')
				continue

			presence=line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").split("\t")
			#print(presence)
			sup=[]
			for i in range(len(presence)):
				
				id=presence[i].split(":")[-1]
				if id!="NULL":
					result=[]
					found=0
					for S in svd:
						if S[0]==filenames[i]:
							for p in presence[i].split(":")[-1].split(","):	
								if p in S[1]:
									found=1
									result.append(p)
									break
						if found==1:
							break
					if found==1:
						sup.append(int(presence[i].split(":")[-2]))
									
# filter for at least support 2.
			
			if any(s > 2 for s in sup):
				#print(sup)
				sections=tuple(line.split("\t"))
				if sections[6]=="PASS" and (sections[7].split(";")[1].lstrip("SVTYPE=") in {"INS","DEL"}):
					print(line,end='')

			
				
def main():
	import gzip
	from os import walk
	import sys
	SVdata=[]
	filenames = next(walk('.'), (None, None, []))[2] 
	for i in filenames:
		if i.endswith(sys.argv[1]):
			with gzip.open(i,'r') as f:
				svs=[]
				for line in f:
					line=line.decode("UTF-8")
					if line.startswith("chr"):
						line=tuple(line.split("\t"))
						svs.append(line[2])
				
			SVdata.append(tuple([sys.argv[2],set(svs)]))
	for i in filenames:
		if i==sys.argv[2]:
			with gzip.open(i,'r') as f:
				svs=[]
				for line in f:
					line=line.decode("UTF-8")
					if line.startswith("chr"):
						line=tuple(line.split("\t"))
						svs.append(line[2])
			SVdata.append(tuple([sys.argv[4],set(svs)]))
	openMerge(sys.argv[5],tuple(SVdata))

if __name__=="__main__":
    main()
