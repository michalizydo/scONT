# Script for determining SNV types for statistical characterization of SNVs and testing presence of MDA footprint.
#
# snvtype.py <REF_allele> <ALT_allele> input.vcf.gz
def checksnv(snvt,file):
	import gzip
	snvt=tuple(snvt)
	with gzip.open(file,'r') as f:
		for line in f:
			line=line.decode("UTF-8")
			if line.startswith("#"):
				continue
			liner=tuple(line.split("\t"))
			if liner[3]==snvt[0] and liner[4][0]==snvt[1]:
				print(line,end='')
def main():
	import sys
	checksnv([sys.argv[1],sys.argv[2]],sys.argv[3])
	

if __name__=="__main__":
	main()
