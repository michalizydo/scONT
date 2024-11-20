# Script to filter SVs by support and coverage from VCF tags, used on blood LR-WGS data.
# Insert vcf file name in line 29.
# Usage python filtersupcov_blood.py 

def filterMerge(mergename):
	import gzip
	with gzip.open(mergename,'r') as f:
		for line in f:
			line=line.decode('UTF-8')
			if line.startswith("#"):
				print(line,end='')
				continue

			presence=line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").split("\t")

			supf=0
			cov=0
			for p in presence[0:-1]:
				cov+=int(p.split(":")[-3])
				cov+=int(p.split(":")[-2])
			for p in presence:
				if int(p.split(":")[-2])>2:
					supf=1
					break
			sections=tuple(line.split("\t"))
			if sections[6]=="PASS" and (sections[7].split(";")[1].lstrip("SVTYPE=") in {"INS","DEL"}) and cov>4 and supf==1:
				print(line,end='')

def main():
	filterMerge('tester.vcf.gz') # give it your sniffles merge vcf file

if __name__=="__main__":
    main()
