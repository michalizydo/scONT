# Script to filter variants based on population frequency from STIX annotation.
# Usage: python filterSTIX.py <stix_annotated_SVs.vcf.gz> <fraction_threshold_to_filter_out>

def main():
	import gzip
	import sys
	namez=[]
	with gzip.open(sys.argv[1],'r') as g:
		for line in g:
			if "Sniffles2.".encode("UTF-8") in line:
				liner=tuple(line.decode("UTF-8").split("\t"))
				if "Sniffles2." in liner[2]:
					data=float(line.decode("UTF-8").split("STIX_FREQ=")[1].split("\t")[0])
					if data<=0.01 int(sys.argv[2]):
						namez.append(liner[2])
			
	#print(namez)
	for line in sys.stdin:
		if line.startswith('Sniffles2.'):
			if tuple(line.split('\t'))[0] in set(namez):
				print(line,end='')

		else:
			print(line,end='')
if __name__ == '__main__':
	main()