# Script to extract SC only, SC+bulk and bulk-only variants from unfiltered vcf.
# Usage cat allMerge.vcf | grep <variant_type> | python nulltest.py <bulk_vcf_index> <sc_vcf_index>
def main():
	import sys
	for line in sys.stdin:
		if line.startswith("#"):
			continue
		liner=tuple(line.split("\t")[9:])

		if "NULL" not in liner[int(sys.argv[1])] and "NULL" not in liner[int(sys.argv[2])]:
			with open(sys.argv[3]+"TP.vcf",'a') as f:
				f.write(line)
			f.close() 
		elif "NULL" in liner[int(sys.argv[1])] and "NULL"  not in liner[int(sys.argv[2])]:
			with open(sys.argv[3]+"FP.vcf",'a') as f:
				f.write(line)
			f.close() 	
		elif "NULL" not  in liner[int(sys.argv[1])] and "NULL" in liner[int(sys.argv[2])]:
			with open(sys.argv[3]+"FN.vcf",'a') as f:
				f.write(line)
			f.close() 	
if __name__=="__main__":
	main()