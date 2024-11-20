# Script to add CIPOS/CIEND tags to vcf, requiered by SVtyper

def main():
	import sys
	for line in sys.stdin:
		if line.startswith("#"):
			print(line,end="")
		else:
			line=line.split("\t")
			line[7]=line[7]+";CIPOS=-100,100;CIEND=-100,100"
			print("\t".join(line),end='')
	
if __name__ == '__main__':
	main()
