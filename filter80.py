# A script that counts up results from "Onecodetofindthemall"
#
# Usage: python filter80.py <fraction_of_target_TE_ref>

def filtereighty(data,threshold):
	if data.startswith("###"): 

		if data.split("\t")[16].rstrip("\n")!="No_ref_available":
			val=float(data.split("\t")[16].rstrip("\n"))
			if val>=threshold:
				print(data,end='')
def main():
	import sys
	for i in sys.stdin:
		filtereighty(i,float(sys.argv[1]))

if __name__=="__main__":
	main()
