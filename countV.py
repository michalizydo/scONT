# Script to count variants in permutation test that were crossed with exons.
# Usage python countV.py

def countTE(TEfi):
	with open(TEfi,'r') as f:
		var=[]
		ins=[]
		dl=[]
		for line in f:
			liner=tuple(line.rstrip("\n").split("\t"))
			var.append(liner[3])
			if "Sniffles2.INS" in line:
				ins.append(liner[3])
			elif "Sniffles2.DEL" in line:
				dl.append(liner[3])
	f.close()
	return set(var),set(ins),set(dl)
def main():
	import glob
	result=[]
	for alf in glob.glob('shuffles_exon/*.bed'):
		var,ins,dl=countTE(alf)
		result.append(alf.lstrip("exon_cross/")+"\t"+str(len(var))+"\t"+str(len(ins))+"\t"+str(len(dl)))
	result=sorted(result)
	for i in tuple(result):
		print(i)
if __name__=="__main__":
	main()
