# Script to count SINE/Alu and LINE/L1 results from permutation test.
# Usage python countTEs.py
# Shuffles must be in shuffles_aluline/ directory

def countTE(TEfi):
	with open(TEfi,'r') as f:
		alu=[]
		l=[]
		for line in f:
			liner=tuple(line.rstrip("\n").split("\t"))
			if 'LINE/L1' in line:
				l.append(liner[3])
			elif 'SINE/Alu' in line:
				alu.append(liner[3])
	f.close()
	return set(alu),set(l)
def main():
	import glob
	result=[]
	for alf in glob.glob('shuffles_aluline/*.bed'):
		alu,ll=countTE(alf)
		result.append(alf.lstrip("shuffles/")+"\tAlu: "+str(len(alu))+"\tLINE: "+str(len(ll)))
	result=sorted(result)
	for i in tuple(result):
		print(i)
if __name__=="__main__":
	main()