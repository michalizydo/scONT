
def filterinsdel(vcf):
	import gzip
	with gzip.open(vcf) as f:
		for line in f:
			if line.startswith(b'#'):
				print(line.decode("UTF-8"),end='')
			elif tuple(tuple(line.decode('UTF-8').split("\t"))[2].split("."))[1] in {"INS","DEL"} and tuple(line.decode('UTF-8').split("\t"))[6]=="PASS":
				sup=0
				for i in tuple(tuple(line.decode('UTF-8').split("\t"))[7].split(";")):
					if i.startswith("SUPPORT="):
						sup=int(i.lstrip("SUPPORT="))
						break
				if sup>2:
					print(line.decode('UTF-8'),end='')
			else:
				continue
	f.close()
def main():
	import sys
	#filter1k(sys.argv[1])
	filterinsdel(sys.argv[1])
if __name__=="__main__":
	main()