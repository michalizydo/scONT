import sys
for i in sys.stdin:
	if i.startswith("#"):
		print(i,end='')
	else:
		if float(i.rstrip("\n").split("\t")[-1].split(":")[-2])>2:
			print(i,end='')
