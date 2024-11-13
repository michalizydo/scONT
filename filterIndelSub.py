import sys
for i in sys.stdin:
	if i.startswith("#"):
		print(i,end='')
	else:
		j=tuple(i.split("\t"))
		# get substitutions
		if len(j[3])==1 and len(j[4])==1:
			print(i,end='')
		# get small InDels - change above to: 
		#	continue
		# else:
		#	print(i,end='')
