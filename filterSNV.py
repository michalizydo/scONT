def countsup(fil,bulkl,scl,scloc,namesc):
	import gzip
	bulk=0
	sc=0
	bulkonly=0
	sconly=0
	bulksc=0
	mosaicTP=0
	mosaicFN=0
	multi=[0,0,0,0,0,0]

	with gzip.open(fil,'r') as f:
		for line in f:
			line=line.decode("UTF-8")
			if line.startswith("#"):
				continue
			else:
				sc=0
				bulk=0
				mosaic=0
				line=line.rstrip("\n").split("\t")[9:18]
				
				countl=0
				ml=-1
				valid=0
				#print(line)
				for i in line:
					if i.startswith("./."):
						continue
					test=tuple(i.split(":"))[2].split(",")
					for t in range(len(test)):
						if test[t]==".":
							test[t]=0
					if any(int(dv)>2 for dv in tuple(test)):
						valid=1
						break
					"""
					loco=tuple(i.split(":"))
					if loco[0]!="0/0" and loco[0]!="./.":
						#print(loco[0])
						for l in loco[2].split(",")[1:]:
							if int(l)>2:
								if countl==bulkl:
									bulk=1
								elif countl==scl:
									sc=1
					countl+=1
					"""
				if valid==1:
					for i in line:
						loco=tuple(i.split(":"))
						for l in loco[2].split(",")[1:]:
							if int(l)!=0:
								if countl==bulkl:
									bulk=1
									if loco[0]=="0/1":
										mosaic=1
								elif countl==scl:
									sc=1
								
						countl+=1

				else:
					continue
				if sc==1 and bulk==0:
					for i in range(len(line)):
						if i in set(scloc):
							loco=tuple(line[i].split(":"))
							if  loco[0]=="./.":
								continue
							if int(loco[2].split(",")[1])>0:	
								ml+=1
					multi[ml]+=1	
				if bulk==1 and sc==1:
					bulksc+=1
					if mosaic==1:
						mosaicTP+=1
						
				elif bulk==1 and sc==0:
					bulkonly+=1
					if mosaic==1:
						mosaicFN+=1
				elif sc==1 and bulk==0:
					sconly+=1		
	#print(bulkonly)
	#print(sconly))
	#print(set(bulksc))
	#print(set(multi))											
	return namesc,bulkonly,sconly,bulksc,mosaicTP,mosaicFN,multi
def main():
	a,b,c,d,e,f,g=countsup("PromMinMerges.SNV.minion2.exons.vcf.gz",4,0,[4,1,2,3,7,8],"minion2")
	print(str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" "+str(e)+" "+str(f)+" "+str(g))
	a,b,c,d,e,f,g=countsup("PromMinMerges.SNV.minion45.exons.vcf.gz",4,1,[4,1,2,3,7,8],"minion45")
	print(str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" "+str(e)+" "+str(f)+" "+str(g))
	a,b,c,d,e,f,g=countsup("PromMinMerges.SNV.minion67c.exons.vcf.gz",5,2,[4,1,2,3,7,8],"minion67C")
	print(str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" "+str(e)+" "+str(f)+" "+str(g))
	a,b,c,d,e,f,g=countsup("PromMinMerges.SNV.minion67m.exons.vcf.gz",6,3,[4,1,2,3,7,8],"minion67M")
	print(str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" "+str(e)+" "+str(f)+" "+str(g))
	a,b,c,d,e,f,g=countsup("PromMinMerges.SNV.promethionC.exons.vcf.gz",5,8,[4,1,2,3,7,8],"promethionC")
	print(str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" "+str(e)+" "+str(f)+" "+str(g))
	a,b,c,d,e,f,g=countsup("PromMinMerges.SNV.promethionM.exons.vcf.gz",6,7,[4,1,2,3,7,8],"promethionM")
	print(str(a)+" "+str(b)+" "+str(c)+" "+str(d)+" "+str(e)+" "+str(f)+" "+str(g))
if __name__=="__main__":
	main()
