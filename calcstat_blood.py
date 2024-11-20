# Script to calculate statistics for blood LR-WGS data
# In line 17, insert vcf file name
# In lines 20 and 21 insert bulk and single-cell index, starting at 0
def main():
	rezM2i=[0,0,0]
	rezM4i=[0,0,0]
	rezM6i=[0,0,0]
	rezM6ci=[0,0,0]
	rezPi=[0,0,0]
	rezPci=[0,0,0]

	rezM2d=[0,0,0]
	rezM4d=[0,0,0]
	rezM6d=[0,0,0]
	rezM6cd=[0,0,0]
	rezPd=[0,0,0]
	rezPcd=[0,0,0]
	with open("promethion.MSA2.filtered.v4.bam.5.vcf",'r') as f: # give it your tested vcf here
		
		bulkind=0 # give it index of bulk file as it appears within the vcf tags, first file is 0
		scind=1 # give it index of single cell file as it appears within the vcf tags, first file is 0


		for line in f:
			if line.startswith("#"):
				pass
			else:
				if line.split("\t")[7].split(";")[1].lstrip("SVTYPE=")=="INS":
					const=[]
					sup=0
					#print(line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").rstrip("\t").split("\t"))
					for i in line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").rstrip("\t").split("\t"):
						if "NULL" not in i:
							const.append(1)
							if int(i.split(":")[-2])>sup:
								sup=int(i.split(":")[-2])
						else:
							const.append(0)
					if sup<3:
						print("INS SUP ERROR!!!")
					# bulk vs minion2	
					#print(const)	
					if const[bulkind]==1 and const[scind]==1:
						rezM2i[0]+=1
					if const[bulkind]==0 and const[scind]==1:
						#print(line,end='')
						rezM2i[1]+=1
					if const[bulkind]==1 and const[scind]==0:
						rezM2i[2]+=1
					
				elif line.split("\t")[7].split(";")[1].lstrip("SVTYPE=")=="DEL":
					const=[]
					sup=0
					#print(line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").rstrip("\t").split("\t"))
					for i in line.split("GT:GQ:DR:DV:ID\t")[1].rstrip("\n").rstrip("\t").split("\t"):
						if "NULL" not in i:
							const.append(1)
							if int(i.split(":")[-2])>sup:
								sup=int(i.split(":")[-2])
						else:
							const.append(0)
					if sup<3:
						print("DEL SUP ERROR!!!")
					# bulk vs minion2		
					if const[bulkind]==1 and const[scind]==1:
						rezM2d[0]+=1
					if const[bulkind]==0 and const[scind]==1:
						rezM2d[1]+=1
						#print(line,end='')
					if const[bulkind]==1 and const[scind]==0:
						rezM2d[2]+=1
					
		print("## INSERTIONS")
		print("## SC+B	SConly	Bonly	Fraction on SC found in bulk	Fraction of bulk found in SC	HM")
		
		print(rezM2i+[(rezM2i[0])/(rezM2i[0]+rezM2i[1]),(rezM2i[0])/(rezM2i[0]+rezM2i[2]),((2*(rezM2i[0])/(rezM2i[0]+rezM2i[1])*(rezM2i[0])/(rezM2i[0]+rezM2i[2]))/((rezM2i[0])/(rezM2i[0]+rezM2i[1])+(rezM2i[0])/(rezM2i[0]+rezM2i[2])))])
		print("## DELETIONS")
		print("## SC+B	SConly	Bonly	Fraction on SC found in bulk	Fraction of bulk found in SC	HM")
		
		print(rezM2d+[(rezM2d[0])/(rezM2d[0]+rezM2d[1]),(rezM2d[0])/(rezM2d[0]+rezM2d[2]),((2*(rezM2d[0])/(rezM2d[0]+rezM2d[1])*(rezM2d[0])/(rezM2d[0]+rezM2d[2]))/((rezM2d[0])/(rezM2d[0]+rezM2d[1])+(rezM2d[0])/(rezM2d[0]+rezM2d[2])))])


if __name__=="__main__":
    main()