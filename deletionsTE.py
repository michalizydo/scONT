# Script to summarize TE-TE deletion events.
#  vcf indexes start with 0 for first sample.
#
# Usage cat <deletions_with_flanks_repeatmasker.out> | python deletionsTE.py <index_bulk_in_vcf> <index_sc_in_vcf> <ins+del_calls_vcf> <insertions_SC_only_repeatmasker.out> <insertions_SC_bulk_repeatmasker.out> <insertions_bulk_only_repeatmasker.out> <0 or 1 for SC-only search (0) or SC+bulk search (1)>

def countSCb(vcffile):
    import sys
    SConly=[]
    SCbulk=[]
    bulk=[]
    bulkins=[]
    with open(vcffile,'r') as vcff:
        total=0
        for line in vcff:
            if line.startswith("#"):
                continue
            else:

                line=tuple(line.split("\t"))
                name=line[2]
                test=tuple(line[9:])
#both       
                if "Sniffles2.DEL" in test[int(sys.argv[1])] and "Sniffles2.DEL" in test[int(sys.argv[2])]:
                        SCbulk.append(name)
                        total+=1
#SC only
                elif "Sniffles2.DEL" in test[int(sys.argv[2])]:
                        SConly.append(name)
                        total+=1
#bulk only
                elif "Sniffles2.DEL" in test[int(sys.argv[1])]:
                        bulk.append(name)
                        total+=1
                
#bulk ins       
                if "Sniffles2.INS" in test[int(sys.argv[1])] and "Sniffles2.INS" not in test[int(sys.argv[2])] :
                        bulkins.append(name)
                        
    vcff.close()
    return set(SConly),set(SCbulk),total,set(bulk),set(bulkins)

def getloc(line):
        line=tuple(line.split())
        #if len(line)==0:
        #   return None 
        #elif line[0].isdigit()==False:
        return line
def main():
    records=[]
    import sys
    for i in sys.stdin:
        i=getloc(i)
        if len(i)>0:
            if i[0][0].isdigit()==1:
                if i[10]=="SINE/Alu" and ((1000 in set(range(int(i[5]),int(i[6])+1))) or (int(i[7].lstrip("(").rstrip(")"))<1000 and 1000 in range(int(i[7].lstrip("(").rstrip(")")),int(i[7].lstrip("(").rstrip(")"))+int(i[6])-int(i[5])))):
                    if len(records)==0:
                        if (int(i[7].lstrip("(").rstrip(")"))<1000 and 1000 in range(int(i[7].lstrip("(").rstrip(")")),int(i[7].lstrip("(").rstrip(")"))+int(i[6])-int(i[5]))):
                            records.append([i[4],[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"],[]])
                        else:
                            records.append([i[4],[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]],[]])
                    else:
                        found=0

                        for j in range(len(records)):
                            if records[j][0]==i[4]:
                                found=1
                                if (int(i[7].lstrip("(").rstrip(")"))<1000 and 1000 in range(int(i[7].lstrip("(").rstrip(")")),int(i[7].lstrip("(").rstrip(")"))+int(i[6])-int(i[5]))):
                                    if i[-1]=="*":
                                        records[j][1]=list(set(records[j][1]+[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]))
                                    else:
                                        records[j][1]=list(set(records[j][1]+[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]))

                                else:
                                    if i[-1]=="*":
                                        records[j][1]=list(set(records[j][1]+[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]]))
                                    else:
                                        records[j][1]=list(set(records[j][1]+[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]]))


                                #   records[j][1]=list(set(records[j][1]+[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]))
                                #else:
                                #   records[j][1]=list(set(records[j][1]+[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]]))
                                break
                            if found==1:
                                break
                        if found==0:
                            if i[-1]=="*":
                                records.append([i[4],[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]],[]])
                            else:
                                records.append([i[4],[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]],[]])
                elif i[10]=="LINE/L1" and ((1000 in set(range(int(i[5]),int(i[6])+1))) or (int(i[7].lstrip("(").rstrip(")"))<1000 and 1000 in range(int(i[7].lstrip("(").rstrip(")")),int(i[7].lstrip("(").rstrip(")"))+int(i[6])-int(i[5])))):
                    if len(records)==0:
                        if (1000 in range(int(i[7].lstrip("(").rstrip(")")),int(i[7].lstrip("(").rstrip(")"))+int(i[6])-int(i[5]))):
                            if i[-1]=="*":
                                records.append([i[4],[],[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]])
                            else:
                                records.append([i[4],[],[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]])
                        else:
                            if i[-1]=="*":
                                records.append([i[4],[],[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]]])
                            else:
                                records.append([i[4],[],[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]])
                    else:
                        found=0

                        for j in range(len(records)):
                            if records[j][0]==i[4]:
                                found=1
                                if (1000 in range(int(i[7].lstrip("(").rstrip(")")),int(i[7].lstrip("(").rstrip(")"))+int(i[6])-int(i[5]))):
                                    if i[-1]=="*":
                                        records[j][2]=list(set(records[j][2]+[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]))
                                    else:
                                        records[j][2]=list(set(records[j][2]+[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]+"*"]))

                                else:
                                    if i[-1]=="*":
                                        records[j][2]=list(set(records[j][2]+[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]]))
                                    else:
                                        records[j][2]=list(set(records[j][2]+[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]]))

                                break
                            if found==1:
                                break
                        if found==0:
                            if i[-1]=="*":
                                records.append([i[4],[],[i[-2].rstrip("\n")+"_"+i[6]+"_"+i[9]]])
                            else:
                                records.append([i[4],[],[i[-1].rstrip("\n")+"_"+i[6]+"_"+i[9]]])
    isSC=countSCb(sys.argv[3])
    bulkins=[]
    with open(sys.argv[4],'r') as f:
        for line in f:
            line=tuple(line.split())
            if len(line)>0:
                if line[0].isdigit()==1:
                    name=line[4].split("_")[2]
                    #print(name)
                    if name in isSC[4]:
                        bulkins.append(line)
    f.close()
    with open(sys.argv[5],'r') as f:
        for line in f:
            line=tuple(line.split())
            if len(line)>0:
                if line[0].isdigit()==1:
                    name=line[4].split("_")[2]
                    #print(name)
                    if name in isSC[4]:
                        bulkins.append(line)
    f.close()
    with open(sys.argv[6],'r') as f:
        for line in f:
            line=tuple(line.split())
            if len(line)>0:
                if line[0].isdigit()==1:
                    name=line[4].split("_")[2]
                    #print(name)
                    if name in isSC[4]:
                        bulkins.append(line)
    f.close()

    alualu=0
    aluinalu=0
    lineline=0
    lineinline=0
    #
    # Change below to isSC[1] to look in bulk+SC, isSC[0] to look SC only
    # 
    TEs=isSC[int(sys.argv[7])]
    for i in records:
        #print(i,end=" ")
        blist=[]
        for b in bulkins:
            if b[4].split("_")[0]==i[0].split("_")[1].split(":")[0] and ( int(i[0].split("_")[1].split(":")[1].split("-")[0])<int(b[4].split("_")[1])<int(i[0].split("_")[1].split(":")[1].split("-")[0])+2000 or int(i[0].split("_")[1].split(":")[1].split("-")[1])-2000<int(b[4].split("_")[1])<int(i[0].split("_")[1].split(":")[1].split("-")[1])):
                blist.append(b)
        #if len(blist)>0:
        #   print(i)
        #   print(blist)
        if i[0].split("_")[0] in TEs:
        #   print("somatic")
            if len(i[1])>1:
                aster=0
                noaster=0
                for j in i[1]:
                    if "*" in j:
                        aster+=1
                    else:
                        noaster+=1
                if aster>0 and noaster>0:
                    alualu+=1
                    print("alu alu\t",end='')
                    print(i)
                elif (aster==0 and noaster>0) or (aster>0 and noaster==0):
                    if len(blist)>0:
                        for b in tuple(blist):
                            if b[10]=="SINE/Alu":
                                aluinalu+=1
                                print("alu in alu\t",end='')
                                print(b,end=' ')
                                print(i)
            elif len(blist)>0:
                for b in tuple(blist):
                    if b[10]=="SINE/Alu" and i[-2]!=[]: 
                        aluinalu+=1
                        print("alu in alu\t",end='')
                        print(b,end=' ')
                        print(i)
            if len(i[2])>1:
                aster=0
                noaster=0
                for j in i[2]:
                    if "*" in j:
                        aster+=1
                    else:
                        noaster+=1
                if aster>0 and noaster>0:
                    lineline+=1
                    print("line line\t",end='')
                    print(i)
                elif (aster==0 and noaster>0) or (aster>0 and noaster==0):
                    if len(blist)>0:
                        for b in tuple(blist):
                            if b[10]=="LINE/L1":
                                lineinline+=1
                                print("line in line\t",end='')
                                print(b,end=' ')
                                print(i)
                                
            elif len(blist)>0:
                for b in tuple(blist):
                    if b[10]=="LINE/L1" and i[-1]!=[]:
                        lineinline+=1
                        print("line in line\t",end='')
                        print(b,end=' ')
                        print(i)
#       else:
#           print("")
    print("Somatic results:\n \t\t\tAlu-Alu on del junctions: "+str(alualu))
    print("\t\t LINE1-LINE1 on del junctions: "+str(lineline))
    print("\t\t Alu-(Novel Alu ins in bulk) on del junctions: "+str(aluinalu))
    print("\t\t LINE1-(Novel LINE1 ins in bulk) on del junctions: "+str(lineinline))
    print("\t\t Total somatic deletions: "+str(len(set(TEs))))
    
    print("\n\n Total deletions: "+str(isSC[2]))
if __name__=="__main__":
    main()