# Script to extract allele frequency information for fig2.
#
# Usage python scatterdata.py <Perbrain_SVs_filtered.vcf> <Brainsample1.sc.SVs.vcf.gz> <Brainsample1.sc.SVs.mosaic.vcf.gz> <Brainsample2.sc.SVs.vcf.gz> <Brainsample2.sc.SVs.mosaic.vcf.gz> Brainsample1.sc.bam Brainsample2.sc.bam <bulk_index> <Brainsample1.sc_index> <Brainsample2.sc_index>

import sys
import gzip
import pysam

def getfilteredSV(vcffile):
    def fil(line):
        liner = tuple(line.rstrip("\n").split("\t"))
        
        ids = []
        ids2=[]
        blk = liner[int(sys.argv[8])]  # Allele frequency information stored here
        sc = liner[int(sys.argv[9])]
        sc2 = liner[int(sys.argv[10])]

        if "Sniffles2.INS" in sc or "Sniffles2.DEL" in sc:
            sc = tuple(sc.split(":"))[-1]
            for j in sc.split(","):
                ids.append(j)
        if "Sniffles2.INS" in sc2 or "Sniffles2.DEL" in sc2:
            sc2 = tuple(sc2.split(":"))[-1]
            for j in sc2.split(","):
                ids2.append(j)
        
        af = None
        if "NULL" not in blk:
            # Calculate allele frequency (AF)
            ref_count = int(blk.split(":")[2])
            alt_count = int(blk.split(":")[3])
            af = alt_count / (ref_count + alt_count) if (ref_count + alt_count) > 0 else 0.0
            line = "!" + line  # Mark valid linesls
        
        return tuple([line, set(ids),set(ids2), af])  # Return AF in the tuple
            
    if vcffile.endswith("gz"):
        with gzip.open(vcffile, 'rt') as f:
            SV = []    
            for line in f:
                if line.startswith("#") or ("Sniffles2.INS" not in line and "Sniffles2.DEL" not in line):
                    continue
                else:
                    out = fil(line)
                    SV.append(out)
    else:
        with open(vcffile, 'r') as f:

            SV = []    
            for line in f:
                if line.startswith("#") or ("Sniffles2.INS" not in line and "Sniffles2.DEL" not in line):
                    continue
                else:
                    out = fil(line)
                    SV.append(out)
    return SV

def traceSV(scvcf,mosaicvcf,svlist):
    result=[]
    with gzip.open(scvcf,'r') as f:
        for line in f:
            line=line.decode("UTF-8")
            for sv in svlist:
                if any(x in line for x in sv[-2]):
                    if len(result)==0:
                        if "Sniffles2.INS" in sv[0]:
                            result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                        elif "Sniffles2.DEL" in sv[0]:
                            result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]+[sv[0].split(";END=")[1].split("\t")[0].split(";")[0]]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                    else:
                        found=0
                        for r in range(len(result)):
                            if result[r][0]==sv[0].split("\t")[2]+" "+str(sv[-1]):
                                result[r][2]=tuple(list(result[r][2])+line.split("RNAMES=")[1].split(";")[0].split(","))
                                found=1
                                break
                        if found==0:
                            if "Sniffles2.INS" in sv[0]:
                                result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                            elif "Sniffles2.DEL" in sv[0]:
                                result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]+[sv[0].split(";END=")[1].split("\t")[0].split(";")[0]]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                   
                    break
    with gzip.open(mosaicvcf,'r') as f:
        for line in f:
            line=line.decode("UTF-8")
            for sv in svlist:
                if any(x in line for x in sv[-2]):
                    if len(result)==0:
                        if "Sniffles2.INS" in sv[0]:
                            result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                        elif "Sniffles2.DEL" in sv[0]:
                            result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]+[sv[0].split(";END=")[1].split("\t")[0].split(";")[0]]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                    else:
                        found=0
                        for r in range(len(result)):
                            if result[r][0]==sv[0].split("\t")[2]+" "+str(sv[-1]):
                                result[r][2]=tuple(list(result[r][2])+line.split("RNAMES=")[1].split(";")[0].split(","))
                                found=1
                                break
                        if found==0:
                            if "Sniffles2.INS" in sv[0]:
                                result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                            elif "Sniffles2.DEL" in sv[0]:
                                result.append([sv[0].split("\t")[2]+" "+str(sv[-1]),tuple(sv[0].split("\t")[0:2]+[sv[0].split(";END=")[1].split("\t")[0].split(";")[0]]),tuple(line.split("RNAMES=")[1].split(";")[0].split(","))])
                   
                    break
    for r in range(len(result)):
        result[r]=tuple(result[r]+[[],[]])
    return tuple(result)

def read_overlaps_region(read_start, read_end, region_start, region_end):
    """Check if the read overlaps with the region."""
    return max(read_start, region_start) < min(read_end, region_end)

def get_reads(input_bam, readlist):
    import pysam
    lastchrom=''
    analysed=[]
    result=[]
    # Open input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as bam_in:
        # Iterate over each read in the input BAM file
        for read in bam_in:
            #found=0
            name=read.query_name
            chrom = bam_in.get_reference_name(read.reference_id)
            if lastchrom!=chrom:
                if analysed!=[]:
                    result+=list(analysed)
                analysed=[]
                lastchrom=chrom
                for i in readlist:

                    if i[1][0].lstrip("!")==chrom:
                        analysed.append(i)
                analysed=tuple(analysed)

            start = read.reference_start
            end = read.reference_end
            for i in range(len(analysed)):
                ##print(readlist[i])
                if name in set(analysed[i][2]):
                    analysed=list(analysed)
                    analysed[i]=list(analysed[i])
                    analysed[i][3].append(read.get_tag('RG')+":"+name)
                    analysed[i]=tuple(analysed[i])
                    analysed=tuple(analysed)
                #    found=1
                elif chrom==analysed[i][1][0].lstrip("!") and ((read_overlaps_region(start,end,int(analysed[i][1][1]),int(analysed[i][1][-1])) and "Sniffles2.DEL" in analysed[i][0]) or (read_overlaps_region(start,end,int(analysed[i][1][1]),int(analysed[i][1][-1])+1) and "Sniffles2.INS" in analysed[i][0])):
                    analysed=list(analysed)
                    analysed[i]=list(analysed[i])
                    analysed[i][4].append(read.get_tag('RG')+":"+name)
                    analysed[i]=tuple(analysed[i])
                    analysed=tuple(analysed)
                #    found=1
                #if found==1:
                #    break
    result+=list(analysed)
    return tuple(result)
            
               





def main():
    svs1=[]
    svs2=[]
    for sv in getfilteredSV(sys.argv[1]):
        if sv[0].startswith("!"): 
            svs1.append(tuple([sv[0], sv[1], sv[3]]))
            svs2.append(tuple([sv[0], sv[2], sv[3]]))
           
    print(svs1)
    print(svs2)

    tracedSV1=traceSV(sys.argv[2],sys.argv[3],tuple(svs1))
    tracedSV2=traceSV(sys.argv[4],sys.argv[5],tuple(svs2))  
    print(tracedSV1)
    print(tracedSV2)
    for i in get_reads(sys.argv[6], tracedSV1):
        print("\t".join(i[0].split(" ")) + "\t" + ",".join(sorted(i[-2])) + "\t" + ",".join(sorted(set(i[-2] + i[-1]))))
    for i in get_reads(sys.argv[7], tracedSV2):
        print("\t".join(i[0].split(" ")) + "\t" + ",".join(sorted(i[-2])) + "\t" + ",".join(sorted(set(i[-2] + i[-1]))))
if __name__ == "__main__":
    main()