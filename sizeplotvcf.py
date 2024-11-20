# Script to create variant size distribution figure
# usage cat variants.vcf | python sizeplotvcf.py <Label> <output_name> 
def plotindel(aa):
    import sys
    lengths=[]
    for line in sys.stdin:
        if line.startswith("#"):
            continue
        else:
            lengths.append(int(line.split("\t")[7].split("SVLEN=")[1].split(";")[0].lstrip("-")))
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots()

    ran = ['100', '200', '300','400', '500', '1000','2000','5000','>5000']
    counts = []
    for i in ran:
        counts.append(0)
    for i in lengths:
        if i<=int(ran[0]):
            counts[0]+=1
        for j in range(1,len(ran[0:-1])):
            if i<=int(ran[j]) and i>int(ran[j-1]):
                counts[j]+=1
        if i>int(ran[-1].lstrip(">")):
            counts[-1]+=1
                
        
    plt.title(aa)
    ax.bar(ran, counts)
    ax.set_xlabel('Variant size [bp]')
    ax.set_ylabel('Variant count')

    plt.savefig(sys.argv[2]+".png",dpi=600)
    print(lengths)
def main():
    import sys
    plotindel(sys.argv[1])
if __name__=="__main__":
    main()