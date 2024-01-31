import sys
import gzip
import subprocess

bgzip_path='' # path to bgzip program

files=tuple([""]) # file list to merge


for fi in files:
	foname=fi.split("/")[-2].rstrip("vcf.gz")+'fixed.vcf'
	with gzip.open(fi,'r') as f:
		for line in f:
			line=line.split("\t")
			if line[0]=="#CHROM":
				line[-1]=fi.split("/")[-2].rstrip("_clair3g")+"\n"
			elif line[0][0]!="#":
				line[3]=line[3].upper()
				line[4]=line[4].upper()
			with open(foname,'a') as fo:
				fo.write("\t".join(line))
	f.close()
	fo.close()
	result = subprocess.run([bgzip_path, foname])

