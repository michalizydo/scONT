# scONT
Single-cell Oxford Nanopore pipeline, supplement of Single cell long read whole genome sequencing to detect somatic mutations in human brain.

HG002 benchmark data links:

https://www.nist.gov/programs-projects/genome-bottle
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/


The HG002 data obtained from Oxford Nanopore Technologies. The analysis pipeline can be summarized as follows:
  1. The HG002 SC samples are analysed with mosdepth to determine regions of >=5 coverage:
       mosdepth --quantize 0:4:5:  <sample.bam>.quantized5.bed <sample.bam>        
  2. SNVs are called using Clair3 (clair3-run.sh).
  3. SNV calls are split into substitutions and small IN/DELs (filterIndelSub.py).
  4. Calls are filtered for those located in regions covered by 5 or more reads:
       bedtools intersect -a <SNVcalls.vcf.gz> -b <sample.bam>.quantized5.bed
  6. Calls supported by at least 3 reads are selected (SNVsup.py). 
