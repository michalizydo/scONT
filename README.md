# scONT
Single-cell Oxford Nanopore pipeline, supplement of Single cell long read whole genome sequencing to detect somatic mutations in human brain.

GIAB HG002 benchmark data links:

https://www.nist.gov/programs-projects/genome-bottle  
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/  
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/  


The HG002 data obtained from Oxford Nanopore Technologies. The analysis pipeline can be summarized as follows:  
  1. The data was obtained in .fastq format and aligned to GRCh38/Hg38 and GRCh37/Hg19 using minimap2 and subsequently converted to binary, sorted and indexed with samtools:
     minimap2 -ax map-ont ../hg38.fa -t 24 <reads.fastq> > <sample.sam>
     samtools view -bS <sample.sam> > <sample.bam>
     samtools sort -@ 24 <sample.bam> > <sample.sorted.bam>
     samtools index -@ 24 <sample.sorted.bam>
  3. The HG002 SC samples are analysed with mosdepth to determine regions of >=5 coverage:  
     mosdepth --quantize 0:4:5:  <sample.bam>.quantized5.bed <sample.bam>    
  4. SNVs are called using Clair3 (clair3-run.sh).  
  5. SNV calls are split into substitutions and small IN/DELs (filterIndelSub.py).  
  6. Calls are filtered for those located in regions covered by 5 or more reads:  
       bedtools intersect -a <SNVcalls.vcf.gz> -b <sample.bam>.quantized5.bed  
  7. Calls supported by at least 3 reads are selected (SNVsup.py).   
  8. Filtering from points 3-5 was applied to GIAB HG002 benchmark .vcf files, using
     <sample.bam>.quantized5.bed regions corresponding to each specific test
     so that the regions within the tested sample and baseline were the same.
  9. Filtered .vcf files were tabix indexed with bcftools:
       bcftools index -t <SNVcalls.filtered.vcf.gz>
  10. Calls are compared using RTGtools:
       rtg vcfeval --baseline  <Baseline.SNVcalls.filtered.vcf.gz> --bed-regions <Baseline.SNVcalls.bed> -c <Sample.SNVcalls.filtered.vcf.gz> -o <output> -t <hg38.sdf/hg19.sdf>
  11. SVs are called with Sniffles2:
       sniffles2 --threads 24 --input <input_bam> --reference <hg38.fa/hg19.fa> --vcf <output_vcf> --snf <output_snf>
  12. SV calls are filtered for those located in regions covered by 5 or more reads:  
       bedtools intersect -a <SVcalls.vcf.gz> -b <sample.bam>.quantized5.bed
  13. SV calls are filtered for those supported by at least 3 reads and containing insertions or deletions (filterIndelSV.py), additionally split between insertions or deletions by using:
      grep SVTYPE=<INS/DEL> <SVcalls.vcf.gz>
  14. SV calls within the tested samples are benchmarked against the Giab truthset using Truvari:
      

Single-cell and corresponding bulk data 
