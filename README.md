# scONT
Single-cell Oxford Nanopore pipeline, supplement of Single cell long read whole genome sequencing to detect somatic mutations in human brain.

GIAB HG002 benchmark data links:

https://www.nist.gov/programs-projects/genome-bottle  
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/  
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/  

Exon locations were obtained using Gencode v44 database:

wget -qO- ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.annotation.gff3.gz | gunzip --stdout | awk '$3 == "exon"' | /hgsc_software/bedops/bedops-2.4.41/bin/convert2bed --input=gff --do-not-sort - | bgzip -c > exons.bed.gz


The HG002 data obtained from Oxford Nanopore Technologies. The analysis pipeline can be summarized as follows:  
  1. The data was obtained in .fastq format, aligned to GRCh38/Hg38 and GRCh37/Hg19 using minimap2 and subsequently converted to binary, sorted and indexed with samtools:
     minimap2 -ax map-ont ../hg38.fa -t 24 <reads.fastq> > <sample.sam>  
     samtools view -bS <sample.sam> > <sample.bam>  
     samtools sort -@ 24 <sample.bam> > <sample.sorted.bam>  
     samtools index -@ 24 <sample.sorted.bam>
  2. The aligned reads were filtered to remove chimeric reads:
     samtools view -h <sample.sorted.bam> -@4 | python MD_1SA5pp.py 1 5| samtools view -bS -@ 5 -) > <sample.chmera.filtered.bam>
  3. Resulting filtered .bam was sorted and indexed.
     samtools sort -@ 24 <sample.chmera.filtered.bam> > <sample.chmera.filtered.sorted.bam>
     samtools index -@ 24 <sample.chmera.filtered.sorted.bam> 
  5. The HG002 samples are analysed with mosdepth to determine regions of >=5 coverage:  
     mosdepth --quantize 0:4:5:  <sample.chmera.filtered.sorted.bam>.quantized5.bed <sample.chmera.filtered.sorted.bam>
  6. Resulting bed files were filtered for regions of >=5x coverage
     grep "5:" <sample.bam>.quantized5.bed > <sample.bam>.quantized5.5.bed
  7. SNVs are called using Clair3 (clair3-run.sh).  
  8. SNV calls are split into substitutions and small IN/DELs (filterIndelSub.py).  
  9. Calls are filtered for those located in regions covered by 5 or more reads:  
     bedtools intersect -a <SNVcalls.vcf.gz> -b <sample.chmera.filtered.sorted.bam>.quantized5.5.bed  
  10. Calls supported by at least 3 reads are selected (SNVsup.py).   
  11. Filtering from points 7-9 was applied to GIAB HG002 benchmark .vcf files, using  
     <sample.chmera.filtered.sorted.bam>.quantized5.bed regions corresponding to each specific test  
     so that the regions within the tested sample and baseline were the same.  
  12. Filtered .vcf files were tabix indexed with bcftools:  
     bcftools index -t <SNVcalls.filtered.vcf.gz>
  13. Calls are compared using RTGtools:  
      rtg vcfeval --baseline  <Baseline.SNVcalls.filtered.vcf.gz> --bed-regions <Baseline.SNVcalls.bed> -c <Sample.SNVcalls.filtered.vcf.gz> -o <output> -t <hg38.sdf/hg19.sdf>
  14. SVs are called with Sniffles2:  
      sniffles2 --threads 24 --input <input_bam> --reference <hg38.fa/hg19.fa> --vcf <output_vcf> --snf <output_snf>
  15. SV calls are filtered for those located in regions covered by 5 or more reads:  
      bedtools intersect -a <SVcalls.vcf.gz> -b <sample.bam>.quantized5.5.bed
  16. SV calls are filtered for those supported by at least 3 reads and containing insertions or deletions (filterIndelSV.py), additionally split between insertions or deletions by using:  
      grep SVTYPE=<INS/DEL> <SVcalls.vcf.gz>
  17. SV calls within the tested samples are benchmarked against the Giab truthset using Truvari:  
      truvari bench --passonly -b <baseline.vcf> --includebed <baseline.bed>  --pick multi -c <sample.vcf> -o	<output>

Single-cell and corresponding bulk data pipeline:
  1. The data was obtained in .fastq format, aligned to GRCh38/Hg38 using minimap2 and subsequently converted to binary, sorted and indexed with samtools:
     minimap2 -ax map-ont ../hg38.fa -t 24 <reads.fastq> > <sample.sam>
     samtools view -bS <sample.sam> > <sample.bam>
     samtools sort -@ 24 <sample.bam> > <sample.sorted.bam>
     samtools index -@ 24 <sample.sorted.bam>
  2. Single-cell samples were grouped by library preparation technique and tissue type. Groups of cell data were merged using Samtools:
     samtools merge -r -o <merge_out.bam> <cell_data1.bam> <cell_data2.bam>...
  3. The aligned reads were filtered to remove chimeric reads:
     samtools view -h <merge_out.bam> -@4 | python MD_1SA5pp.py 1 5| samtools view -bS -@ 5 -) > <merge.chmera.filtered.bam>
  4. The single cell assemblies and merges were analysed with mosdepth to determine regions of >=5 coverage:  
     mosdepth --quantize 0:4:5:  <sample.bam>.quantized5.bed <sample.bam>
  5. Resulting bed files were filtered for regions of >=5x coverage
     grep "5:" <sample.bam>.quantized5.bed > <sample.bam>.quantized5.5.bed
  7. SNVs are called using Clair3 (clair3-run.sh).
  8. Resulting vcf were corrected to use uppercase letters in sequence and report correct sample names (fixglnexus.py).
  9. SNV calls from groups of single cells and corresponding bulks were merged using Glnexus with jemalloc:
       LD_PRELOAD=libjemalloc.so glnexus_cli_1.4.1 -c clair3.yml --mem-gbytes 32 --threads 12 <input.gvcf.gz> > <output.bcf>
  10. SNV calls are split into substitutions and small IN/DELs (filterIndelSub.py).  
  11. Calls are filtered for those located in regions covered by 5 or more reads within single cells, both on single cell and corresponding bulk:  
     bedtools intersect -a <SNVcalls.vcf.gz> -b <sample.bam>.quantized5.5.bed
  13. Calls are intersected with exons:  
     bedtools intersect -a <SNVcalls.5.vcf.gz> -b exons.bed.gz -u
  14. Calls supported by at least 3 reads are selected and a summary of Bulk only (TN)), Single-cell only (FP) and shared calls is generated (TP) (filterSNV.py). Further, input in lines (91-101) in the script was changed for calls intersected with exons and analysed accordingly.
      Resulting summary used to calculate F1 scores, fractions and absolute variant counts.
  15. Similiarly to point 14, calls intersecting exons associated with MSA (gene list compiled from literature) supported by at least 3 reads are selected and a summary of Bulk only (TN)), Single-cell only (FP) and shared calls is generated (TP) (filterSNVgeneMSA.py).
      Resulting summary used to calculate F1 scores, fractions and absolute variant counts.
  16. SNV types (substitutions) were determined using snvtype.py script:
      snvtype.py <REF_allele> <ALT_allele> input.vcf.gz
  17. SNV locations were extracted from .VCF and converted into .bed format using snv2bed.py script.
  18. SNVs locations were extracted into bed format and shuffled 10,000 times for each experiment; the results were summarized and used in Z-test to determine results of permutation test (hg38-N.bed - list of unmapped regions represented as N in GRCh38, hg38len.bed - lengths of GRCh38 chromosomes):

       for j in {1..10000}; do bedtools shuffle -i <input.SNV.bed> -chromFirst -excl hg38-N.bed -noOverlapping -g hg38len.bed | bgzip -c > shuffles/$(basename -s .bed $i).$j.bed.gz ; done
     
  19. SVs are called with Sniffles2:  
      sniffles2 --threads 24 --input <input_bam> --reference <hg38.fa> --vcf <output_vcf> --snf <output_snf>
  20. SV calls are filtered for those located in regions covered by 5 or more reads:  
      bedtools intersect -a <SVcalls.vcf.gz> -b <sample.bam>.quantized5.5.bed
  21. SV calls are filtered for those with PASS filter, supported by at least 3 reads in any of the merged single cells and containing insertions or deletions (filterSVmerge.py).
  22. Basic statistics of SVs were calculated using stats_SV.py script. 
  23. Filtered SV insertions are converted to .fasta, while deletions are converted to .bed (vcf2fasta.py) and subsequently extracted from reference genome (extractfromref.py). 
