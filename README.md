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
     samtools view -h <sample.sorted.bam> -@4 | python MD_1SA5pp.py 1 5| samtools view -bS -@ 5 - > <sample.chmera.filtered.bam>
  3. Resulting filtered .bam was sorted and indexed. 
     samtools sort -@ 24 <sample.chmera.filtered.bam> > <sample.chmera.filtered.sorted.bam> 
     samtools index -@ 24 <sample.chmera.filtered.sorted.bam>  
  5. The HG002 samples are analysed with mosdepth to determine regions of >=5 coverage:   
     mosdepth --quantize 0:4:5:  <sample.chmera.filtered.sorted.bam>.quantized5.bed <sample.chmera.filtered.sorted.bam> 
  6. Resulting bed files were filtered for regions of >=5x coverage 
     grep "5:" <sample.bam>.quantized5.bed > <sample.bam>.quantized5.5.bed
  7. Coverage was summarized using getCoverage.py. 
  8. SNVs are called using Clair3 (clair3-run.sh).   
  9. SNV calls are split into substitutions and small IN/DELs (filterIndelSub.py).   
  10. Calls are filtered for those located in regions covered by 5 or more reads:   
     bedtools intersect -a <SNVcalls.vcf.gz> -b <sample.chmera.filtered.sorted.bam>.quantized5.5.bed   
  11. Calls supported by at least 3 reads are selected (SNVsup.py).    
  12. Filtering from points 7-9 was applied to GIAB HG002 benchmark .vcf files, using   
     <sample.chmera.filtered.sorted.bam>.quantized5.bed regions corresponding to each specific test   
     so that the regions within the tested sample and baseline were the same.   
  13. Filtered .vcf files were tabix indexed with bcftools:   
     bcftools index -t <SNVcalls.filtered.vcf.gz> 
  14. Calls are compared using RTGtools:   
      rtg vcfeval --baseline  <Baseline.SNVcalls.filtered.vcf.gz> --bed-regions <Baseline.SNVcalls.bed> -c <Sample.SNVcalls.filtered.vcf.gz> -o <output> -t <hg38.sdf/hg19.sdf> 
  15. SVs are called with Sniffles2:   
      sniffles2 --threads 24 --input <input_bam> --reference <hg38.fa/hg19.fa> --vcf <output_vcf> --snf <output_snf> 
  16. SV calls are filtered for those located in regions covered by 5 or more reads:   
      bedtools intersect -a <SVcalls.vcf.gz> -b <sample.bam>.quantized5.5.bed 
  17. SV calls are filtered for those supported by at least 3 reads and containing insertions or deletions (filterIndelSV.py), additionally split between insertions or deletions by using:   
      grep SVTYPE=<INS/DEL> <SVcalls.vcf.gz> 
  18. SV calls within the tested samples are benchmarked against the Giab truthset using Truvari:   
      truvari bench --passonly -b <baseline.vcf> --includebed <baseline.bed>  --pick multi -c <sample.vcf> -o	<output> 

Single-cell and corresponding bulk data pipeline:
  1. The data was obtained in .fastq format, aligned to GRCh38/Hg38 using minimap2 and subsequently converted to binary, sorted and indexed with samtools: 
     minimap2 -ax map-ont ../hg38.fa -t 24 <reads.fastq> > <sample.sam> 
     samtools view -bS <sample.sam> > <sample.bam> 
     samtools sort -@ 24 <sample.bam> > <sample.sorted.bam> 
     samtools index -@ 24 <sample.sorted.bam> 

  2. Single-cell samples were grouped by library preparation technique and tissue type. Groups of cell data were merged using Samtools: 
     samtools merge -r -o <merge_out.bam> <SC_1.sorted.bam> <SC_2.sorted.bam>...
  3. The aligned reads were filtered to remove chimeric reads: 
     samtools view -h <merge_out.bam> -@4 | python MD_1SA5pp.py 1 5| samtools view -bS -@ 5 - > <merge.chmera.filtered.bam> 
  4. The single cell assemblies and merges were analysed with mosdepth to determine regions of >=5 coverage:   
     mosdepth --quantize 0:4:5:  <sample.bam>.quantized5.bed <sample.bam> 
  5. Resulting bed files were filtered for regions of >=5x coverage 
     grep "5:" <sample.bam>.quantized5.bed > <sample.bam>.quantized5.5.bed
  6. Coverage was summarized using getCoverage.py. 
  7. SNVs are called using Clair3 (clair3-run.sh). 
  8. Resulting vcf were corrected to use uppercase letters in sequence and report correct sample names (fixglnexus.py). 
  9. SNV calls from groups of single cells and corresponding bulks were merged using Glnexus with jemalloc: 
       LD_PRELOAD=libjemalloc.so glnexus_cli_1.4.1 -c clair3.yml --mem-gbytes 32 --threads 12 <input.gvcf.gz> > <output.bcf> 
  10. SNV calls are split into substitutions and small IN/DELs (filterIndelSub.py).   
  11. Calls are filtered for those located in regions covered by 5 or more reads within single cells, both on single cell and corresponding bulk:   
     bedtools intersect -a <SNVcalls.vcf.gz> -b <sample.bam>.quantized5.5.bed 
  12. Calls are intersected with exons:   
     bedtools intersect -a <SNVcalls.5.vcf.gz> -b exons.bed.gz -u 
  13. Calls supported by at least 3 reads are selected and a summary of Bulk only (TN)), Single-cell only (FP) and shared calls is generated (TP) (filterSNV.py). Further, input in lines (91-101) in the script was changed for calls intersected with exons and analysed accordingly. 
      Resulting summary used to calculate F1 scores, fractions and absolute variant counts. 
  14. Similiarly to point 14, calls intersecting exons associated with MSA (gene list compiled from literature) supported by at least 3 reads are selected and a summary of Bulk only (TN)), Single-cell only (FP) and shared calls is generated (TP) (filterSNVgeneMSA.py). 
      Resulting summary used to calculate F1 scores, fractions and absolute variant counts. 
  15. SNV types (substitutions) were determined using snvtype.py script: 
      snvtype.py <REF_allele> <ALT_allele> input.vcf.gz 
  16. SNV locations were extracted from .VCF and converted into .bed format using snv2bed.py script.
  17. SNV In/Del ratio was determined using compareindels.py script.
  18. SNVs locations were extracted into bed format and shuffled 10,000 times for each experiment; the results were summarized and used in Z-test to determine results of permutation test (hg38-N.bed - list of unmapped regions represented as N in GRCh38, hg38len.bed - lengths of GRCh38 chromosomes): 

       for j in {1..10000}; do bedtools shuffle -i <input.SNV.bed> -chromFirst -excl hg38-N.bed -noOverlapping -g hg38len.bed | bgzip -c > shuffles/$(basename -s .bed $i).$j.bed.gz ; done 
     
  19. SVs were called with Sniffles2:   
      sniffles2 --threads 24 --input <input_bam> --reference <hg38.fa> --vcf <output_vcf> --output-rnames --snf <output_snf> 
      sniffles2 --threads 24 --input <input_bam> --reference <hg38.fa> --vcf <output_vcf> --output-rnames --noqc --mosaic 
  20. SV calls were merged with Sniffles2: 
      sniffles2 --threads 24 --input <SC+bulk.snf.list.tsv> --reference <hg38.fa> --vcf <output_vcf>
  21. Initial, unfiltered SV calls were split into duplications (cat svcalls.vcf | grep DUP) and inversions (cat svcalls.vcf | grep INV). Each of the files were split into SC+bulk, SC only and bulk only datasets using nulltest.py. Resulting files were used as input for sizeplotvcf.py to plot INV and DUP sizes and as
      input to medmax.sh script to get median and maximum sizes of these SVs.
  22. SV calls from merges are filtered for those located in regions covered by 5 or more reads:   
      bedtools intersect -a <SVcalls.vcf.gz> -b <sample.bam>.quantized5.5.bed 
  23. SV calls are filtered for those with PASS filter, supported by at least 3 reads in any of the merged single cells and containing insertions or deletions (filterSVmerge.py). 
  24. SV calls were split into SC-only,  bulk-only and SC+bulk (getFP.py). 
  25. SC+bulk and bulk-only 0/1 bulk genotypes were counted using awk and shell commands: 
      grep -v "#" <SC_bulk.vcf/Bulk_only.vcf>  | grep <INS/DEL> | awk '$<bulk_location> /0\/1/ {print $0}' | wc -l 
  26. Basic statistics of SVs were calculated using stats_SV.py script.  
  27. Reads that contain SVs are selected, traced back to their source SC .bam files and statistics on how many variants are present in how many single cells are produced (getReadname.py). 
  28. Filtered SV insertions are converted to .fasta, while deletions (+flanks) are converted to .bed (vcf2fasta.py) and subsequently extracted from reference genome (extractfromref.py). Insertion loci are also extracted from reference (extractfromref.py). 
  29. Fasta files of insertions and deletions (+flanks), as well as insertion loci extracted from reference are analysed with Repeatmasker: 
      RepeatMasker -dir $(basename -s .fasta <input_fasta>)_RMout -species human -s -e hmmer -pa 24 <input_fasta> 
  30. Statistics about Alu/Line insertions are calculated (alulineins.py)
  31. Numbers of INS into TE loci are calculated (insTE.py)
  32. Types of recombinational deletions are calculated (deletionsTE.py)
  33. Summary of TEs within deletions was calculated (checkAluLine.py)
  34. Illumina PicoPLEX bam files were tagged (RG) and CIPOS/CIEND tag was added to long-read deletions:
      samtools addreplacerg -r ID:<IO> -r SM:<ID> <INPUT.picoplex.merge.bam> -o <OUTPUT.picoplex.merge.rgtag.bam> -@ 7
      cat <long_read_del.vcf> | python fixcip[os.py > <long_read_del_cipos.vcf>
  35. SVtyper was used to genotype long-read deletions in PicoPLEX samples:
      svtyper -i <long_read_del_cipos.vcf> -B <picoplex.merge.rgtag.bam> -l <picoplex.merge.rgtag.bam>>.json > picoplex.merge.rgtag.GT.vcf> 
  36. Script filterill.py was used to get all genotyped deletions in each category, output was counted using wc -l.
  37. Locations of SVs  were shuffled 10000 times and intersected with exons for each of the shuffles, separately for INS and DEL.
      for j in {1..10000}; do bedtools shuffle -i <input.SV.bed> -chromFirst -excl hg38-N.bed -noOverlapping -g hg38len.bed | bgzip -c > shuffles/$(basename -s .bed $i).$j.bed.gz ; done
      bedtools intersect -a <SVcalls.X.vcf.gz> -b exons.bed.gz -u
  38. To summarize, countV.py script was used. Averages and standard deviations were calculated in MS Excel for each experiment and used to calculate z-score: (real_observed_value - mean_randomized) / stdev_randomized.
  39. SV vcf files were uploaded to online AnnotSV service (https://www.lbgi.fr/AnnotSV/), ran with default settings and resulting XML files were searched for genes related with MSA/neurodegeneration.
  40. Repeatmasker .out files for insertions were used to determine sizes of detected TEs using Onecodetofindthemall as per provided instructions (https://doua.prabi.fr/software/one-code-to-find-them-all). A filter80.py python script was used to filter output, getting counts of elements that spanned preset % TE reference:
      cat <Onecodetofindthemall_output_elem_sorted.csv> | python filter80.py <fraction_of_target_TE_ref> | grep <SINE/Alu or LINE/L1> | wc -l
  41. Locations of deletions  were shuffled 10000 times and intersected with reference SINE/Alu and LINE/L1 locations (extracted from https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.out.gz and selected with grep <SINE/Alu\|LINE/L1>, converted to bed with rmtobed.py script) for each of the shuffles, for DEL. Summarized with countTEs.py script.
      for j in {1..10000}; do bedtools shuffle -i <input.SV.bed> -chromFirst -excl hg38-N.bed -noOverlapping -g hg38len.bed | bgzip -c > shuffles/$(basename -s .bed $i).$j.bed.gz ; done
      bedtools intersect -a <SVcalls.X.vcf.gz> -b repeatmasker_hg38_out.bed.gz -u -wa -wb
42. Plots for sizes of INS and DEL were generated using sizeplotvcf.py:
    cat <source.vcf> | grep <INS/DEL> | python sizeplotvcf.py <Title_text> <output_file_name>
43. Plots for read size distributions were generated using chart_violin.py script.
44. Mosdepth coverage files, after filtering for +5x regions, were compared using bedtools intersect for MSA2 and control brains for PromethION and MiniON coverage in overlapping regions. Output was summarized using countcoveragebreadth.py and summary used to draw charts of coverage in MS Excel.
    bedtools intersect -a <MSA2_promethion.5.bed> -b ../<MSA2_minion.5.bed>  -wa -u > <MSA2.promethion.minion.intersect.regions.bed>
    cat <MSA2.promethion.minion.intersect.regions.bed> | python countcoveragebreadth.py
    
