# scONT
Single-cell Oxford Nanopore pipeline, supplement of Single cell long read whole genome sequencing to detect somatic mutations in human brain.

HG002 benchmark data links:

https://www.nist.gov/programs-projects/genome-bottle
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/CMRG_v1.00/

The HG002 data obtained from Oxford Nanopore Technologies was analysed for SNVs using Clair3 (clair3-run.sh). Calls were split into substitutions and small IN/DELs (filterIndelSub.py). Further, calls were filtered for those encompassing 5 or more
