# Step1: for each shard, run stix individually

idxdir=$1
binddir=$2
gidx=$3
sidx=$4
invcf=$5
out=$6

sif_path=$(readlink -f sifs/stix-suite-1.0.0.sif)


cd ${idxdir}


singularity exec -B ${binddir} ${sif_path} stix \
    -i ${gidx} \
    -d ${sidx} \
    -f ${invcf} \
    -s 100 \
    -T 5 | tee ${out} 1>/dev/null
	
	
# Step2: merge the individual results 

python3 src/STIX_combine_individual_vcf_v2.py  -i shards*.vcf  -o merged.stix-anno.vcf;
