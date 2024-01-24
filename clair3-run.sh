INPUT_DIR="/stornext/snfs4/next-gen/scratch/mike/rob"       
OUTPUT_DIR="/stornext/snfs5/next-gen/scratch/mike/rob/clair3_output$3"      # e.g. /home/user1/output (absolute path needed)
REF_DIR="/stornext/snfs5/next-gen/scratch/mike"
THREADS="24"            # e.g. 8
MODEL_NAME="r941_prom_hac_g360+g422"         # e.g. r941_prom_hac_g360+g422
SINGULARITY_LOCALCACHEDIR=/stornext/snfs5/next-gen/scratch/mike/temp
export SINGULARITY_LOCALCACHEDIR

mkdir -p "${OUTPUT_DIR}"
BIND_TMPDIR="${TMPDIR}"
unset TMPDIR
set -euo pipefail

singularity run --bind "${INPUT_DIR}:/mnt/input,${REF_DIR}:/mnt/reference,${OUTPUT_DIR}:/mnt/output" clair3_latest.sif /opt/bin/run_clair3.sh --bam_fn=/mnt/input/$1 --ref_fn=/mnt/reference/$2.fa --threads=${THREADS} --platform="ont" --model_path="/opt/models/${MODEL_NAME}" --output="/mnt/output"  --gvcf
