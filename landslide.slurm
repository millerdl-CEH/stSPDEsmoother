#!/bin/bash
# Slurm job for LandslideSTspde2026.R on the Crop Diversity HPC
# (https://help.cropdiversity.ac.uk).
#
# Submit from the project directory (the one containing
# LandslideSTspde2026.R, smooth.construct.spdeST.smooth.spec.R and data/),
# ideally under /mnt/shared/projects or /mnt/shared/scratch:
#
#   sbatch landslide.slurm
#
# Override the image location with:  IMAGE=/path/to/landslide.sif sbatch landslide.slurm
#
# Resources are a first guess. After the job, check what it used with
#   seff <jobid>      (or: sacct -j <jobid> --format=JobID,Elapsed,MaxRSS,ReqMem,State)
# and trim --mem/--time. Jobs that exceed --mem are killed.
#
# Queue: `long` (14-day limit, up to 384G per user). Use `medium` (24h) if
# the fit turns out to be quicker; use `himem` (512G-4T) if 384G isn't enough.

#SBATCH --job-name=landslide-dl
#SBATCH --partition=long
#SBATCH --cpus-per-task=4
# memory: the spdeST basis is ~28 GB dense with discrete=98^3 for that term,
# and model setup can hold 2-3 copies at once
#SBATCH --mem=192G
#SBATCH --time=7-00:00:00
#SBATCH --output=%x-%j.out
#SBATCH --open-mode=append
#SBATCH --export=ALL

set -euo pipefail

IMAGE="${IMAGE:-$APPS/landslide/landslide.sif}"

# the cluster provides `singularity`; newer installs also have `apptainer`
CONTAINER_CMD="$(command -v apptainer || command -v singularity)"

cd "$SLURM_SUBMIT_DIR"

echo "==== $(date) job $SLURM_JOB_ID on $(hostname)"
echo "restart count: ${SLURM_RESTART_COUNT:-0}"
echo "cpus: $SLURM_CPUS_PER_TASK  mem: ${SLURM_MEM_PER_NODE:-?}M"
echo "image: $IMAGE"
echo "workdir: $PWD"

if [[ ! -f "$IMAGE" ]]; then
  echo "Container image not found: $IMAGE" >&2
  exit 1
fi

# --cleanenv: don't leak host R/conda settings into the container
# binds: /mnt/shared (projects + scratch) and the job's local $TMPDIR;
#        $HOME is bound by default
# SLURM_CPUS_PER_TASK sets bam(nthreads=) in the R script
"$CONTAINER_CMD" exec \
  --cleanenv \
  --bind /mnt/shared:/mnt/shared \
  --bind "${TMPDIR:-/tmp}" \
  --env "TMPDIR=${TMPDIR:-/tmp},SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK}" \
  "$IMAGE" \
  Rscript --no-init-file LandslideSTspde2026.R

echo "==== $(date) finished"
