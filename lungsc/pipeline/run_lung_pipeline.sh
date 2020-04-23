#!/bin/sh
# Launch script that calls the Snakefile
# This script is only here for convenience in testing and to keep records of what was done

if [ $# -lt 1 ]; then
  echo "Arguments: {dry|cluster|selfbatch}"
  exit 1
fi

if [ $# -lt 2 ]; then
  # Skip run 29 for now as it needs separate merging
  dcs='DC20_DC21_DC23_DC25_DC28_DC30_DC34_DC35_DC36'
else
  dcs=$2
fi

case $1 in
  'dry')
    snakemake --dryrun --quiet --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition=normal,owners,quake --mem={params.mem} --time={params.time} --output={params.log}" --keep-target-files -j 100 -w 100 -k -p -r --rerun-incomplete --config "dcs=$dcs"
    ;;
  'cluster')
    snakemake --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition=normal,owners,quake --mem={params.mem} --time={params.time} --output={params.log}" --keep-target-files -j 100 -w 100 -k -p -r --rerun-incomplete --config "dcs=$dcs"
    ;;
  'unlock')
    snakemake --unlock
    ;;
  'unlock_cluster')
    snakemake --unlock && snakemake --cluster "sbatch --ntasks=1 --job-name={params.name} --cpus-per-task={threads} --partition=normal,owners,quake --mem={params.mem} --time={params.time} --output={params.log}" --keep-target-files -j 100 -w 100 -k -p -r --rerun-incomplete --config "dcs=$dcs"
    ;;
  'selfbatch')
    sbatch --ntasks=1 --job-name=LPipe --cpus-per-task=1 --partition=quake --mem=8000 --time='3-0:0:0' --output log/LPipe.log /oak/stanford/groups/quake/fzanini/postdoc/lung_development/lungsc/run_lung_pipeline.sh unlock_cluster $dcs
    ;;
  *)
    echo "Argument not understood: $1"
    exit 2
    ;;
esac
