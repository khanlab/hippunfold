#!/bin/bash
#SBATCH --account=rrg-lpalaniy
#SBATCH --ntasks=1
#SBATCH --gres=gpu:v100:8
#SBATCH --exclusive
#SBATCH --cpus-per-task=28
#SBATCH --mem=86000M
#SBATCH --time=24:00:00

module load arch/avx512 StdEnv/2018.3
nvidia-smi

singularity exec --nv hippocampal_autotop_latest.sif bash /src/resources/fineTune_UNet.sh mynewdataset/training mynewdataset/newCNNmodel
