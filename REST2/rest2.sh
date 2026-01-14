#!/bin/bash
#SBATCH --job-name=namd3TbxtRest2
#SBATCH --nodes=1                # 4 replicas
#SBATCH --gres=gpu:4              # 1 GPU per replica
#SBATCH --ntasks-per-node=1 # 1 process per node
#SBATCH --cpus-per-task=4 # 4 threads mapping to 4 cores per node
#SBATCH --partition=dept_gpu
#SBATCH --constraint="A4500"
#SBATCH --mail-user=rop174@pitt.edu
#SBATCH --mail-type=ALL
#SBATCH --time=144:00:00          # adjust based on estimated runtime


# -----------------------------
# Working directories
# -----------------------------
SCRDIR=/scr/${SLURM_JOB_ID}
mkdir -p $SCRDIR
cd $SCRDIR

# Copy the whole input directory
cp -r $SLURM_SUBMIT_DIR/* ./
cp $SLURM_SUBMIT_DIR/job0.conf ./
cp $SLURM_SUBMIT_DIR/stretchTbxtAf3g177d.conf ./
cd $SCRDIR/

# Output directory for results
OUTDIR=$SCRDIR/output
mkdir -p $OUTDIR

# Pre-create replica subdirectories for NAMD output
for i in $(seq 0 3); do
    mkdir -p $OUTDIR/$i
    cp equinpt.coor $OUTDIR/$i
    cp equinpt.xsc $OUTDIR/$i
    cp equinpt.vel $OUTDIR/$i
done

# Clean exit: move results back safely
trap "cp -r $OUTDIR /net/dali/home/mscbio/rop174/tbxt2035_mut/tbxtAf3g177dNamd3Imp 2>/dev/null" EXIT

NAMD_HOME="/net/dali/home/mscbio/rop174/tbxt2035/NAMD_3.0.2_Linux-x86_64-netlrts-smp-CUDA/"
PPN=$SLURM_CPUS_PER_TASK
$NAMD_HOME/charmrun $NAMD_HOME/namd3 ++local +p $PPN +replicas 4 +devicesperreplica 1 job0.conf +stdout /net/dali/home/mscbio/rop174/tbxt2035_mut/tbxtAf3g177dNamd3Imp/job0.%d.log


# -----------------------------
# Launch REST2 simulation
# -----------------------------
# mpirun -np 12 /net/dali/home/mscbio/rop174/tbxt2035/NAMD_3.0.2_Linux-x86_64-netlrts-smp-CUDA/namd3 +ppn 3 \
#     +replicas 4 job0.conf \
#     +stdout $OUTDIR/%d/job0.%d.log

