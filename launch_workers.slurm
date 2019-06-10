#!/bin/bash
#
# Example shell script for running job that runs off the Wilmerlab subjobserver.
# $Revision: 1.0 $
# $Date:  2016-03-21 $
# $Author: paulboone $

#SBATCH --job-name=BAZGAM_test2
#SBATCH --output=BAZGAM.out
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=3
#SBATCH --cpus-per-task=1
#SBATCH --time=0-00:20:00

# accepts a parameter stay_alive if you don't want the worker to exit immediately after all jobs
# are complete
# use like `qsub -v stay_alive=1`

echo JOB_ID: $SLURM_JOBID JOB_NAME: $SLURM_JOBNAME HOSTNAME: $SLURM_O_HOST
echo start_time: `date`

# dependencies
module purge
module load python/3.5.1
source ~/venv/mof_sensor_array_ads/bin/activate

#cd /ihome/cwilmer/jag227/RASPA-2.0/jenna/Ads/test_h2p
cp -pR $SLURM_SUBMIT_DIR/* $SLURM_SCRATCH
cd $SLURM_SCRATCH

run_on_exit(){
  cp -pR $SLURM_SCRATCH/output_* $SLURM_SUBMIT_DIR/
  cp -pR $SLURM_SCRATCH/logs* $SLURM_SUBMIT_DIR/
  cp -pR $SLURM_SCRATCH/*.out $SLURM_SUBMIT_DIR/
}
trap run_on_exit EXIT

sjs_launch_workers.sh $SLURM_NTASKS_PER_NODE $stay_alive
