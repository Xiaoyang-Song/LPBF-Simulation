#!/bin/bash
#SBATCH --job-name=matlab_job
#SBATCH --account=jhjin1
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --output=/nfs/turbo/coe-sunwbgt/xysong/LPBF-Simulation/jobs/generation.log

module load matlab

matlab -nodisplay -nosplash -r "try, run('generate_trajectories.m'); catch e, disp(getReport(e)); exit(1); end; exit(0);"
