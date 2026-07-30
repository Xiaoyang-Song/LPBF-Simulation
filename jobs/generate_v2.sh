#!/bin/bash
#SBATCH --job-name=matlab_job_v2
#SBATCH --account=jhjin1
#SBATCH --partition=standard
#SBATCH --time=72:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8G
#SBATCH --output=/nfs/turbo/coe-sunwbgt/xysong/LPBF-Simulation/jobs/generation_v2.log

module load matlab

matlab -nodisplay -nosplash -r "try, run('simulation_v2/generate_trajectories_v2.m'); catch e, disp(getReport(e)); exit(1); end; exit(0);"
