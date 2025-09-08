#!/bin/bash  
#PBS -N %TASK%
#PBS -r n
#PBS -j oe
#PBS -o %LOGDIR%/%TASK%.out
#PBS -S /bin/bash
#PBS -l select=1:ncpus=1:mpiprocs=1:mem=20GB
#PBS -l job_priority=regular
#PBS -q develop
#PBS -A %PROJECT%
#PBS -l walltime=02:00:00
#PBS -V
