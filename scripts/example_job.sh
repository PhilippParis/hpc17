#!/usr/bin/env bash

#SBATCH -p q_student
#SBATCH -N 32
#SBATCH --cpu-freq=HIGH
#SBATCH --job-name=scatter_32x16
#SBATCH --ntasks-per-node=16

srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output/$
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output/$
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output/$
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output/$
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output/$







