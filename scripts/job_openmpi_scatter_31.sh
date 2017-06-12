#!/usr/bin/env bash

#SBATCH -p q_student
#SBATCH -N 31
#SBATCH --cpu-freq=HIGH
#SBATCH --job-name=scatter_32x16
#SBATCH --ntasks-per-node=16

srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output31/scatter_out1.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output31/scatter_out2.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output31/scatter_out3.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output31/scatter_out4.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > /home/students/e1325664/output31/scatter_out5.dat






