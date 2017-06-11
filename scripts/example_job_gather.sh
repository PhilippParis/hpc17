#!/usr/bin/env bash

#SBATCH -p q_student
#SBATCH -N 32
#SBATCH --cpu-freq=HIGH
#SBATCH --job-name=gather_32x16
#SBATCH --ntasks-per-node=16

srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_gather.txt > /home/students/e1325664/output/gather_out1.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_gather.txt > /home/students/e1325664/output/gather_out2.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_gather.txt > /home/students/e1325664/output/gather_out3.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_gather.txt > /home/students/e1325664/output/gather_out4.dat
srun ~/hpc17/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_gather.txt > /home/students/e1325664/output/gather_out5.dat






