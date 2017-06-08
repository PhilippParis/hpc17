#!/usr/bin/env bash

#SBATCH -p q_student
#SBATCH -N 32
#SBATCH --cpu-freq=HIGH
#SBATCH --job-name=job_32x16
#SBATCH --ntasks-per-node=16

srun /home/hunold/exp/sprojects/reprompi-1.0.1-src/bin/mpibenchmark --input-file=/home/hunold/exp/sprojects/input_gather.txt > /home/hunold/exp/sprojects/output/gather_out1.dat
srun /home/hunold/exp/sprojects/reprompi-1.0.1-src/bin/mpibenchmark --input-file=/home/hunold/exp/sprojects/input_gather.txt > /home/hunold/exp/sprojects/output/gather_out2.dat
srun /home/hunold/exp/sprojects/reprompi-1.0.1-src/bin/mpibenchmark --input-file=/home/hunold/exp/sprojects/input_gather.txt > /home/hunold/exp/sprojects/output/gather_out3.dat
srun /home/hunold/exp/sprojects/reprompi-1.0.1-src/bin/mpibenchmark --input-file=/home/hunold/exp/sprojects/input_gather.txt > /home/hunold/exp/sprojects/output/gather_out4.dat
srun /home/hunold/exp/sprojects/reprompi-1.0.1-src/bin/mpibenchmark --input-file=/home/hunold/exp/sprojects/input_gather.txt > /home/hunold/exp/sprojects/output/gather_out5.dat

