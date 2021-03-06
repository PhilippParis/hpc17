#!/usr/bin/env bash

#SBATCH -p q_student
#SBATCH -N 31
#SBATCH --cpu-freq=HIGH
#SBATCH --job-name=31scatter

BENCHDIR=/home/students/e1325664/hpc17
OUTDIR=/home/students/e1325664/output31
MACHINE_FILE_NAME=nodelist

module load mpi/mvapich2-2.2

scontrol show hostnames ${SLURM_NODELIST} > ${MACHINE_FILE_NAME}

mpirun -np 496 -ppn 16 -bind-to core -machinefile ${MACHINE_FILE_NAME} ${BENCHDIR}/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > ${OUTDIR}/scatter_out1.dat
mpirun -np 496 -ppn 16 -bind-to core -machinefile ${MACHINE_FILE_NAME} ${BENCHDIR}/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > ${OUTDIR}/scatter_out2.dat
mpirun -np 496 -ppn 16 -bind-to core -machinefile ${MACHINE_FILE_NAME} ${BENCHDIR}/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > ${OUTDIR}/scatter_out3.dat
mpirun -np 496 -ppn 16 -bind-to core -machinefile ${MACHINE_FILE_NAME} ${BENCHDIR}/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > ${OUTDIR}/scatter_out4.dat
mpirun -np 496 -ppn 16 -bind-to core -machinefile ${MACHINE_FILE_NAME} ${BENCHDIR}/bin/mpibenchmark --input-file=/home/students/e1325664/hpc17/scripts/input_scatter.txt > ${OUTDIR}/scatter_out5.dat
