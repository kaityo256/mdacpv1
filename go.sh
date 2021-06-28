#!/bin/bash
#PBS -l nodes=1:ppn=20

cd $PBS_O_WORKDIR

time mpirun -np 20 ./mdacp L064_074.cfg
