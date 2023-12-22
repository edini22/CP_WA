#!/bin/bash
#SBATCH --ntasks=40
#SBATCH --time=00:10:00
#SBATCH --partition=cpar
#SBATCH --exclusive


threads=(1 2 4 8 12 15 18 20 22 25 28 30 33 35 38 40)


for run in {1..3}
do
    for nthreads in "${threads[@]}"
    do
        export OMP_NUM_THREADS=${nthreads}
        echo ${OMP_NUM_THREADS}
        time `./MDpar.exe <inputdata.txt > output.txt`
    done
done
