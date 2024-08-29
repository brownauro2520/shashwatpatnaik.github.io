#!/bin/bash
# JOB HEADERS HERE
#SBATCH --job-name=liftcase1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=astrob@umich.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000m
#SBATCH --time=48:00:00
#SBATCH --account=aerosp623w23_class
#SBATCH --partition=standard
#SBATCH --export=ALL
#SBATCH --open-mode=truncate
#SBATCH --output=/home/astrob/liftcase1.log

folder_path="pictures"
rm -rf ${folder_path}/*
rm sols/COEFFICIENTSOUTPUTS.txt

# set the number of times to run the script
n=100

module load python/3.9.12

g++ -O3 1st_order.cpp -o main
g++ -O3 finres.cpp -o finres
g++ -O3 tapenade/adjoint.cpp -o adjoint -I/home/astrob/adolc_base/include -I/home/astrob/eigen -L/home/astrob/adolc_base/lib64 -ladolc

# loop n times, passing in a different number parameter each time
for (( i=0; i<=$n; i++ ))
do
    python3 python/plotgri.py $i
    
    python3 python/makemats.py $i coarse

    ./main $i load coarse
    
    python3 python/calc.py $i

    python3 python/inject.py $i

    ./finres $i

    ./main $i load fine

    python3 python/ojac.py
    
    ./adjoint

    python3 python/localref.py $i
    
done
