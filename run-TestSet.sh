#!/bin/sh
#PBS -l nodes=1:ppn=15
# specify the time you expect the job to run hh:mm:ss
#PBS -l walltime=00:59:00
#specify the amount of memory needed
#PBS -l mem=4G
# output and error files
#PBS -o myout.o$PBS_JOBID
#PBS -e myout.e$PBS_JOBID

# load paths
source /home/neher/.bashrc

# move to current working directory
cd $PBS_O_WORKDIR

source activate panX

./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -t 15 > TestSet.log 
## example for using divide-and-conquer algorithm on large datasets (#strains>50) (use parameters -dmdc and -dcs) (-sitr will not use treetime to compute mutations on branches of each gene tree)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -dmdc -sitr -t 32 > TestSet.log

## example for diverse datasets: (use parameter -cg)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -cg 0.7 -t 32 > TestSet.log

## example for using soft core_gene and core_gene strain constraint list: (use parameters -cg and -csf)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -cg 0.7 -csf ./data/TestSet/core_strain_list.txt -t 32 > TestSet.log
