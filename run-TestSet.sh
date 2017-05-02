#!/bin/sh
#
#  Reserve 4 CPUs for this job
#$ -pe smp 4
#
#  Request it to run this for HH:MM:SS with ?G per core
#$ -l runtime=00:59:00,membycore=2G
#
#  Run job from current working directory
#$ -cwd
#
ml Python MAFFT FastTree RAxML MCL

./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -t 4 > TestSet.log
## example for using divide-and-conquer algorithm on large datasets (use parameters -dmdc and -dcs, maybe also -sitr )
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -dmdc 1 -dcs 50 -sitr 1 -t 32 > TestSet.log

## example for diverse datasets: (use parameter -cg)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -cg 0.7 -t 32 > TestSet.log

## example for using soft core_gene and core_gene strain constraint list: (use parameters -cg and -csf)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -cg 0.7 -csf ./data/TestSet/core_strain_list.txt -t 32 > TestSet.log
