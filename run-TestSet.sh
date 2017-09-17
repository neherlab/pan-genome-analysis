#!/bin/sh
#
#  Reservation desired
#$ -R y
#
#  Reserve 8 CPUs for this job
#$ -pe parallel 64
#
#  Request 8G of RAM
#$ -l h_vmem=1G
#
#  Request it to run this long HH:MM:SS
#$ -l h_rt=00:59:00
#
#  Use /bin/bash to execute this script
#$ -S /bin/bash
#
#  Run job from current working directory
#$ -cwd
#
#  Send email when the job begins, ends, aborts, or is suspended
#$ -m beas

./panX.py -fn ./data/TestSet -sl TestSet -t 64 > TestSet.log
## example for using divide-and-conquer algorithm on large datasets (#strains>50) (use parameters -dmdc and -dcs) (-sitr will not use treetime to compute mutations on branches of each gene tree)
#./panX.py -fn ./data/TestSet -sl TestSet -dmdc -sitr -t 32 > TestSet.log

## example for diverse datasets: (use parameter -cg)
#./panX.py -fn ./data/TestSet -sl TestSet -cg 0.7 -t 32 > TestSet.log

## example for using soft core_gene and core_gene strain constraint list: (use parameters -cg and -csf)
#./panX.py -fn ./data/TestSet -sl TestSet -cg 0.7 -csf ./data/TestSet/core_strain_list.txt -t 32 > TestSet.log
