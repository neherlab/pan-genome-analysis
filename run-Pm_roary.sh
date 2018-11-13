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

./panX.py -st 1 3 4 5 6 7 8 9 10 11 -rp ./data/Pm_roary/clustered_proteins_roary -fn ./data/Pm_roary -sl Pm_roary -t 64 > Pm_roary.log
## example for using divide-and-conquer algorithm on large datasets (#strains>50) (use parameters -dmdc and -dcs) (-sitr will not use treetime to compute mutations on branches of each gene tree)
#./panX.py -fn ./data/Pm_roary -sl Pm_roary -dmdc -sitr -t 32 > Pm_roary.log

## example for diverse datasets: (use parameter -cg)
#./panX.py -fn ./data/Pm_roary -sl Pm_roary -cg 0.7 -t 32 > Pm_roary.log

## example for using soft core_gene and core_gene strain constraint list: (use parameters -cg and -csf)
#./panX.py -fn ./data/Pm_roary -sl Pm_roary -cg 0.7 -csf ./data/Pm_roary/core_strain_list.txt -t 32 > Pm_roary.log
