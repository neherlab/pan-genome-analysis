#!/bin/sh
#
#  Reserve 32 CPUs for this job
#$ -pe smp 16
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#
#  Request it to run this for HH:MM:SS with ?G per core
#$ -l runtime=00:59:00,membycore=2G
#SBATCH --time=05:59:00
#
#  Run job from current working directory
#$ -cwd
#
export PYTHONPATH=/scicore/home/neher/neher/.local/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/scicore/home/neher/GROUP/bin:$PATH

#/etc/slurm/scripts/prolog.sh

ml purge
ml Python MAFFT FastTree RAxML
ml MCL/14.137-goolf-1.7.20
echo ${TMPDIR}

cp -r ./data/P_aeruginosa ${TMPDIR}/P_aeruginosa
./panX.py -fn ${TMPDIR}/P_aeruginosa -sl P_aeruginosa-RefSeq.txt -st 1 2 3 4 5 -t 16 --scratch ${TMPDIR} > P_aeruginosa1.log

ls ${TMPDIR}/P_aeruginosa/geneCluster > test.log

ml purge
ml Python MAFFT FastTree RAxML
./panX.py -fn ${TMPDIR}/P_aeruginosa -sl P_aeruginosa-RefSeq.txt -st 6 7 8 9 10 11 -t 16  --scratch ${TMPDIR}  > P_aeruginosa2.log

tar czf ${TMPDIR}/P_aeruginosa/vis.tar.gz ${TMPDIR}/P_aeruginosa/vis

cp ${TMPDIR}/P_aeruginosa/vis.tar.gz ./data/P_aeruginosa/

#/etc/slurm/scripts/epilog.sh
## example for using divide-and-conquer algorithm on large datasets (use parameters -dmdc and -dcs, maybe also -sitr )
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -dmdc 1 -dcs 50 -sitr 1 -t 32 > TestSet.log

## example for diverse datasets: (use parameter -cg)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -cg 0.7 -t 32 > TestSet.log

## example for using soft core_gene and core_gene strain constraint list: (use parameters -cg and -csf)
#./panX.py -fn ./data/TestSet -sl TestSet-RefSeq.txt -cg 0.7 -csf ./data/TestSet/core_strain_list.txt -t 32 > TestSet.log
