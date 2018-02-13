#! /bin/bash
#This script is generated by rudzinski@pckr177
#on Mon Jan 12 17:02:31 CET 2015
#by the command "q2start -8p --walltime 00:15:00 mdrun_mpi"
#using q2start, version 0.5.5

# usage options
#$ -cwd
#$ -j y
#$ -N KINZIP
#$ -e $JOB_ID.log
#$ -o $JOB_ID_2.log
#$ -l h_rt=36:00:00
#$ -m bes
#$ -M rudzinski@mpip-mainz.mpg.de
#$ -S /bin/bash


#BEGIN COMPILER and MPI VARIABLES for gnu
#export LD_LIBRARY_PATH="/sw/linux/mpi/gcc/openmpi/lib:$LD_LIBRARY_PATH";
#export CC="gcc";
#export F77="gfortran";
#export CXX="g++";
#source /home/theorie/rudzinski/soft/anaconda/bin
#export python="/home/theorie/rudzinski/soft/anaconda/bin/python";
#export MPIEXEC="/sw/linux/mpi/gcc/openmpi/bin/mpiexec --prefix /sw/linux/mpi/gcc/openmpi/";
#LD_LIBRARY_PATH="/home/theorie/rudzinski/soft/anaconda/lib/:$LD_LIBRARY_PATH"
#LD_LIBRARY_PATH="/home/theorie/rudzinski/soft/anaconda/lib/python2.7:$LD_LIBRARY_PATH"
#LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
#LD_LIBRARY_PATH="/usr/local/lib64:$LD_LIBRARY_PATH"
export PATH="/sw/linux/mpi/gcc/openmpi/bin:$PATH";
export PATH="/data/isilon/rudzinski/soft_backup/anaconda-cluster/bin:$PATH"
#END  COMPILER and MPI VARIABLES for gnu

#BEGIN useful script variables
#walltime=00:15:00
#wallsecs=900
#wallhours=0
#ncpus=8
#END useful script variables

echo Hi, I am job $JOB_ID on $HOSTNAME in $PWD

echo Starting simulation in $PWD

# JFR - use my own script
#/home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/mpiexec -n 16 /home/theorie/rudzinski/soft/anaconda/envs/PyEmma-new/bin/python /home/theorie/rudzinski/soft/biased-MSM-2/biased_msm/biased-MSM.py $PWD
#/data/isilon/rudzinski/soft_backup/anaconda-cluster/bin/mpiexec -n 16 /data/isilon/rudzinski/soft_backup/anaconda-cluster/bin/python /home/theorie/rudzinski/soft/biased-MSM-2/biased_msm/biased-MSM.py $PWD

/data/isilon/rudzinski/soft_backup/anaconda-cluster/bin/python Get_LR_sampled_params_nocap.py

result=$?
[ $result -ne 0 ] && echo "$JOB_ID finished unhappy!"
