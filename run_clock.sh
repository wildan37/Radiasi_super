#!/bin/bash

#PBS -l nodes=1:ppn=1,walltime=24:00:00
#PBS -q laser                         #name of the queue (default, laser, etc.)
#PBS -V
#PBS -m ae -M w.abdussalam@hzdr.de    #your email address
#PBS -N sr_clock                      #name of job 
#PBS -e sr_clock.$PBS_JOBID.err      #name of error file
#PBS -o sr_clock.$PBS_JOBID.out      #name of out file

. /etc/profile.modules
#module purge
module load gcc/4.8.2 
module load libxc/2.0.2
module load openmpi/1.8.0 
module load gsl/2.3 

HOMEDIR=/bigdata/hplsim/production/nu_urang/superradiance_atom_clock
TARGETDIR=/bigdata/hplsim/production/nu_urang/superradiance_atom_clock/out_clock
PARAMETERS=(0 0 0)

cd ${HOMEDIR}
EXEC=./a.out

echo ${HOMEDIR}
hostname

${EXEC} ${PARAMETERS[*]}
mv *dat ${TARGETDIR}
