#!/bin/bash -x
#PJM --name "Test.rst"
#PJM -o "00stdout-%j.log"
#PJM -j
#PJM --norestart
#PJM --rsc-list "node=16x16"
#PJM --rsc-list "elapse=24:00:00"
#PJM --rsc-list "rscgrp=small"
#PJM --mpi "shape=16x16, proc=2048"
#PJM -m b,e
#PJM --mail-list daisuke.namekata@riken.jp
#PJM --stg-transfiles all
#PJM --stgin "./fujitsu.tar.gz ./"
##PJM --stgin-dir "./tables ./tables"
##PJM --stgin "./make_restart_data.pl ./"
#PJM --stgin "./nbody.x.rst ./"
#PJM --stgin "./output/[a-zA-Z]*.dat ./output/"
#PJM --stgin "./job_history.dat ./ confirm=no"
#PJM --stgout-dir "./output ./output"
#PJM --stgout "./job_history.dat ./"
#PJM --stgout-dir "./profile_data ./profile_data"
#PJM -s
#PJM --spath "stat-%j.log"
#
# [Notes for staging]
#   1) By the instruction `--stgin "./output[a-zA-Z]*.dat ./output/"`,
#      all the files required to a restart job are transfered to 
#      computing nodes.
#   2) The file `job_history.dat` is created by make_restart_data.pl.
#
. /work/system/Env_base
export GSL_LIB=./fujitsu/gsl-2.1/lib
export QD_LIB=./fujitsu/lib
export LD_LIBRARY_PATH=${GSL_LIB}:${QD_LIB}:${LD_LIBRARY_PATH}
#
tar -xvzf fujitsu.tar.gz
#./arrange_calc_info.pl $PJM_JOBID
./make_restart_data.pl $PJM_JOBID
mpiexec -n 2048 ./nbody.x.rst -Wl,-T
#mpiexec -n 512 ./nbody.x.rst -debuglib -Wl,-T > 00stdout.log 2>&1
#fipp -C -Icall,hwm -d profile_data -m 40000 mpiexec -n 512 ./nbody.x.rst -Wl,-T
#fapp -C -Impi,hwm -Hevent=Cache,mode=sys,method=normal -d profile_data mpiexec -n 512 ./nbody.x.rst -Wl,-T

# [Notes]
#   1) In the K-computer, 1 node have 1 CPU which has 8 cores.
#   2) "-j" option merges STDERR to STDOUT.
#   3) "--norestart" option stops the system restarting this job. 
#   4) The group of resource is specified by the `--rsc-list "rscgrp="` option.
#      We have to select a resource group from the following:
#      --------------------------------------
#         Resource group  | Max. # of nodes 
#      --------------------------------------
#         small/interact  |   384        
#         micro           |   1152       
#         large           |   36864      
#         huge            |   82944      
#      --------------------------------------
#   5) If "-m" option is given, the system informs the user abount job status
#      by e-mails. The e-mail address should be specified by "--mail-list"
#      option. 
#   6) "-s" option activates the output of statistical information on job. 
#      We can specify the name of output file by "--spath" option.
#   7) With "-Wl,-T" option, we can perfrom IO with little endian.
#   8) The profiling data created by fipp can be converted into 
#      a text form by the command:
#      $ fipppx -Icpu,balance,call,hwm -A -d profile_data
#      On the other hand, that created by fapp can be converted by the command:
#      $ fapppx -Impi,hwm -ttext -o profile.txt -A -d profile_data
#
#
