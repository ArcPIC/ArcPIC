# Source this script to setup an LCG environment on lxplus7.cern.ch
# capable of building and running ArcPIC
# K. Sjobak, 11/02-2019

export LCGENV_PATH=/cvmfs/sft.cern.ch/lcg/releases

eval "`/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc8-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc8-opt gcc`"
eval "`/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc8-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc8-opt blas`"
eval "`/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc8-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc8-opt CMake`"
eval "`/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc8-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc8-opt hdf5`"
eval "`/cvmfs/sft.cern.ch/lcg/releases/LCG_94python3/lcgenv/1.3.6/x86_64-centos7-gcc8-opt/lcgenv -p LCG_94python3 x86_64-centos7-gcc8-opt GSL`"
