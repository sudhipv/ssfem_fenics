#!bin/bash
module load CCEnv
module load StdEnv
module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6

export CMAKE_PREFIX_PATH=/home/$USER/software/dolfin_17/share/dolfin/cmake/:$CMAKE_PREFIX_PATH

git clone https://bitbucket.org/fenics-project/mshr.git
cd mshr
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/software/mshr -DCMAKE_BUILD_TYPE=Release -DCMAKE_SKIP_RPATH=ON -DCMAKE_SYSTEM_PREFIX_PATH=$NIXUSER_PROFILE -DCMAKE_PREFIX_PATH=$NIXUSER_PROFILE:$CMAKE_PREFIX_PATH -DEIGEN3_INCLUDE_DIR=$EBROOTEIGEN/include
make
make install