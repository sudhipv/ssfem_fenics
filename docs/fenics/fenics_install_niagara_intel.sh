module load NiaEnv
module load intel/2018.2  
module load intelmpi/2018.2
module load python/3.6.4-anaconda5.1.0
module load cmake/3.11.1
module load boost/1.66.0 
module load eigen/3.3.4
module load hdf5-mpi/1.8.20 
module load swig/3.0.12 
module load netcdf-mpi/4.6.1 
module load trilinos/12.12.1
module load gmp/6.1.2
module load mpfr/4.0.1 
module load petsc/3.8.4

mkdir fenics_intel && cd fenics_intel
wget https://bitbucket.org/fenics-project/fiat/downloads/fiat-2017.1.0.tar.gz
wget https://bitbucket.org/fenics-project/instant/downloads/instant-2017.1.0.tar.gz
wget https://bitbucket.org/fenics-project/dijitso/downloads/dijitso-2017.1.0.tar.gz
wget https://bitbucket.org/fenics-project/ufl/downloads/ufl-2017.1.0.tar.gz
wget https://bitbucket.org/fenics-project/ffc/downloads/ffc-2017.1.0.tar.gz
wget https://bitbucket.org/fenics-project/dolfin/downloads/dolfin-2017.1.0.tar.gz
tar xvfz fiat-2017.1.0.tar.gz
mv fiat-2017.1.0 fiat
tar xvfz instant-2017.1.0.tar.gz
mv instant-2017.1.0 instant
tar xvfz dijitso-2017.1.0.tar.gz
mv dijitso-2017.1.0 dijitso
tar xvfz ufl-2017.1.0.tar.gz
mv ufl-2017.1.0 ufl
tar xvfz ffc-2017.1.0.tar.gz
mv ffc-2017.1.0 ffc
tar xvfz dolfin-2017.1.0.tar.gz
mv dolfin-2017.1.0 dolfin_intel

python3.6 -m venv ~/packages/fenics_intel
source ~/packages/fenics_intel/bin/activate
cd fiat    && pip3 install . && cd -
cd instant && pip3 install . && cd -
cd dijitso && pip3 install . && cd -
cd ufl     && pip3 install . && cd -
cd ffc     && pip3 install . && cd -

pip3 install --upgrade sympy==1.1.1 
pip3 install scipy
pip3 install ply
pip3 install mpi4py
pip3 install lxml
pip3 install meshio
cd dolfin_intel
mkdir build && cd build

cmake .. -DDOLFIN_SKIP_BUILD_TESTS=true -DEIGEN3_INCLUDE_DIR="$SCINET_EIGEN_ROOT/include" -DBOOST_ROOT="$SCINET_BOOST_ROOT" -DCMAKE_INSTALL_PREFIX=$HOME/software/dolfin_intel -DCMAKE_SKIP_RPATH=ON -DHDF5_LIBRARIES=$SCINET_HDF5_MPI_ROOT/lib -DHDF5_INCLUDE_DIRS=$SCINET_HDF5_MPI_ROOT/include -DDOLFIN_ENABLE_TRILINOS:BOOL=ON -DTRILINOS_LIBRARIES=$SCINET_TRILINOS_ROOT/lib -DTRILINOS_INCLUDE_DIRS=$SCINET_TRILINOS_ROOT/include -DDOLFIN_ENABLE_CHOLMOD=OFF -DDOLFIN_ENABLE_PYTHON:BOOL=ON -DPYTHON_EXECUTABLE:FILEPATH=$(which python) 

 
#-DHDF5_C_LIBRARY_dl=$EBROOTNIXPKGS/lib64/libdl.so -DHDF5_C_LIBRARY_m=$EBROOTNIXPKGS/lib64/libm.so -DHDF5_C_LIBRARY_pthread=$EBROOTNIXPKGS/lib64/libpthread.so -DHDF5_C_LIBRARY_z=$EBROOTNIXPKGS/lib/libz.so 
#-DLIB_ifcore_pic=$EBROOTIFORT/lib/intel64/libifcore.so -DLIB_ipgo=$EBROOTIPP/lib/intel64/libipgo.a -DLIB_decimal=$EBROOTIPP/lib/intel64/libdecimal.a -DLIB_irc_s=$EBROOTIPP/lib/intel64/libirc_s.a 
nice make -j 8 install && cd -
sed -i s'^export LD_LIBRARY_PATH=/lib^#export LD_LIBRARY_PATH=/lib^' ~/software/dolfin_intel/share/dolfin/dolfin.conf







