module load CCEnv
module load StdEnv

module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6

mkdir fenics_17 && cd fenics_17
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
mv dolfin-2017.1.0 dolfin_17_1
chmod u+w ~/fenics_17/*/.git/objects/pack/*

pyvenv ~/packages/fenics_17
source ~/packages/fenics_17/bin/activate
cd fiat    && pip3 install . && cd -
cd instant && pip3 install . && cd -
cd dijitso && pip3 install . && cd -
cd ufl     && pip3 install . && cd -
cd ffc     && pip3 install . && cd -
pip3 install ply
pip3 install mpi4py
pip3 install lxml
pip3 install meshio
cd dolfin_17
mkdir build && cd build

cmake .. -DDOLFIN_SKIP_BUILD_TESTS=true -DEIGEN3_INCLUDE_DIR=$EBROOTEIGEN/include -DBOOST_ROOT="$SCINET_BOOST_ROOT" -DCMAKE_INSTALL_PREFIX=$HOME/software/dolfin_17 -DCMAKE_SKIP_RPATH=ON -DRT_LIBRARY=$EBROOTNIXPKGS/lib64/librt.so -DHDF5_C_LIBRARY_dl=$EBROOTNIXPKGS/lib64/libdl.so -DHDF5_C_LIBRARY_m=$EBROOTNIXPKGS/lib64/libm.so -DHDF5_C_LIBRARY_pthread=$EBROOTNIXPKGS/lib64/libpthread.so -DHDF5_C_LIBRARY_z=$EBROOTNIXPKGS/lib/libz.so
nice make -j 8 install && cd -
sed -i s'^export LD_LIBRARY_PATH=/lib^#export LD_LIBRARY_PATH=/lib^' ~/software/dolfin_17/share/dolfin/dolfin.conf





