module load hdf5-mpi/1.8.18 boost eigen python/3 python35-scipy-stack/2017a petsc/3.7.5 fftw-mpi/3.3.6

mkdir fenics && cd fenics
git clone https://bitbucket.org/fenics-project/fiat.git
git clone https://bitbucket.org/fenics-project/instant.git
git clone https://bitbucket.org/fenics-project/dijitso.git
git clone https://bitbucket.org/fenics-project/ufl.git
git clone https://bitbucket.org/fenics-project/ffc.git
git clone https://bitbucket.org/fenics-project/dolfin.git
chmod u+w ~/fenics/*/.git/objects/pack/*

pyvenv ~/fenics   ## niagara python3.6 -m venv ~/fenics
source ~/fenics/bin/activate
cd fiat    && pip3 install . && cd -
cd instant && pip3 install . && cd -
cd dijitso && pip3 install . && cd - 
cd ufl     && pip3 install . && cd - 
cd ffc     && pip3 install . && cd - 
pip3 install ply
pip3 install mpi4py
pip3 install meshio
pip3 install lxml
cd dolfin  
mkdir build && cd build  

cmake .. -DDOLFIN_SKIP_BUILD_TESTS=true -DEIGEN3_INCLUDE_DIR=$EBROOTEIGEN/include -DCMAKE_INSTALL_PREFIX=$HOME/software/dolfin -DCMAKE_SKIP_RPATH=ON -DRT_LIBRARY=$EBROOTNIXPKGS/lib64/librt.so -DHDF5_C_LIBRARY_dl=$EBROOTNIXPKGS/lib64/libdl.so -DHDF5_C_LIBRARY_m=$EBROOTNIXPKGS/lib64/libm.so -DHDF5_C_LIBRARY_pthread=$EBROOTNIXPKGS/lib64/libpthread.so -DHDF5_C_LIBRARY_z=$EBROOTNIXPKGS/lib/libz.so


nice make -j 8 install && cd -
sed -i s'^export LD_LIBRARY_PATH=/lib^#export LD_LIBRARY_PATH=/lib^' ~/software/dolfin/share/dolfin/dolfin.conf

## Niagara
cmake .. -DDOLFIN_SKIP_BUILD_TESTS=true -DEIGEN3_INCLUDE_DIR=$EBROOTEIGEN/include -DBOOST_ROOT="$SCINET_BOOST_ROOT" -DCMAKE_INSTALL_PREFIX=$HOME/software/dolfin -DCMAKE_SKIP_RPATH=ON -DRT_LIBRARY=$EBROOTNIXPKGS/lib64/librt.so -DHDF5_C_LIBRARY_dl=$EBROOTNIXPKGS/lib64/libdl.so -DHDF5_C_LIBRARY_m=$EBROOTNIXPKGS/lib64/libm.so -DHDF5_C_LIBRARY_pthread=$EBROOTNIXPKGS/lib64/libpthread.so -DHDF5_C_LIBRARY_z=$EBROOTNIXPKGS/lib/libz.so
