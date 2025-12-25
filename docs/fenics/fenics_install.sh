
#### Assuming all the default modules are loaded
# module load nixpkgs/16.09 imkl/11.3.4.258 intel/2016.4 openmpi/2.1.1 StdEnv/2016.4

module load hdf5-mpi/1.8.18 boost eigen python/3.5 scipy-stack/2017b petsc/3.7.5 fftw-mpi/3.3.6
### Loading cmake to detect Boost
module load cmake
### Loading ipp - intel parallel processosing library
module load ipp/9.0.4

mkdir fenics_17_1 && cd fenics_17_1
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

pyvenv ~/fenics_17_1
source ~/fenics_17_1/bin/activate
cd fiat    && pip3 install . && cd -
cd instant && pip3 install . && cd -
cd dijitso && pip3 install . && cd -
cd ufl     && pip3 install . && cd -
cd ffc     && pip3 install . && cd -
pip3 install ply
pip3 install mpi4py
pip3 install meshio
pip3 install lxml
### Error corrected checking online answer
pip3 install --upgrade sympy==1.1.1

cd dolfin_17_1
mkdir build && cd build

#-DCMAKE_SKIP_INSTALL_RPATH=ON
cmake .. -DDOLFIN_SKIP_BUILD_TESTS=true -DEIGEN3_INCLUDE_DIR=$EBROOTEIGEN/include  -DCMAKE_INSTALL_PREFIX=$HOME/software/dolfin_17_1 -DCMAKE_SKIP_RPATH=ON -DRT_LIBRARY=$EBROOTNIXPKGS/lib64/librt.so -DHDF5_C_LIBRARY_dl=$EBROOTNIXPKGS/lib64/libdl.so -DHDF5_C_LIBRARY_m=$EBROOTNIXPKGS/lib64/libm.so -DHDF5_C_LIBRARY_pthread=$EBROOTNIXPKGS/lib64/libpthread.so -DHDF5_C_LIBRARY_z=$EBROOTNIXPKGS/lib/libz.so -DLIB_ifcore_pic=$EBROOTIFORT/lib/intel64/libifcore.so -DLIB_ipgo=$EBROOTIPP/lib/intel64/libipgo.a -DLIB_decimal=$EBROOTIPP/lib/intel64/libdecimal.a -DLIB_irc_s=$EBROOTIPP/lib/intel64/libirc_s.a

nice make -j 8 install && cd -
sed -i s'^export LD_LIBRARY_PATH=/lib^#export LD_LIBRARY_PATH=/lib^' ~/software/dolfin_17_1/share/dolfin/dolfin.conf










