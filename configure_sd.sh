#!/bin/sh

# this is an example file...please DO NOT MODIFY if you don't want to do it for everyone
#to use it, copy it to another file and run it

# additional compiler flags could be added customizing the corresponding var, for example
# -DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -msse3 ". Note that the "defaults are already correctly coded"
#so we should ass here only machine specific stuff

#an effort is made to autodetect all of the libraries needed HOWEVER:
#METIS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
#TRILINOS_APPLICATION needs the var PARMETIS_ROOT_DIR to be specified by the user (not needed if the app is set to OFF)
#MKL_SOLVERS_APPLICATION needs the var MKLSOLVER_INCLUDE_DIR and MKLSOLVER_LIB_DIR to be specified by the user (not needed if the app is set to OFF)
#note that the MKLSOLVER_LIB_DIR should include /lib/em64t. This is needed as intel is changing location of mkl at every update of the compiler!!

#the user should also note that the symbol "\" marks that the command continues on the next line. IT SHOULD ONLY BE FOLLOWED
#BY the "ENTER" and NOT by any space!!

reset

#you may want to decomment this the first time you compile
rm CMakeCache.txt
#rm *.cmake
#rm -rf CMakeFiles\

#cd ..
#rm applications
#ln -s applications_sd applications
#cd -

/opt/cmake-3.4.2/bin/cmake .. \
-DCMAKE_CXX_COMPILER=g++ \
-DCMAKE_C_COMPILER=gcc \
-DCMAKE_INSTALL_RPATH="${HOME}/kratos_bcn3/libs" \
-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=TRUE \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -std=c++11 " \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} " \
-DCMAKE_BUILD_TYPE=Release \
-DKRATOS_SD_REF_NUMBER=3 \
-DMPI_C_COMPILER="$MPI_ROOT/bin/mpicc" \
-DMPI_CXX_COMPILER="$MPI_ROOT/bin/mpicxx" \
-DPYTHON_INCLUDE_DIR="/usr/include/python2.7" \
-DPYTHON_LIBRARY="/usr/lib/x86_64-linux-gnu/libpython2.7.so" \
-DPYTHON_EXECUTABLE="/usr/bin/python2.7" \
-DSTRUCTURAL_APPLICATION=ON \
-DEXTERNAL_SOLVERS_APPLICATION=ON \
-DMKL_SOLVERS_APPLICATION=ON \
-DMKLSOLVER_INCLUDE_DIR="/opt/intel/mkl/include" \
-DMKLSOLVER_LIB_DIR="/opt/intel/mkl/lib/intel64" \
-DUSE_INTEL_GREATER_THAN_13=TRUE \
-DEKATE_AUXILIARY_APPLICATION=OFF \
-DEXTERNAL_CONSTITUTIVE_LAWS_APPLICATION=OFF \
-DFREEZING_SOIL_APPLICATION=OFF \
-DMORTAR_APPLICATION=OFF \
-DISOGEOMETRIC_APPLICATION=ON \
-DUSE_METIS_5=ON \
-DMETIS_ROOT_DIR="${HOME}/opt/parmetis-4.0.3" \
-DMETIS_APPLICATION=ON \
-DPETSC_DIR="${HOME}/opt/petsc-dev" \
-DTRILINOS_ROOT="${HOME}/opt/trilinos-12.10.1" \
-DUSE_TRILINOS_EPETRA=ON \
-DUSE_TRILINOS_TPETRA=ON \
-DDISTRIBUTED_BUILDERS_APPLICATION=ON \

make install -j8

