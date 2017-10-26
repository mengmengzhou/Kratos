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


#clear
reset

#you may want to decomment this the first time you compile
rm CMakeCache.txt


/opt/cmake-3.4.2/bin/cmake .. \
-DCMAKE_C_COMPILER=gcc \
-DCMAKE_CXX_COMPILER=g++ \
-DCMAKE_Fortran_COMPILER=gfortran \
-DCMAKE_CXX_FLAGS="${CMAKE_CXX_FLAGS} -lgfortran --std=c++11" \
-DCMAKE_C_FLAGS="${CMAKE_C_FLAGS} -lgfortran" \
-DCMAKE_Fortran_FLAGS="${CMAKE_Fortran_FLAGS} -lgfortran" \
-DCMAKE_BUILD_TYPE=RelWithDebInfo \
-DKRATOS_SD_REF_NUMBER=3 \
-DPYTHON_INCLUDE_DIR="/usr/include/python2.7" \
-DPYTHON_LIBRARY="/usr/lib/x86_64-linux-gnu/libpython2.7.so" \
-DPYTHON_EXECUTABLE="/usr/bin/python2.7" \
-DMPI_COMPILER="/opt/openmpi-1.10.1/bin/mpicxx" \
-DMPI_LIBRARY="MPI_LIBRARY-NOTFOUND" \
-DMETIS_APPLICATION=ON \
-DUSE_METIS_5=ON \
-DMETIS_ROOT_DIR="$HOME/opt/parmetis-4.0.3" \
-DSTRUCTURAL_APPLICATION=ON \
-DMKL_SOLVERS_APPLICATION=ON \
-DMKLSOLVER_INCLUDE_DIR="/opt/intel/mkl/include" \
-DMKLSOLVER_LIB_DIR="/opt/intel/mkl/lib/intel64" \
-DUSE_INTEL_GREATER_THAN_13=TRUE \


#decomment this to have it verbose
# make VERBOSE=1 -j4
make install -j`nproc`

