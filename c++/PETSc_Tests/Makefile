# set compilers
CXX = g++
# set flags
CPPFLAGS = -Wall -O3 -std=c++11

# set paths for petsc and mpi libs
# PETSC_LIB = ${PETSC_DIR}/lib/petsc
PETSC_INC = ${PETSC_DIR}/include
PETSC_ARCH_INC = ${PETSC_DIR}/${PETSC_ARCH}/include
X11_INC = /opt/X11/include
MPI_LIB = /usr/local/Cellar/open-mpi/4.0.2/lib
GCC_LIB = /usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9
ARCH_DEBUG = /Users/tonyalberti/Desktop/petsc/arch-darwin-c-debug/lib
X11_LIB = /opt/X11/lib
LIBVENT_LIB = /usr/local/opt/libevent/lib
GCC_DARWIN = /usr/local/Cellar/gcc/9.2.0_3/lib/gcc/9/gcc/x86_64-apple-darwin19/9.2.0

# define library paths in addition to /usr/lib
LFLAGS = -L${MPI_LIB} -L${GCC_LIB} -L${GCC_DARWIN} -L${ARCH_DEBUG} -L${X11_LIB} -L${LIBVENT_LIB}

RPATH = -Wl,-rpath,${MPI_LIB} -Wl,-rpath,${GCC_LIB} -Wl,-rpath,${GCC_DARWIN} -Wl,-rpath,${ARCH_DEBUG} -Wl,-rpath,${X11_LIB} -Wl,-rpath,${LIBVENT_LIB}

# define header file directories
IFLAGS = -I${PETSC_INC} -I${PETSC_ARCH_INC} -I${X11_INC} -I./include

# define source file directories
SRCS = Main.cpp ./src/*

# define libraries to link into executable
LIBS = -lpetsc -llapack -lblas -lX11 -lc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lquadmath -lm

# define the default executable file name
EXEC = petsc_test.out

# define the main compilation scheme
default: main

# commands to compile and create executable
main:
	$(CXX) $(CPPFLAGS) $(LFLAGS) $(LIBS) $(IFLAGS) $(SRCS) -g3 -o $(EXEC)

# .PHONY ensures that clean will always work, even if there is a file called "clean"
.PHONY:
	clean cleanall

clean:
	rm $(EXEC)