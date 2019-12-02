# set compilers
CXX = g++
# set flags
CPPFLAGS = -Wall -O3 -std=c++11

# set paths for petsc and mpi libs
PETSC_LIB = ${PETSC_DIR}/lib
PETSC_INC = ${PETSC_DIR}/include
MPI_LIB = /usr/local/Cellar/open-mpi/4.0.2/lib

# define library paths in addition to /usr/lib
LFLAGS = -L$(PETSC_LIB) -L$(MPI_LIB)

# define header file directories
IFLAGS = -I$(PETSC_INC) -I ./include

# define source file directories
SRCS = Main.cpp ./src/*

# define libraries to link into executable
LIBS = -lpetsc -lmpi

# define the default executable file name
EXEC = petsc_test.out

# define the main compilation scheme
default: main

# commands to compile and create executable
main:
	$(CXX) $(CPPFLAGS) $(LFLAGS) $(LIBS) $(IFLAGS) $(SRCS) -o $(EXEC)

# .PHONY ensures that clean will always work, even if there is a file called "clean"
.PHONY:
	clean cleanall

clean:
	rm $(EXEC)