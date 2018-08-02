# This is a makefile.
# This makes a parallel simulation for different dielectric problem
# Use option -p in CC for profiling with gprof

# nanoHUB flags. 
nanoHUBCC = mpicxx -O3 -g -Wall -fopenmp -std=c++11 -Wunused-variable -Wunknown-pragmas
nanoHUBLFLAG = -lgsl -lgslcblas -lm -L${BOOST_LIBDIR} -lboost_filesystem -lboost_system -lboost_program_options -lboost_mpi -lboost_serialization
nanoHUBCFLAG = -I${BOOST_INCDIR}
nanoHUBOFLAG = -o
# BigRed2 flags. 
BigRed2CC = CC -O3 -g -Wall -fopenmp -std=c++11 -Wunused-variable -Wunknown-pragmas
BigRed2LFLAG = -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -lboost_program_options -lboost_mpi -lboost_serialization
BigRed2CFLAG = -c
BigRed2OFLAG = -o
# General purpose flags.
CC = mpicxx -O3 -g -Wall -fopenmp -std=c++11 -Wunused-variable -Wunknown-pragmas
LFLAG = -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -lboost_program_options -lboost_mpi -lboost_serialization
LFLAG = -lgsl -lgslcblas -lm -lboost_filesystem -lboost_system -L${BOOST_LIBDIR} -lboost_program_options -lboost_mpi -lboost_serialization
CFLAG = -c
OFLAG = -o


PROG = np_shape_lab
OBJ = main.o newmd.o newforces.o newenergies.o interface.o functions.o vertex.o edge.o face.o vector3d.o

all: $(PROG)

install: all
	@echo "compiling the np_electrostatics code in local computer"

nanoHUB-install:
	. /etc/environ.sh; use -e -r boost-1.62.0-mpich2-1.3-gnu-4.7.2; make CCF=nanoHUB all
	@echo "compiling the np_electrostatics code on Nanohub"

cluster-install:
	module swap PrgEnv-cray PrgEnv-gnu && module load boost/1.65.0 && module load gsl; make CCF=BigRed2 all
	@echo "compiling the np_electrostatics code on BIGRED 2"

$(PROG) : $(OBJ)
ifeq ($(CCF),BigRed2)	
	$(BigRed2CC) $(BigRed2OFLAG) $(PROG) $(OBJ) $(BigRed2LFLAG)
%.o : %.cpp
	$(BigRed2CC) -c $(BigRed2CFLAG) $< -o $@
else ifeq ($(CCF),nanoHUB)
	$(nanoHUBCC) $(nanoHUBOFLAG) $(PROG) $(OBJ) $(nanoHUBLFLAG)
%.o : %.cpp
	$(nanoHUBCC) -c $(nanoHUBCFLAG) $< -o $@
else
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LFLAG)
%.o : %.cpp
	$(CC) -c $(CFLAG) $< -o $@	
endif
