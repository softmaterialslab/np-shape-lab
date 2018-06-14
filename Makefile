# This is a makefile.
# This makes a parallel simulation for different dielectric problem
# Use option -p in CC for profiling with gprof

PROG = oncluster_simulate_membrane

OBJ = main.o newmd.o newforces.o newenergies.o interface.o functions.o vertex.o edge.o face.o vector3d.o

CC = g++ -O3 -g -Wall

LFLAG = -lgsl -lgslcblas -lboost_filesystem -lboost_system -lboost_program_options

CFLAG = -c

OFLAG = -o

$(PROG) : $(OBJ)
	$(CC) $(OFLAG) $(PROG) $(OBJ) $(LIBS) $(LFLAG)

main.o: utility.h interface.h vertex.h edge.h face.h control.h functions.h thermostat.h
newmd.o: utility.h interface.h thermostat.h functions.h newforces.h newenergies.h
newforces.o: newforces.h functions.h
newenergies.o: newenergies.h
interface.o: interface.h functions.h newenergies.h
functions.o: functions.h
vertex.o: vertex.h functions.h
edge.o: edge.h
face.o: face.h functions.h
vector3d.o: vector3d.h

clean:
	rm -f *.o

dataclean:
	rm -f outfiles/*

full:
	rm -f *.o
	rm -f outfiles/*
	make -j 10
	echo 900 0.01 3 30 30 | ./oncluster_simulate_membrane
	
	# d1 d2 qstr alpha c_s k_b k_s t_tot C T Q dt anneal_F off_F
