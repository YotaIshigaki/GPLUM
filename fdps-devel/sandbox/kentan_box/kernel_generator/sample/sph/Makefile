PS_PATH = ../../../src/
INC = -I$(PS_PATH)

CC = time g++
#CC = time mpicxx
CFLAGS = -std=c++11 -O3
#CFLAGS += -Wall
#CFLAGS += -ffast-math
#CFLAGS += -funroll-loops
#CFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

RUBY=/usr/local/opt/ruby/bin/ruby

all:hydro_force.h density_pressure.h

PG_BUILD =
PG_CLEAN = 

hydro_force.h:hydro_force.txt
	$(RUBY) -I ../../src ../../src/parserdriver.rb -i $< > $@
density_pressure.h:density_pressure.txt
	$(RUBY) -I ../../src ../../src/parserdriver.rb -i $< > $@

clean:
	rm -f *.o *~

distclean: clean
	rm -f nbody.out
	rm -rf result
	rm kernel.h

test: 
# This command is only for FDPS developers.
	./test.py

# fdps-autotest-run (DO NOT CHANGE THIS LINE)
