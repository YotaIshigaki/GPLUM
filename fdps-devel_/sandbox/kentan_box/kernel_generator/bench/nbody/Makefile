CC = g++
CFLAGS += -O2
INC = -I$(HOME)/FDPS/src
#CFLAGS		= -O0 -DUSE_ARM_SVE -Nlst=t

RUBY=/usr/local/opt/ruby/bin/ruby
MONOPOLE_PICG_FLAGS= --kernel-name CalcForceLongEpEp --epi-name EPI --epj-name EPJ --force-name Force --kernel-name CalcForceLongEpEp --strip-mining 64
QUADRUPOLE_PICG_FLAGS= --kernel-name CalcForceLongEpSp --epi-name EPI --epj-name SPJ --force-name Force --kernel-name CalcForceLongEpSp

#user defined option
#CFLAGS+= -Kfast,ocl,swp_strong,loop_nofusion,nounroll -K__control=0x2 -W0,-zie=Xdiv-fx100 -W0,-zie=Xsqrt-fx100 -W0,-zswpl=insert-move-into-prolog

main.out: main.cpp kernel.cpp
	$(CC) $(CFLAGS) $(INC) -DLOOP_DIVISION -DCHECK_FORCE_DETAIL -I./benchmark/include $< -o $@
kernel.cpp:	kernel.txt
	$(RUBY) -I ../../src ../../src/parserdriver.rb -i kernel.txt $(MONOPOLE_PICG_FLAGS) > kernel.cpp
	$(RUBY) -I ../../src ../../src/parserdriver.rb -i quadrupole.txt $(QUADRUPOLE_PICG_FLAGS) >> kernel.cpp
quadrupole.cpp:	quadrupole.txt
	$(RUBY) -I ../../src ../../src/parserdriver.rb -i quadrupole.txt $(QUADRUPOLE_PICG_FLAGS) > quadrupole.cpp
profile: main.cpp Makefile
	$(CC) $(CFLAGS) $(INC) -Nlinkprof -DLOOP_DIVISION main.cpp -o prof.out
clean:
	rm -rf *.out m5out *.lst *.o
