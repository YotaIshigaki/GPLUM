CC = sw5cc.new
CXX = sw5CC
LINKER = mpiCC
all: main.out

CFLAGS = -O3 -msimd
LIBFLAGS = -hybrid #-lslave -allshare

main.out:master.o slave.o
	$(LINKER) $(LIBFLAGS) master.o slave.o -o main.out

master.o:master.cpp
	$(CXX) $(CFLAGS) -host -E master.cpp > master.e
	#	$(CXX) $(CFLAGS) -host -s master.cpp -o master.s
	$(CXX) $(CFLAGS) -host -c master.cpp -o master.o

slave.o:slave.c
	$(CC) $(CFLAGS) -slave -E slave.c > slave.e
	#$(CC) $(CFLAGS) -slave -s slave.c -o slave.s
	$(CC) $(CFLAGS) -slave -c slave.c -o slave.o

run:
	bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 4096 -host_stack 128 ./main.out
