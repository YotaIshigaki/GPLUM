CC = sw5cc.new
CXX = sw5CC
LINKER = mpiCC

.PHONY:all run

all: main-gcc.out

CFLAGS = -O3 -msimd
LIBFLAGS = -hybrid #-lslave -allshare

main.out:master.o slave.o
	$(LINKER) $(LIBFLAGS) master.o -o main.out

main-gcc.out:master-gcc.o slave-gcc.o
	g++ master-gcc.o slave-gcc.o -o main.out


master.o:master.cpp
	$(CXX) $(CFLAGS) -host -E master.cpp > master.e
	$(CXX) $(CFLAGS) -host -s master.cpp -o master.s
	$(CXX) $(CFLAGS) -host -c master.cpp -o master.o

master-gcc.o:master.cpp
	g++ -std=c++11 -c \
  -O3 -funroll-loops -march=core-avx2 -ftree-vectorize -ftree-vectorizer-verbose=1 -mcmodel=large \
  -Wa,-adhln=master-gcc.s -g -fverbose-asm -masm=intel \
  master.cpp -o master-gcc.o


slave-gcc.o:slave.c
	gcc -std=c11 -c \
	-O0 \
  -Wa,-adhln=slave-gcc.s -fverbose-asm \
	slave.c -o slave-gcc.o

slave.o:slave.c
	$(CC) $(CFLAGS) -slave -E slave.c > slave.e
	$(CC) $(CFLAGS) -slave -s slave.c -o slave.s
	$(CC) $(CFLAGS) -slave -c slave.c -o slave.o

run:main.out Makefile
	bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 4096 -host_stack 128 ./main.out
