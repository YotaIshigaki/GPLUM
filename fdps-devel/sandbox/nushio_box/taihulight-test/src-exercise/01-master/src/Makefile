0MCC = swcc
CC = sw5cc.new
.PHONY:all run

all: main.out

CFLAGS = -O3 -msimd
LIBFLAGS = -hybrid #-lslave -allshare

main.out:master.o
	$(CC) $(LIBFLAGS) master.o -o main.out

master.o:master.c
	$(CC) $(CFLAGS) -host -E master.c > master.e
	$(CC) $(CFLAGS) -host -s master.c -o master.s
	$(CC) $(CFLAGS) -host -c master.c -o master.o

run:main.out Makefile
	bsub -I -b -q q_sw_expr -n 1 -cgsp 64 -share_size 4096 -host_stack 128 ./main.out
