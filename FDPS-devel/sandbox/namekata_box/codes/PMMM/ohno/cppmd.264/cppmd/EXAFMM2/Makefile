CXX	= mpicxx -ggdb3 -Wall -Wextra -Wshadow -Wuninitialized -O3 -fopenmp -I../include
#CXX	= mpicxx -xHOST -O3 -openmp -funroll-loops -finline-functions -fPIC -ansi-alias
#CXX	= FCCpx -Kfast,openmp

LFLAGS =

.cxx.o  :
	$(CXX) -c $? -o $@ $(LFLAGS)

a	: main.cxx
	$(CXX) $? $(LFLAGS) -DDO_P2P
	mpirun -np 16 ./a.out
#	mpirun -np 12 -machinefile ~/hosts ./a.out

b	: main.cxx
	$(CXX) $? $(LFLAGS) -DDO_P2P -DSerial
	./a.out
#	scp a.out mg:
#	ssh mg fpcoll -A -Icall,hwm ./a.out
#	ssh mg ./a.out

cleanall:
	find . -name "*.o" -o -name "*.out*" -o -name "*.dat" | xargs rm -rf
commit  :
	hg commit
	hg push
	hg pull -u
save    :
	make cleanall
	cd .. && tar zcvf cppmdfmm.tgz cppmdfmm
revert  :
	hg revert --all
	rm -rf `find . -name "*.orig"`
