CCC = g++

CFLAGS = -O3 -ffast-math -funroll-loops #-Wall

OBJ = main.o
ETC = Makefile

run: $(OBJ) $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -o run $(OBJ)

main.o: main.cc $(HDR) $(ETC)
	$(CCC) $(PS_PATH) $(CFLAGS) -c main.cc

clean:
	rm -f *~
distclean:
	make clean
	rm -f *.o
allclean:
	make distclean
	rm -f run
