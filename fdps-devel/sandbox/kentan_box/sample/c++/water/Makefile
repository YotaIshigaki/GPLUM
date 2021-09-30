all:	rigid flex

rigid:
	make -C ./src/rigid/
	cp ./src/rigid/main.out rigid.out
flex:
	make -C ./src/flex/
	cp ./src/flex/main.out flex.out
gpu:
	make -C ./src/cuda/
	cp ./src/cuda/gpu.out .
reuse:
	make reuse -C ./src/cuda/
	cp ./src/cuda/reuse.out .
pezy:
	make -C ./src/pezy/
	cp ./src/pezy/pezy.out .
clean:
	make clean -C ./src/rigid/
	make clean -C ./src/flex/
	make clean -C ./src/cuda/
	make clean -C ./src/pezy/
	rm *.out
