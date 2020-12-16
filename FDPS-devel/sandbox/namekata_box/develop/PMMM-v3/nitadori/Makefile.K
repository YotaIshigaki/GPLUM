# (i) compiler, compile options
CXX = FCCpx -std=c++11 -Kfast -Ksimd=2 -Krestp=all 
LIBFFTW = -I/opt/klocal/fftw-3.3.6pl2-fujitsu/include -L/opt/klocal/fftw-3.3.6pl2-fujitsu/lib64  -lfftw3
HEADERS = vector3.h fmm.h ewald.h cell.h particle.h 

# (ii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.sh
WORKDIR = /data/ra000008/namekata/PMMM/nitadori/

pbc: pbc.cpp $(HEADERS)
	$(CXX) $< $(LIBFFTW) -o $@
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

obc: obc.cpp $(HEADERS)
	$(CXX) $< $(LIBFFTW) -o $@
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

ewald: ewald.cpp $(HEADERS)
	$(CXX) $< -o $@
	mkdir -p $(WORKDIR)
	mv $@ $(WORKDIR)
	cp $(JOB_FILE_NAME) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME) 

clean:
	rm -f obc pbc ewald ?diff*.dat

distclean:
	make clean
	rm -f *.txt
