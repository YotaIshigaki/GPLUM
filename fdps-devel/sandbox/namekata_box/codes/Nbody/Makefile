#include ./mk_files/armagnac0.mk
include ./mk_files/K-computer.mk

#----------------------------------------------------------------------
#   Source and Object files
#----------------------------------------------------------------------
SRCS	= Physical_Constants.cpp \
	  Simulation_Units.cpp \
	  Calculation_Conditions.cpp \
	  Particle_Class.cpp \
	  Nbody_Objects_Class.cpp \
	  Kernel_Functions.cpp \
	  Time_Integration.cpp \
	  IO.cpp \
	  Initialize.cpp \
	  Initial_Conditions.cpp \
	  main.cpp
OBJS	= $(SRCS:.cpp=.o)
HDRS	= $(SRCS:.cpp=.h)

TARGET = nbody.x
RESTART	= nbody.x.rst

$(TARGET): Makefile
	tar --directory=$(LIBpDIR) -cvzf $(LIBCOMP) ./$(notdir $(LIBDIR))
	./create_preprocess_keywords.pl initial_run
	./check_undefined_macro.pl
	make obj
	mv $(LIBCOMP) $(WORKDIR)
	mv $(TARGET) $(WORKDIR)
	cp $(JOB_FILE) $(WORKDIR)
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME)

restart: Makefile
	tar --directory=$(LIBpDIR) -cvzf $(LIBCOMP) ./$(notdir $(LIBDIR))
	./create_preprocess_keywords.pl restart_run
	./check_undefined_macro.pl
	make obj_rst
	mv $(LIBCOMP) $(WORKDIR)
	mv $(RESTART) $(WORKDIR)
	cp $(RESTART_FILE) $(WORKDIR)
	cp ./scripts/make_restart_data.pl $(WORKDIR)
	cp ./scripts/arrange_calc_info.pl $(WORKDIR)

obj: Makefile $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET) $(LDFLAGS)

obj_rst: Makefile $(OBJS) $(HDRS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(RESTART) $(LDFLAGS)

%.o: %.cpp $(HDRS)
	$(CXX) -c $< $(CXXFLAGS)

clean:
	rm -f $(TARGET) $(RESTART) *.o fujitsu.tar.gz

