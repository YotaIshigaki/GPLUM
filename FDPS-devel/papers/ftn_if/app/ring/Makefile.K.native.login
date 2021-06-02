#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../../../../../../../../
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_FTN_MOD_DIR = $(FDPS_LOC)/src/fortran_interface/modules
FDPS_FTN_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_ftn_if.py

# (ii) Makefile used in computing nodes
MAKEFILE = Makefile.K.native.computing

# (iii) Job submission variables
QSUB = pjsub
JOB_FILE_NAME = job.K.native.computing.sh
WORKDIR = /data/ra000008/namekata/ftn_if/nbody/fortran


#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------
SRC_USER_DEFINED_TYPE = user_defined.F90
SRC_USER = f_main.F90
SRC_FDPS_MOD = $(wildcard $(FDPS_FTN_MOD_DIR)/*.F90)
SRC_FTN = $(SRC_FDPS_MOD) \
	  $(SRC_USER_DEFINED_TYPE) \
	  FDPS_module.F90 \
	  $(SRC_USER)

all:
	# Make a working directory
	mkdir -p $(WORKDIR)
	# Copy FDPS
	cd $(FDPS_LOC); tar -cvzf src.tar.gz src; mv -f src.tar.gz $(WORKDIR); cd -
	# Copy auto-generation script
	cp -Rf $(FDPS_LOC)/scripts   $(WORKDIR)
	# Copy user's files
	cp -f $(SRC_USER_DEFINED_TYPE) $(WORKDIR)
	cp -f $(SRC_USER)              $(WORKDIR)
	# Copy Makefile
	cp -f $(MAKEFILE)              $(WORKDIR)
	# Copy job file
	cp -f $(JOB_FILE_NAME)         $(WORKDIR)
	# Job submission
	cd $(WORKDIR); $(QSUB) $(JOB_FILE_NAME)
