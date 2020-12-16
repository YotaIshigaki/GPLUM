# Target objects
SRC_CPP_SMPL = src_cpp_smpl
SRC_FTN_SMPL = src_ftn_smpl
SRC_CPP_NBODY_KERNEL = src_cpp_nbody_kernel
SRC_FTN_NBODY_KERNEL = src_ftn_nbody_kernel

# Shell commands
LATEX = platex -shell-escape
PDF2EPS = pdf2eps
DVIPDFMX = dvipdfmx


all:
	make cpp_smpl
	make ftn_smpl
	make cpp_nbody_kernel
	make ftn_nbody_kernel

cpp_smpl:
	$(LATEX) $(SRC_CPP_SMPL)
	$(DVIPDFMX) $(SRC_CPP_SMPL)
	$(PDF2EPS) 1 $(SRC_CPP_SMPL)

ftn_smpl:
	$(LATEX) $(SRC_FTN_SMPL)
	$(DVIPDFMX) $(SRC_FTN_SMPL)
	$(PDF2EPS) 1 $(SRC_FTN_SMPL)
	
cpp_nbody_kernel:
	$(LATEX) $(SRC_CPP_NBODY_KERNEL)
	$(DVIPDFMX) $(SRC_CPP_NBODY_KERNEL)
	$(PDF2EPS) 1 $(SRC_CPP_NBODY_KERNEL)

ftn_nbody_kernel:
	$(LATEX) $(SRC_FTN_NBODY_KERNEL)
	$(DVIPDFMX) $(SRC_FTN_NBODY_KERNEL)
	$(PDF2EPS) 1 $(SRC_FTN_NBODY_KERNEL)

clean:
	rm -f *.aux *.dvi *.log *.toc *.out

distclean:
	make clean
	rm -f *.pdf
