#
# Makefile for Cray XC30 and XC50.
#

#For XC30
#CC = /home/saitoutk/bin/bin/ccache cc
CC = cc
LINKER = cc

#FOR GNU
OPTIMIZE = -std=c99
OPTIMIZE += -O3 -finline-functions -ffast-math -funroll-loops -frerun-cse-after-loop -fexpensive-optimizations -fomit-frame-pointer -fprefetch-loop-arrays
OPTIMIZE += -Winline -Wundef -Wmissing-noreturn -Wstrict-prototypes
OPTIMIZE += -mtune=native -march=native
#OPTIMIZE += -g

#FOR Intel
#OPTIMIZE = -std=c99  -O3 -fast
#OPTIMIZE = -std=c99  -O3 -ipo -no-prec-div -static -fp-model fast=2 -xHost 
#OPTIMIZE = -std=c99 -g
#OPTIMIZE = -std=c99 -O3 
#OPTIMIZE += -ipo 
#OPTIMIZE += -Ofast -fp-model fast=2

#FOR Cray
#OPTIMIZE = -O3 
#OPTIMIZE += -hfp3
#OPTIMIZE += -ra
## OPTIMIZE += -h list=a
#OPTIMIZE += -h bounds
# #OPTIMIZE += -fprofile-arcs -ftest-coverage  # for gcov
# #OPTIMIZE += -pg 
#OPTIMIZE += -K trap=divz,fp,inv
#OPTIMIZE += -G0 


GRAPE="XC30"
ifeq ($(UNAME),saitoh)
INCLUDEPATH += -I/home/saitoh/lib/avx_phantom
INCLUDEPATH += -I/home/saitoh/lib/gsl/include
#INCLUDEPATH += -I/home/saitoh/lib/CloudyCooling
INCLUDEPATH += -I/home/saitoh/lib/asrch
INCLUDEPATH += -I/home/saitoh/lib/CELib
#INCLUDEPATH += -I/home/saitoh/lib/asrsic
INCLUDEPATH += -I/home/saitoh/lib/asrflx
OPTIONS += -DHAVE_AVX_PHANTOM_GRAPE
LDFLAGS += -L/home/saitoh/lib/avx_phantom -lpg5
LDFLAGS += -L/home/saitoh/lib/gsl/lib -lgsl -lgslcblas
#LDFLAGS += -L/home/saitoh/lib/CloudyCooling -lCloudyCooling
LDFLAGS += -L/home/saitoh/lib/asrch -lASRCH
LDFLAGS += -L/home/saitoh/lib/CELib -lCELib
#LDFLAGS += -L/home/saitoh/lib/asrsic -lASRSIC
LDFLAGS += -L/home/saitoh/lib/asrflx -lASRFLX
endif
ifeq ($(UNAME),saitoutk)
INCLUDEPATH += -I/home/saitoutk/lib/saitoh_phantom
INCLUDEPATH += -I/home/saitoutk/lib/gsl/include
#INCLUDEPATH += -I/home/saitoutk/lib/CloudyCooling
INCLUDEPATH += -I/home/saitoutk/lib/asrch
INCLUDEPATH += -I/home/saitoutk/lib/CELib
INCLUDEPATH += -I/home/saitoutk/lib/asrsic
INCLUDEPATH += -I/home/saitoutk/lib/asrflx
OPTIONS += -DHAVE_PHANTOM_GRAPE
LDFLAGS += -L/home/saitoutk/lib/saitoh_phantom -lpg5pot
LDFLAGS += -L/home/saitoutk/lib/gsl/lib -lgsl -lgslcblas
#LDFLAGS += -L/home/saitoutk/lib/CloudyCooling -lCloudyCooling
LDFLAGS += -L/home/saitoutk/lib/asrch -lASRCH
LDFLAGS += -L/home/saitoutk/lib/CELib -lCELib
LDFLAGS += -L/home/saitoutk/lib/asrsic -lASRSIC
LDFLAGS += -L/home/saitoutk/lib/asrflx -lASRFLX
endif
