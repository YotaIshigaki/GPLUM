#
# Makefile for ASURA.
#

include ./machines/Makefile.switch
include makeflags


CFLAGS = $(OPTIMIZE) $(OPTIONS) $(INCLUDEPATH) 
CFLAGS += -MMD -MP 
LDFLAGS += -lm

SRCDIR    = ./
ifeq "$(strip $(SRCDIR))" ""
  SRCDIR  = .
endif
SOURCES   = $(wildcard $(SRCDIR)/*.c)
OBJDIR    = ./obj
ifeq "$(strip $(OBJDIR))" ""
  OBJDIR  = .
endif
OBJECTS   = $(addprefix $(OBJDIR)/, $(notdir $(SOURCES:.c=.o)))
DEPENDS   = $(OBJECTS:.o=.d)

$(OBJDIR)/%.o: $(SRCDIR)/%.c
	@-mkdir -p $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ -c $<

.PHONY: ALL
ALL: preprocess start

.PHONY: starte
starte: $(OBJECTS) 
	@$(LINKER) $(CPPFLAGS) $(CFLAGS) -o asura.out $^ $(LDFLAGSCPP) $(LDFLAGS)

.PHONY: preprocess
preprocess:
	./AutoLog.rb makeflags;\
	./AutoMakeIO.rb ;\

.PHONY: start
start: $(OBJECTS) 
	@echo '';\
	echo '=============================================================================';\
	echo ' make ASURA';\
	echo ' OPTIMIZATION PARAMETERS = ';\
	echo '   $(OPTIMIZE)';\
	echo ' OPTIONS = ';\
	echo '   $(OPTIONS)';\
	echo ' INCLUDEPATH = ';\
	echo '   $(INCLUDEPATH)';\
	echo ' LIBRARIES =';\
	echo '   $(LDFLAGS)';\
	echo ' Compile mode = $(GRAPE)';\
	echo '=============================================================================';\
	echo '';\
	echo '';\
	make starte ;\
	echo '';\
	echo '=============================================================================';\
	echo ' make ASURA';\
	echo ' OPTIMIZATION PARAMETERS = ';\
	echo '   $(OPTIMIZE)';\
	echo ' OPTIONS = ';\
	echo '   $(OPTIONS)';\
	echo ' INCLUDEPATH = ';\
	echo '   $(INCLUDEPATH)';\
	echo ' LIBRARIES =';\
	echo '   $(LDFLAGS)';\
	echo ' Compile mode = $(GRAPE)';\
	echo ' System = $(SYSTEM)';\
	echo '=============================================================================';\
	echo '';\
	echo ''


.PHONY: clean
tag:
	ctags *.[c,h] AutoMakeIO.rb makeflags

.PHONY: clean
clean:
	-rm -f $(OBJECTS) $(DEPENDS) *.out *.swp 

DATE = $(shell date '+%Y%m%d%H%M')
ARCDIR    = ./archives
.PHONY: pkg
pkg:
	-mkdir -p $(ARCDIR)
	#echo '$(DATE)' 
	tar cvf $(ARCDIR)/asurapkg.$(DATE).tar ./*.[ch] ./*.rb ./CoolingCurveSpaans2008.dat Makefile ./machines/* makeflags 
	gzip $(ARCDIR)/asurapkg.$(DATE).tar

.PHONY: strip
strip:
	strip asura.out

-include $(DEPENDS)
