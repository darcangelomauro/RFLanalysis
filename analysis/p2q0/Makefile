################################################################################
#
# Makefile to compile and link C++ programs
#
# Version valid for Linux machines
#
# "make" compiles and links the specified main programs and modules
# using the specified libraries (if any), and produces the executables
# 
# "make clean" removes all files created by "make"
#
# Dependencies on included files are automatically taken care of
#
################################################################################

all: rmxeq mkdep mkxeq
.PHONY: all


# main programs and required modules 

MAIN = Aij2Bij2 Aij2Bkl2 ABii2 ABij2 ABij4 ABij2il2 ABij2kl2 ABkllmmnnk AB_aggregate AB2 A2B2 AB4 anticomm_AB r2AB2 rA3AB2 r2 AB24_dim_manip A2B2_dim_manip anticomm_AB_dim_manip rAB_dim_manip r2_dim_manip 

SOURCE = params utils geometry clifford statistics

# search path for modules

VPATH = ../../../RFLmain/source ../../../RFL/source ../../source


# additional include directories

INCPATH = /home/pmxmd10/gsl/include ../../../RFLmain/include ../../../RFL/include ../../include


# additional libraries to be included 
 
LIBS = gsl openblas armadillo

LIBPATH = /home/pmxmd10/gsl/lib


# scheduling and optimization options (such as -DSSE -DSSE2 -DP4)
 
CFLAGS = -std=c++14 -pedantic -fstrict-aliasing -Wall -Wno-long-long -Werror 

############################## do not change ###################################

SHELL=/bin/bash

CC=g++

PGMS= $(MAIN) $(SOURCE)

INCDIRS = $(addprefix -I,$(INCPATH))

OBJECTS = $(addsuffix .o,$(SOURCE))

LDFLAGS = $(addprefix -L,$(LIBPATH)) $(addprefix -l,$(LIBS))

-include $(addsuffix .d,$(PGMS))


# rule to make dependencies

$(addsuffix .d,$(PGMS)): %.d: %.cpp Makefile
	@ $(CC) -MM $(INCDIRS) $< -o $@


# rule to compile source programs

$(addsuffix .o,$(PGMS)): %.o: %.cpp Makefile
	$(CC) $< -c $(CFLAGS) $(INCDIRS) -o $@


# rule to link object files

$(MAIN): %: %.o $(OBJECTS) Makefile
	$(CC) $< $(OBJECTS) $(CFLAGS) $(LDFLAGS) -o $@


# produce executables

mkxeq: $(MAIN)


# remove old executables and old error log file

rmxeq:
	@ -rm -f $(MAIN); \
        echo "delete old executables"		


# make dependencies

mkdep:  $(addsuffix .d,$(PGMS))
	@ echo "generate tables of dependencies"


# clean directory 

clean:
	@ -rm -rf *.d *.o .tmp $(MAIN)
.PHONY: clean

################################################################################
