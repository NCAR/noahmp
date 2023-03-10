# Makefile 
#
.SUFFIXES:
.SUFFIXES: .o .F90

include ../../hrldas/user_build_options

OBJS = Machine.o \
       CheckNanMod.o
#      ErrorHandleMod.o

all:	$(OBJS)

.F90.o:
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(*).F90

#
# This command cleans up object (etc) files:
#

clean:
	rm -f *.o *.mod *.exe

#
# Dependencies:
#
