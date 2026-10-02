# Makefile 

.SUFFIXES:
.SUFFIXES: .o .F90

# Detect MPAS build from CPPFLAGS
ifneq (,$(findstring -DNOAHMP_MPAS,$(CPPFLAGS)))
NOAHMP_MPAS_BUILD = yes
endif

# Only HRLDAS/community compilation needs this file.
# It is not needed for "make clean".
ifneq ($(MAKECMDGOALS),clean)
ifndef NOAHMP_MPAS_BUILD
include ../../hrldas/user_build_options
endif
endif


OBJS = Machine.o \
       CheckNanMod.o \
       PiecewiseLinearInterp1dMod.o
#      ErrorHandleMod.o

ifdef NOAHMP_MPAS_BUILD

.PHONY: utility utility_lib
all: dummy utility
dummy:
	@echo "****** compiling physics_noahmp/utility ******"
utility: $(OBJS)
utility_lib:
	ar -ru ./../../libphys.a $(OBJS)

else

all: $(OBJS)

endif

#
# Dependencies:
#
CheckNanMod.o: Machine.o
PiecewiseLinearInterp1dMod.o: Machine.o


ifdef NOAHMP_MPAS_BUILD

.F90.o:
	$(FC) $(CPPFLAGS) $(COREDEF) $(FFLAGS) -c $*.F90 $(CPPINCLUDES) $(FCINCLUDES) -I../../../../framework

else

.F90.o:
	@echo ""
	$(RM) $(*).f90
	$(CPP) $(CPPFLAGS) $(*).F90 > $(*).f90
	$(COMPILERF90) -c $(F90FLAGS) $(FREESOURCE) $(NETCDFMOD) $(*).f90
	@echo ""

endif

clean:
	$(RM) *.o *.mod *.exe *.f90 *.i

