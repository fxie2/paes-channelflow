include ./configure.user

.SUFFIXES:
.SUFFIXES: .o .F90 .f90 

LIBS = $(FLIBS) $(MPILIBS)
INCLUDES = $(MPIINCLUDES)

MODULES = modules.o parameters.o compbuf.o

LESOBJS = compsgs1.o compsgs3.o comp1.o \
       comp2.o compadv1.o les.o comprhs.o compmn.o abort.o \
       allocate.o default.o lowerbc.o derivs1.o derivs2.o tridv.o \
       upperbc.o read_write.o complen.o buffers.o init.o press.o #\
	#spectra1.o spectra2.o 

OBJS = $(LESOBJS) $(FFTLIB) $(MPIOBJS)

REDIMOBJS = redim.o $(FFTLIB)


all: 	
	@$(MAKE) $(MODULES)
	@$(MAKE) model

model: $(OBJS)
	$(FC) -o model $(FFLAGS) $(OBJS) $(MODULES) $(LIBS)

$(OBJS): $(MODULES)

.f90.o:	
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F90.o:
	$(FC) $(FFLAGS) $(INCLUDES) -c $< $(DEFS)

redim:	$(REDIMOBJS)
	$(FC) -o redim $(FFLAGS) $(REDIMOBJS) $(FLIBS)

clean:
	rm -f *.o *.mod model

