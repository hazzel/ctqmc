DEFINES+= -DMCL_DUMP_BUFFER=0
DEFINES+= -DMCL_MEASUREMENTS_APPEND=0
DEFINES+= -DMCL_MCL_RNG_MT

MODE=MPI
#MODE=SINGLE
#MODE=PT

OBJS = dump.o parser.o measurements.o evalable.o observable.o random.o mc.o ConfigSpace.o main.o
OBJSLN = dump.LN.o parser.LN.o measurements.LN.o evalable.LN.o observable.LN.o random.LN.o mc.LN.o ConfigSpace.LN.o runner_single.LN.o merge.LN.o

ifeq ($(MODE),MPI)
  OBJS+=runner.o
  ifeq ($(MPICC),)
    MPICC = /usr/lib64/mpi/gcc/openmpi/bin/mpiCC
    #MPICC = /usr/lib64/openmpi/bin/mpiCC
  endif
  CC=$(MPICC)
  LD=$(MPICC)
endif

ifeq ($(MODE),SINGLE)
  OBJS+=runner_single.o
  DEFINES+= -DMCL_SINGLE
  ifneq ($(MCLL_SYSTEM_INFO), rwthcluster)
    CC=g++
    LD=g++
  endif
endif

ifeq ($(MODE),PT)
  OBJS+=runner.o
  DEFINES+= -DMCL_PT
endif

MCLL  = $(HOME)/mc/load_leveller/trunk
APPMCLL = $(HOME)/mc/ctqmc/

ifeq ($(MCLL_SYSTEM_INFO), rwthcluster)
	CFLAGS  = $(FLAGS_FAST) -Wno-deprecated -ansi -std=c++11 -DDEBUG_CXXBLAS -DNDEBUG $(FLAGS_OPENMP) $(DEFINES)
	INCLUDE = $(FLAGS_MATH_INCLUDE) -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/FLENS -DWITH_MKL -DALWAYS_USE_CXXLAPACK
	LDFLAGS = $(FLAGS_MATH_LINKER) $(FLAGS_OPENMP)
	SUPERLP = 
else
        CFLAGS  = -O3 -Wno-deprecated -ansi -ffast-math -std=c++11 -fopenmp $(DEFINES)
        INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/FLENS -DUSE_CXXLAPACK
        LDFLAGS = -fopenmp -llapack
        SUPERLP =
endif

CCLN = g++
LDLN = g++
CFLAGSLN  = $(CFLAGS) -DMCL_SINGLE
INCLUDELN = $(INCLUDE)
LDFLAGSLN = $(LDFLAGS)
SUPERLPLN = $(SUPERLP)

RM = /bin/rm -f

#all: mc merge cleano
all: mc cleano

mc : $(OBJS)
	$(LD) $(LDFLAGS) -o $@ $(OBJS) $(SUPERLP)

merge : $(OBJSLN)
	$(LDLN) $(LDFLAGSLN) -o $@ $(OBJSLN) $(SUPERLPLN)

%.o : $(APPMCLL)/%.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

%.LN.o : $(APPMCLL)/%.cpp
	$(CCLN) -c $(CFLAGSLN) $(INCLUDELN)  $< -o $@

%.o : $(MCLL)/%.cpp
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

%.LN.o : $(MCLL)/%.cpp
	$(CCLN) -c $(CFLAGSLN) $(INCLUDELN) $< -o $@ 

clean:
	$(RM) mc merge $(OBJS) $(OBJSLN)

cleano:
	$(RM) $(OBJS) $(OBJSLN)



