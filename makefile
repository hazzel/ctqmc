DEFINES+= -DMCL_DUMP_BUFFER=0
DEFINES+= -DMCL_MEASUREMENTS_APPEND=0
DEFINES+= -DMCL_MCL_RNG_MT

#MODE=MPI
MODE=SINGLE
#MODE=PT

OBJS = dump.o parser.o measurements.o evalable.o observable.o random.o mc.o ConfigSpace.o main.o
OBJSLN = dump.LN.o parser.LN.o measurements.LN.o evalable.LN.o observable.LN.o random.LN.o mc.LN.o ConfigSpace.LN.o runner_single.LN.o merge.LN.o

ifeq ($(MODE),MPI)
  OBJS+=runner.o
endif

ifeq ($(MODE),SINGLE)
  OBJS+=runner_single.o
  DEFINES+= -DMCL_SINGLE
endif

ifeq ($(MODE),PT)
  OBJS+=runner.o
  DEFINES+= -DMCL_PT
endif

MCLL  = $(HOME)/mc/load_leveller/trunk
APPMCLL = $(HOME)/mc/ctqmc/trunk

ifeq ($(MPICC),)
  MPICC = /usr/lib64/mpi/gcc/openmpi/bin/mpiCC
  #MPICC = /usr/lib64/openmpi/bin/mpiCC
endif
CC = $(MPICC)
LD = $(MPICC)
ifeq ($(MODE),SINGLE)
  CC=g++
  LD=g++
endif
CFLAGS  = -O3 -Wno-deprecated -ansi -ffast-math -std=c++11 -fopenmp $(DEFINES)
INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/libs/eigen/ -I$(HOME)/libs/FLENS -DWITH_ATLAS -DALWAYS_USE_CXXLAPACK
LDFLAGS = -L/usr/lib64/atlas/ -fopenmp -llapack -lf77blas -lcblas -latlas
SUPERLP = 

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



