DEFINES+= -DMCL_DUMP_BUFFER=0
DEFINES+= -DMCL_MEASUREMENTS_APPEND
DEFINES+= -DMCL_MCL_RNG_MT

MODE=MPI
#MODE=SINGLE
#MODE=PT

OBJS = dump.o parser.o measurements.o evalable.o observable.o random.o mc.o main.o
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
	OBJS+=runner_pt.o
	DEFINES+= -DMCL_PT
	ifeq ($(MPICC),)
		MPICC = /usr/lib64/mpi/gcc/openmpi/bin/mpiCC
		#MPICC = /usr/lib64/openmpi/bin/mpiCC
	endif
	CC=$(MPICC)
	LD=$(MPICC)
endif

MCLL  = $(HOME)/mc/load_leveller/trunk
APPMCLL = $(CURDIR)

ifeq ($(MCLL_SYSTEM_INFO), rwthcluster)
	CFLAGS  = -O3 -DMKL_DIRECT_CALL -ip -axCORE-AVX2,AVX,SSE4.2,SSE4.1 -fp-model fast=2 -Wno-deprecated -std=c++11 -pipe $(DEFINES)
	INCLUDE = -openmp -I${MKLROOT}/include -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/armadillo-5.200.2/include -I$(HOME)/gperftools-2.4/install/include
	LDFLAGS = -L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -openmp
	SUPERLP = 
else ifeq ($(MCLL_SYSTEM_INFO), juqueen)
	CFLAGS  = -O3 -std=c++0x -DMCL_JUQUEEN $(DEFINES)
	INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/gperftools-2.4/install/include
	LDFLAGS = 
	SUPERLP =
else ifeq ($(MCLL_SYSTEM_INFO), desktop_home)
	CFLAGS  = -O3 -Wno-deprecated -std=c++11 -pipe $(DEFINES)
	INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/libs/eigen/ -I$(HOME)/armadillo-5.200.2/include -I$(HOME)/libs/gperftools-2.4/install/include
	LDFLAGS = -L/usr/lib64/atlas -latlas -lf77blas -llapack
	SUPERLP =
else
	CFLAGS  = -Ofast -ffast-math -march=native -flto -fwhole-program -Wno-deprecated -pipe -g -std=c++11 $(DEFINES)
	INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/armadillo-5.200.2/include -I$(HOME)/gperftools-2.4/install/include
	#LDFLAGS = -L$(HOME)/gperftools-2.4/install/lib -Wl,-rpath=$(HOME)/OpenBLAS/lib/ -L$(HOME)/OpenBLAS/lib/ -L$(HOME)/armadillo-5.200.2/lib -lblas -llapack -lprofiler
	LDFLAGS = -L$(HOME)/gperftools-2.4/install/lib -L$(HOME)/armadillo-5.200.2/lib -lblas -llapack -lprofiler
	#LDFLAGS = -L$(HOME)/gperftools-2.4/install/lib -Wl,-rpath=$(HOME)/OpenBLAS/lib/ -L$(HOME)/OpenBLAS/lib/ -lblas -llapack -lprofiler
	SUPERLP =
endif

ifeq ($(MCLL_SYSTEM_INFO), rwthcluster)
	CCLN = icc
	LDLN = icc
else
	CCLN = g++
	LDLN = g++
endif
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


