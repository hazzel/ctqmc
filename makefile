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

mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
current_dir := $(notdir $(patsubst %/,%,$(dir $(mkfile_path))))
MCLL  = $(HOME)/mc/load_leveller/trunk
APPMCLL = $(CURDIR)
USE_MKL = FALSE

ifeq ($(MCLL_SYSTEM_INFO), rwthcluster)
	CFLAGS  = $(FLAGS_FAST) -Wno-deprecated -std=c++11 -DNDEBUG $(DEFINES)
	INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/gperftools-2.4/install/include
	LDFLAGS = 
	SUPERLP = 
	ifeq ($(USE_MKL), TRUE)
		CFLAGS += -DEIGEN_USE_MKL_ALL
		INCLUDE += $(FLAGS_MATH_INCLUDE)
		LDFLAGS += $(FLAGS_MATH_LINKER)
	endif
else ifeq ($(MCLL_SYSTEM_INFO), desktop_home)
	CFLAGS  = -O3 -Wno-deprecated -std=c++11 $(DEFINES)
	INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/libs/eigen/
	LDFLAGS =
	SUPERLP =
else
#	CFLAGS  = -O3 -ffast-math -march=native -pipe -g -Wall $(DEFINES)
	CFLAGS  = -Ofast -ffast-math -march=native -flto -fwhole-program -Wno-deprecated -pipe -std=c++11 -g $(DEFINES)
#	CFLAGS  = -Ofast -ffast-math -flto -fwhole-program -Wno-deprecated -pipe -std=c++11 -g $(DEFINES)
	INCLUDE = -I$(MCLL) -I$(APPMCLL) -I$(HOME)/eigen/ -I$(HOME)/gperftools-2.4/install/include
	LDFLAGS = -L$(HOME)/gperftools-2.4/install/lib -lprofiler
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



