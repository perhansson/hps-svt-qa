# Variables
#CFLAGS  := -g -Wall `xml2-config --cflags` `root-config --cflags` -I$(PWD)/../tracker -I$(PWD)/../generic -I$(PWD)/../t0fit -I$(PWD)/../evio -I.
CFLAGS  := -g -Wall `xml2-config --cflags` `root-config --cflags` -I$(PWD)/../tracker -I$(PWD)/../generic -I$(PWD)/../t0fit 
LFLAGS  := `xml2-config --libs` `root-config --libs` -lMinuit -lbz2 -lrt -lgsl -lgslcblas
CC      := g++
GCC      := gcc
BIN     := $(PWD)/../bin
OBJ     := $(PWD)/.obj

# Generic Sources
GEN_DIR := $(PWD)/../generic
GEN_SRC := $(GEN_DIR)/Data.cpp $(GEN_DIR)/DataRead.cpp $(GEN_DIR)/XmlVariables.cpp
#GEN_HDR := $(GEN_DIR)/Data.h   $(GEN_DIR)/DataRead.h $(GEN_DIR)/XmlVariables.h
GEN_OBJ := $(patsubst $(GEN_DIR)/%.cpp,$(OBJ)/%.o,$(GEN_SRC))

# Tracker Sources
TRK_DIR := $(PWD)/../tracker
#TRK_SRC := $(TRK_DIR)/DevboardEvent.cpp $(TRK_DIR)/DevboardSample.cpp $(TRK_DIR)/DataReadEvio.cpp $(TRK_DIR)/TrackerEvent.cpp 
TRK_SRC := $(TRK_DIR)/DevboardEvent.cpp $(TRK_DIR)/DevboardSample.cpp $(TRK_DIR)/DataReadEvio.cpp $(TRK_DIR)/TrackerEvent.cpp $(TRK_DIR)/TrackerSample.cpp $(TRK_DIR)/TriggerEvent.cpp $(TRK_DIR)/TriggerSample.cpp
#TRK_HDR := $(TRK_DIR)/DevboardEvent.h   $(TRK_DIR)/DevboardSample.h $(TRK_DIR)/DataReadEvio.h
TRK_OBJ := $(patsubst $(TRK_DIR)/%.cpp,$(OBJ)/%.o,$(TRK_SRC))

# Fitter Sources
FIT_DIR := $(PWD)/../t0fit
FIT_SRC := $(FIT_DIR)/ShapingCurve.cc $(FIT_DIR)/Samples.cc $(FIT_DIR)/Event.cc $(FIT_DIR)/Fitter.cc $(FIT_DIR)/AnalyticFitter.cc $(FIT_DIR)/LinFitter.cc $(FIT_DIR)/SmoothShapingCurve.cc
#FIT_HDR := $(FIT_DIR)/ShapingCurve.hh $(FIT_DIR)/Samples.hh $(FIT_DIR)/Event.hh $(FIT_DIR)/Fitter.hh $(FIT_DIR)/AnalyticFitter.hh $(FIT_DIR)/LinFitter.hh $(FIT_DIR)/SmoothShapingCurve.hh
FIT_OBJ := $(patsubst $(FIT_DIR)/%.cc,$(OBJ)/%.o,$(FIT_SRC))

# EVIO Sources
#EVIO_DIR := $(PWD)/../evio
#EVIO_SRC := $(EVIO_DIR)/evio.c $(EVIO_DIR)/evioswap.c $(EVIO_DIR)/eviofmt.c $(EVIO_DIR)/eviofmtswap.c
#EVIO_HDR := $(EVIO_DIR)/evio.h $(EVIO_DIR)/evioswap.h
#EVIO_OBJ := $(patsubst $(EVIO_DIR)/%.c,$(OBJ)/%.o,$(EVIO_SRC))
CODA     := /u1/coda_rth
EVIO_DIR := $(CODA)/coda/src/evio-4.0
EVIO_LIB := $(EVIO_DIR)/Linux_i686/obj
#EVIO_SRC := $(addprefix $(EVIO_DIR)/evio.c $(EVIO_DIR)/evioswap.c $(EVIO_DIR)/eviofmt.c $(EVIO_DIR)/eviofmtswap.c
EVIO_SRC := evio.c evioswap.c eviofmt.c eviofmtswap.c
EVIO_SRC2 := $(addprefix $(EVIO_LIB)/,$(EVIO_SRC))
EVIO_OBJ := $(patsubst %.c,%.o,$(EVIO_SRC2))


#EVIO_HDR := $(EVIO_DIR)/evio.h $(EVIO_DIR)/evioswap.h
#EVIO_OBJ := $(patsubst $(EVIO_DIR)/%.c,$(OBJ)/%.o,$(EVIO_SRC))
EVIO_INC := $(EVIO_DIR)/src.s
EVIO_OBJS2 := $(wildcard $(EVIO_LIB)/*.o)
EVIO_OBJS := $(addprefix -l,$(notdir $(wildcard $(EVIO_LIB)/*.o)))
#EVIO_OBJ := $(patsubst %.o,%,$(EVIO_OBJS))
#EVIO_OBJ := $(EVIO_LIB)/%.o

# Root Sources
ROOT_DIR := $(PWD)
#ROOT_SRC := $(wildcard $(ROOT_DIR)/*.cpp)
ROOT_SRC := $(ROOT_DIR)/prbs_verify.cpp
ROOT_BIN := $(patsubst $(ROOT_DIR)/%.cpp,$(BIN)/%,$(ROOT_SRC))
#ROOT_OBJ := $(OBJ)/meeg_utils.o $(OBJ)/cosmic_utils.o
ROOT_OBJ := 


# Default
#all: dir $(GEN_OBJ) $(OFF_OBJ) $(TRK_OBJ) $(FIT_OBJ) $(EVIO_OBJ) $(ROOT_OBJ) $(ROOT_BIN)
all: dir $(ROOT_OBJ) $(ROOT_BIN)
#all: dir $(ROOT_BIN)

# Object directory
dir:
	test -d $(OBJ) || mkdir $(OBJ)

# Clean
clean:
	rm -rf $(OBJ)

# Compile Generic Sources
#$(OBJ)/%.o: $(GEN_DIR)/%.cpp $(GEN_DIR)/%.h
#	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile Tracker Sources
#$(OBJ)/%.o: $(TRK_DIR)/%.cpp $(TRK_DIR)/%.h
#	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

## Compile Fitter Sources
#$(OBJ)/%.o: $(FIT_DIR)/%.cc $(FIT_DIR)/%.hh
#	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile EVIO Sources
#$(OBJ)/%.o: $(EVIO_DIR)/%.c
#	$(GCC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile ROOT Sources
#$(OBJ)/%.o: $(ROOT_DIR)/%.cc $(ROOT_DIR)/%.hh
#	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile root
#$(BIN)/%: $(ROOT_DIR)/%.cpp $(OFF_OBJ) $(GEN_OBJ) $(TRK_OBJ) $(FIT_OBJ) $(ROOT_OBJ)
#$(BIN)/%: $(ROOT_DIR)/%.cpp $(OFF_OBJ) $(GEN_OBJ) $(TRK_OBJ) $(ROOT_OBJ)
$(BIN)/%: $(ROOT_DIR)/%.cpp $(ROOT_OBJ)
	@echo "Compiling $(BIN)/$@"
#$(CC) $(CFLAGS) $(DEF) $(OBJ)/* -o $@ $< $(LFLAGS)
#	$(CC) $(CFLAGS) $(DEF) -I$(EVIO_INC) -o $@ $< $(LFLAGS) -L$(EVIO_LIB) $(EVIO_OBJ)
	$(CC) $(CFLAGS) $(DEF) -I$(EVIO_INC) -o $@ $< $(LFLAGS) $(EVIO_OBJ)

