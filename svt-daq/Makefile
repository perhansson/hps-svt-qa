OS := $(shell uname -s)

# Variables
CFLAGS  := -fpermissive -g -Wall `xml2-config --cflags` `root-config --cflags` -I$(BASE)/tracker -I$(BASE)/generic -I$(BASE)/evio -I.
LFLAGS  := `xml2-config --libs` `root-config --libs` -lMinuit -lbz2 -lgsl -lgslcblas

ifeq ($(OS),Linux) #	hack to make this compile on OS X
	LFLAGS += -lrt
endif

CC      := g++
GCC      := gcc
BIN     := $(BASE)/bin
OBJ     := $(BASE)/.obj

# Generic Sources
GEN_DIR := $(BASE)/generic
GEN_SRC := $(GEN_DIR)/Data.cpp $(GEN_DIR)/DataRead.cpp $(GEN_DIR)/XmlVariables.cpp
GEN_OBJ := $(patsubst $(GEN_DIR)/%.cpp,$(OBJ)/%.o,$(GEN_SRC))

# Tracker Sources
TRK_DIR := $(BASE)/tracker
TRK_SRC := $(TRK_DIR)/DevboardEvent.cpp $(TRK_DIR)/DevboardSample.cpp $(TRK_DIR)/DataReadEvio.cpp $(TRK_DIR)/TrackerEvent.cpp $(TRK_DIR)/TrackerSample.cpp $(TRK_DIR)/TriggerEvent.cpp $(TRK_DIR)/TriggerSample.cpp $(TRK_DIR)/TiTriggerEvent.cpp
TRK_OBJ := $(patsubst $(TRK_DIR)/%.cpp,$(OBJ)/%.o,$(TRK_SRC))

# EVIO Sources
EVIO_DIR := $(BASE)/evio
EVIO_SRC := $(EVIO_DIR)/evio.c $(EVIO_DIR)/evioswap.c $(EVIO_DIR)/eviofmt.c $(EVIO_DIR)/eviofmtswap.c
#EVIO_HDR := $(EVIO_DIR)/evio.h $(EVIO_DIR)/evioswap.h
EVIO_OBJ := $(patsubst $(EVIO_DIR)/%.c,$(OBJ)/%.o,$(EVIO_SRC))

# Root Sources
ROOT_DIR := $(BASE)/root
ROOT_SRC := $(wildcard $(ROOT_DIR)/*.cpp)
ROOT_BIN := $(patsubst $(ROOT_DIR)/%.cpp,$(BIN)/%,$(ROOT_SRC))
ROOT_OBJ := $(pathsubst $(ROOT_DIR)/%.cpp,$(OBJ)/%.o $(ROOT_DIR))

# Default
all: test_base dir $(GEN_OBJ) $(TRK_OBJ) $(EVIO_OBJ) $(ROOT_OBJ) $(ROOT_BIN)

# Object directory
dir:
	test -d $(OBJ) || mkdir $(OBJ)

# Object directory
test_base:
ifeq ($(BASE),)
	echo "Need to source setup file."
	exit 1
endif

# Clean
clean:
	rm -rf $(OBJ)

# Compile Generic Sources
$(OBJ)/%.o: $(GEN_DIR)/%.cpp $(GEN_DIR)/%.h
	@echo "$(GEN_DIR)"
	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile Tracker Sources
$(OBJ)/%.o: $(TRK_DIR)/%.cpp $(TRK_DIR)/%.h
	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile EVIO Sources
$(OBJ)/%.o: $(EVIO_DIR)/%.c
	$(GCC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile ROOT Sources
$(OBJ)/%.o: $(ROOT_DIR)/%.cc $(ROOT_DIR)/%.hh
	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile root
$(BIN)/%: $(ROOT_DIR)/%.cpp $(OFF_OBJ) $(GEN_OBJ) $(TRK_OBJ) $(FIT_OBJ) $(ROOT_OBJ)
	$(CC) $(CFLAGS) $(DEF) $(OBJ)/* -o $@ $< $(LFLAGS)

