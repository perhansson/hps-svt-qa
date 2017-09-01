OS := $(shell uname -s)

# Variables
CFLAGS  := -fpermissive -g -Wall `xml2-config --cflags` `root-config --cflags` -I$(BASE)/tracker -I$(BASE)/generic -I$(BASE)/evio -I.
LFLAGS  := `xml2-config --libs` `root-config --libs` -lMinuit -lbz2 -lgsl -lgslcblas

ifeq ($(OS),Linux) #	hack to make this compile on OS X
	LFLAGS += -lrt
endif

CC      := g++
GCC      := gcc
BIN     := $(PWD)/bin
OBJ     := $(BASE)/.obj

# Example Sources
EXAMPLES_DIR := $(PWD)/examples
EXAMPLES_SRC := $(wildcard $(EXAMPLES_DIR)/*.cpp)
EXAMPLES_BIN := $(patsubst $(EXAMPLES_DIR)/%.cpp,$(BIN)/%,$(EXAMPLES_SRC))
EXAMPLES_OBJ := $(pathsubst $(EXAMPLES_DIR)/%.cpp,$(OBJ)/%.o $(EXAMPLES_DIR))

# Default
all: test_base dir $(EXAMPLES_OBJ) $(EXAMPLES_BIN)

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

# Compile EXAMPLES Sources
$(OBJ)/%.o: $(EXAMPLES_DIR)/%.cc $(EXAMPLES_DIR)/%.hh
	$(CC) -c $(CFLAGS) $(DEF) -o $@ $<

# Compile root
$(BIN)/%: $(EXAMPLES_DIR)/%.cpp $(OFF_OBJ) $(GEN_OBJ) $(TRK_OBJ) $(FIT_OBJ) $(EXAMPLES_OBJ)
	$(CC) $(CFLAGS) $(DEF) $(OBJ)/* -o $@ $< $(LFLAGS)

