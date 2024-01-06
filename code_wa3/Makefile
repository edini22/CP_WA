SHELL = /bin/sh
BIN_NAME = MDcuda

CXX = nvcc
LD  = nvcc

CC = gcc

CXXFLAGS   = -O3 -g -std=c++11 -arch=sm_35 -Wno-deprecated-gpu-targets  
# CXXFLAGS   = -O3 -g -std=c++11 -arch=sm_70 -Wno-deprecated-gpu-targets  # local machine

SRC_DIR = src
BIN_DIR = bin
BUILD_DIR = build
SRC = $(wildcard $(SRC_DIR)/MDcuda.cu)
OBJ = $(patsubst src/%.cu,build/%.o,$(SRC))
BIN = $(BIN_NAME)

vpath %.cu $(SRC_DIR)

.DEFAULT_GOAL = all

$(BUILD_DIR)/%.o: %.cu
	$(CXX) -c $(CXXFLAGS) $(INCLUDES) $< -o $@ $(LIBS)

$(BIN_DIR)/$(BIN_NAME): $(OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ $(OBJ) $(LIBS)

checkdirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

all: checkdirs $(BIN_DIR)/$(BIN_NAME)

clean:
	rm -f $(BUILD_DIR)/* $(BIN_DIR)/* 

run: clean all
	sbatch run.sh