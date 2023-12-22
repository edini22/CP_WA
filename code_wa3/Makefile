CC = gcc
SRC = src/
CFLAGS =  -pg -g -fno-omit-frame-pointer -O2 -ftree-vectorize -msse4 -mavx

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -lstdc++ -o MDpar.exe

clean:
	rm ./MD*.exe

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=25;\
	./MDpar.exe < inputdata.txt
