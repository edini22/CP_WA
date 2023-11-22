CC = gcc
SRC = src/
CFLAGS = -pg -g -fno-omit-frame-pointer -O2 -ftree-vectorize -msse4 -mavx

# .DEFAULT_GOAL = all
.DEFAULT_GOAL = gerpar

all: MDseq.exe MDpar.exe

gerpar: MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
			module load gcc/11.2.0;
			$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
			module load gcc/11.2.0;
			$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -lstdc++ -o MDpar.exe

clean:
	rm ./MD*.exe

runseq:
	srun --partition=cpar ./MDseq.exe < inputdata.txt

runpar:
	srun --partition=cpar --cpus-per-task=32 ./MDpar.exe < inputdata.txt
