CC = gcc
SRC = src/
CFLAGS = -pg -g -fno-omit-frame-pointer -O2 -ftree-vectorize -msse4 -mavx

.DEFAULT_GOAL = all

all: MDseq.exe MDpar.exe

local: MDseq_pc.exe MDpar_pc.exe

localpar: MDpar_pc.exe

MDseq.exe: $(SRC)/MDseq.cpp
			module load gcc/11.2.0;
			$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe		

MDpar.exe: $(SRC)/MDpar.cpp
			module load gcc/11.2.0;
			$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -lstdc++ -o MDpar.exe
			
MDseq_pc.exe: $(SRC)/MDseq.cpp
			$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe	
			
MDpar_pc.exe: $(SRC)/MDpar.cpp
			$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -lstdc++ -o MDpar.exe

clean:
	rm ./MD*.exe

runseq:
	srun --partition=cpar perf stat ./MDseq.exe < inputdata.txt

runpar:
	srun --partition=cpar --cpus-per-task=32 perf stat  ./MDpar.exe < inputdata.txt

runparpc:
	perf stat ./MDpar.exe < inputdata.txt
	
runseqpc:
	perf stat ./MDseq.exe < inputdata.txt
