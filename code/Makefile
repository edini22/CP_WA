CC = gcc
SRC = src/
CFLAGS = -pg -g -fno-omit-frame-pointer -O2 -ftree-vectorize -msse4 -mavx

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm ./MD.exe

run:
	./MD.exe < inputdata.txt

