CC = gcc
SRC = src/
CFLAGS = -pg -g -fno-omit-frame-pointer

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)/MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

clean:
	rm ./MD.exe

run:
	./MD.exe < inputdata.txt
