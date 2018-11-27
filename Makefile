SRC = ./src/
BIN = ./bin/
CC  = gcc
RM = /bin/rm
CFLAGS = -O3 -lm -Wall -Wextra -llapacke

all: $(BIN)poisson

$(BIN)poisson: $(SRC)main.c $(SRC)util.c
	$(CC) $(CFLAGS) -o $(BIN)poisson $(SRC)main.c $(SRC)util.c $(SRC)triangle/triangle.o

clean:
	$(RM) -f $(BIN)poisson
