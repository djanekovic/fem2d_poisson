SRC = ./src/
BIN = ./bin/
CC  = gcc
RM = /bin/rm
CFLAGS = -g -lm -Wall -Wextra

all: $(BIN)poisson

$(BIN)poisson: $(SRC)main.c
	$(CC) $(CFLAGS) -o $(BIN)poisson $(SRC)main.c $(SRC)triangle/triangle.o

clean:
	$(RM) -f $(BIN)poisson
