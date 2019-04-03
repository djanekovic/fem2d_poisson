PROJECT = poisson
SRC = $(wildcard ./src/*.c)
OBJ = $(SRC:.c=.o)
CC  = gcc
RM = /bin/rm
CFLAGS = -g -lm -Wall -Wextra
LDFLAGS = -llapacke -ltriangle

all: $(PROJECT)

$(PROJECT): $(OBJ)
	@mkdir -p ./bin
	$(CC) $(CFLAGS) -o ./bin/$@ $^ $(LDFLAGS)

clean:
	$(RM) -f $(OBJ)
