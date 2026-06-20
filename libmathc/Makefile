CC = gcc
CFLAGS = -std=c99 -fgnu89-inline -Wall -Wextra -O3 -march=native
LDFLAGS = -lm

TARGET = test
SRC = test.c

all: $(TARGET)

$(TARGET): $(SRC) *.h
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
