CC= g++
CFLAGS= -Wall -O2 -I/usr/local/include -I/usr/local/include/opencv
LIBS= -lfftw3 -lm `pkg-config --libs opencv`
OBJ= main.o

all: $(OBJ)
		$(CC) $(CFLAGS) -o comsim $< $(LIBS)

fftw_test: fftw_test.o
		$(CC) $(CFLAGS) -o $@ $< $(LIBS)

%.o: %.c
		$(CC) $(CFLAGS) -o $@ $<

clean:
		-rm *.o
