CC= g++
CFLAGS= -Wall -O2 -I/usr/local/include
LIBS= -lfftw3 -lm -L/usr/local/lib

fftw_test: fftw_test.o
		$(CC) $(CFLAGS) -o $@ $< $(LIBS)

%.o: %.c
		$(CC) $(CFLAGS) -o $@ $<

clean:
		-rm *.o
