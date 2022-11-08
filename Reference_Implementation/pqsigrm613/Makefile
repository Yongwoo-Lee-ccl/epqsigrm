LDFLAGS =  -L/usr/local/lib
CFLAGS = -I/usr/local/include
LIBFLAGS = -lcrypto -lssl -lm
CC = gcc
CFILES := $(shell find src -name '*.c' | sed -e 's/\.c/\.o/')

OBJS = ${CFILES}

sign: ${OBJS} ./PQCgenKAT_sign.c
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $^  $(LIBFLAGS)

runtime: ${OBJS} ./runtime.c
	$(CC) $(LDFLAGS) $(CFLAGS) -o $@ $^  $(LIBFLAGS)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< 

rng.o : rng.h
api.o : api.h

clean:
	rm -f *.o
	rm -f sign runtime