ifeq ($(RELEASE), 1)
DBG_FLAGS = -O2 -DNDEBUG
else
DBG_FLAGS = -g -O0 -DDEBUG
endif

LDFLAGS =  -L/usr/local/lib
CFLAGS = -I/usr/local/include
LIBFLAGS = -lcrypto -lssl -lm
CC = gcc
CFILES := $(shell find src -name '*.c' | sed -e 's/\.c/\.o/')

OBJS = ${CFILES}

sign: ${OBJS} ./PQCgenKAT_sign.c
	$(CC) $(LDFLAGS) $(CFLAGS) $(DBG_FLAGS) -o $@ $^  $(LIBFLAGS)

runtime: ${OBJS} ./runtime.c
	$(CC) $(LDFLAGS) $(CFLAGS) $(DBG_FLAGS) -o $@ $^  $(LIBFLAGS)

%.o : %.c
	$(CC) $(CFLAGS) -o $@ -c $< 

rng.o : rng.h
api.o : api.h

clean:
	cd src; rm -f *.o; cd ..
	rm -f *.o
	rm -f sign runtime