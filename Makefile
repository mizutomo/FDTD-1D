	# Makefile for 2D FDTD Program

CC = /usr/bin/gcc
#CC = icc
OBJS = fdtd.o utility.o
##CFLAGS = -g -Wall -O2 -pedantic
CFLAGS = -g -Wall -O3 -fopenmp
DEBUGFLAG = -pg
INCLUDES = -I.
LIBS = -L.
TARG = fdtd

all: $(TARG) utility.so
$(TARG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -lm -g -o $@ $(OBJS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -fPIC -g -o $@ -c $<

utility.so: utility.o
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -shared -fPIC -o $@ $<

.PHONY: clean
clean:
	rm $(OBJS) $(TARG)

test_c:
	./$(TARG)

test_py:
	./$(TARG).py
