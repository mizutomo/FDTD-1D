	# Makefile for 2D FDTD Program

CC = /usr/bin/gcc
#CC = icc
OBJS = fdtd.o utility.dylib
##CFLAGS = -g -Wall -O2 -pedantic
CFLAGS = -g -Wall -O2 -fopenmp
DEBUGFLAG = -pg
INCLUDES = -I. -I/opt/local/include
LIBS = -L. -L/opt/local/lib 
TARG = fdtd

all: $(TARG)
$(TARG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -m32 -g -o $@ $(OBJS)

.c.o: 
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -m32 -g -o $@ -c $<

utility.dylib: utility.c
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -dynamiclib -m32 -o $@ $<

.PHONY: clean
clean:
	rm $(OBJS) $(TARG)

test_c:
	./$(TARG)

test_py:
	./$(TARG).py
