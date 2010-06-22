	# Makefile for 2D FDTD Program

CC = /usr/bin/gcc
#CC = icc
OBJS = fdtd.o 
##CFLAGS = -g -Wall -O2 -pedantic
CFLAGS = -g -Wall -O2 -fopenmp
DEBUGFLAG = -pg
INCLUDES = -I. -I/opt/local/include
LIBS = -L. -L/opt/local/lib 
TARG = fdtd

all: $(TARG)
$(TARG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -g -o $@ $(OBJS)

.c.o: 
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -g -o $@ -c $<

.PHONY: clean
clean:
	rm $(OBJS) $(TARG)

test:
	./$(TARG)
