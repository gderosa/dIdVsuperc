CC = gcc
CFLAGS = \
  -Wall -W -g -O2 -Wmissing-prototypes -Wshadow -Wpointer-arith \
  -Wcast-qual -Wcast-align -Wwrite-strings -fno-common \
  -ansi -pedantic 
EXE = dIdVsuperc
OBJECTS = misc.o globals.o init.o ui.o functions.o fit.o plot.o splot.o main.o 
LDFLAGS = 
#LIBS = -lm -lgsl -lcblas -latlas
LIBS = -lm -lgsl -lgslcblas

all: dIdVsuperc

dIdVsuperc: $(OBJECTS)
	$(CC) -o $(EXE) $(OBJECTS) $(LDFLAGS) $(LIBS)

strip: $(EXE)
	strip $(EXE)

.SUFFIXES: .c .o

.c.o:  
	$(CC) -o $@ -c $(CFLAGS) $<

clean:
	rm -f *.plt *.o $(EXE)

.PHONY: all
.PHONY: clean
