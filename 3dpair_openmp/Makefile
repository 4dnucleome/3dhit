CC = gcc
CPPFLAGS = -DNOPROF
CFLAGS = -Wall -O3 -ffast-math -fopenmp
LDFLAGS = -lm

PROGRAM = 3dpair
OBJECTS = mod_align.o mod_cont.o mod_func.o mod_move.c

all: $(PROGRAM)

$(PROGRAM): modul.h $(PROGRAM).c $(OBJECTS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) -o $(PROGRAM) $(OBJECTS) $(PROGRAM).c

%.o : %.c modul.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

clean:
	rm -rf $(PROGRAM) *.bak *.o *~

