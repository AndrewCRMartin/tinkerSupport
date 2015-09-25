LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include

CC = cc
OFILES1 = tinkerpatch.o
LIBS   = -lbiop -lgen -lm -lxml2
#CFLAGS = -g -ansi -Wall -DDEBUG=1
#CFLAGS = -g -ansi -Wall
CFLAGS = -O3 -ansi -Wall

EXE = tinkerpatch

all : $(EXE)

tinkerpatch : $(OFILES1)
	$(CC) $(CFLAGS) -o $@ $(OFILES1) -L $(LIBDIR) $(LIBS)

.c.o :
	$(CC) $(CFLAGS) -c $< -I $(INCDIR)

clean :
	\rm -f $(OFILES1)

distclean: clean
	\rm -f $(EXE)
