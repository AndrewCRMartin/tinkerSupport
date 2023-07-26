LIBDIR = $(HOME)/lib
INCDIR = $(HOME)/include

CC = cc
OFILES1 = tinkerpatch.o
OFILES2 = fixoverlap.o
LIBS   = -lbiop -lgen -lm -lxml2
#CFLAGS = -g -ansi -Wall -DDEBUG=1
#CFLAGS = -g -ansi -Wall
CFLAGS = -O3 -ansi -Wall

EXE = tinkerpatch fixoverlap

all : $(EXE)

tinkerpatch : $(OFILES1)
	$(CC) $(CFLAGS) -o $@ $(OFILES1) -L $(LIBDIR) $(LIBS)

fixoverlap : $(OFILES2)
	$(CC) $(CFLAGS) -o $@ $(OFILES2) -L $(LIBDIR) $(LIBS)

.c.o :
	$(CC) $(CFLAGS) -c $< -I $(INCDIR)

clean :
	\rm -f $(OFILES1) $(OFILES2)

distclean: clean
	\rm -f $(EXE)
