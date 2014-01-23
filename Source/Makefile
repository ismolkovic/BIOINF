CC = g++
CFLAGS = -Wall -O3

INCS = -I ./inc
SRCS = lib/dm.cc lib/nj.cc lib/pt.cc
OBJS = $(SRCS:.cc=.o)
OUTDIR = ../Bin/

#project name
MAIN = ../Bin/LAB

all:	$(MAIN)

$(MAIN): $(OBJS)
	test -d $(OUTDIR) || mkdir $(OUTDIR)
	$(CC) $(CFLAGS) $(INCS) -o $(MAIN) $(OBJS) main.cc

.cc.o:
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@ 

clean:
	rm -rf $(MAIN) lib/*.o
