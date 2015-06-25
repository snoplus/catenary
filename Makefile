OPTS = -DLINUX -g -Wno-write-strings -Wno-deprecated

CC = g++

all:	catenary parser

zip:
	tar -cf catenary.tar *.cpp *.h Makefile catenary parser
	gzip catenary.tar

clean:
	rm *.o catenary parser

parser:	parser.o
	$(CC) -o parser parser.o -lm

catenary: catenary.o
	$(CC) -o catenary catenary.o -lm

parser.o: parser.cpp catenary.h
	$(CC) $(OPTS) -c parser.cpp

catenary.o: catenary.cpp catenary.h
	$(CC) $(OPTS) -c catenary.cpp
