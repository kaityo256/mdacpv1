#------------------------------------------------------
# Makefile for C/C++ Program
#------------------------------------------------------

TARGET=mdacp

LDFLAGS=

#------------------------------------------------------
# Default Parameters
#------------------------------------------------------

CC=mpic++
CPPFLAGS=-std=c++11 -O3

#------------------------------------------------------
# Compile Option
#------------------------------------------------------

-include makefile.opt

#------------------------------------------------------
# Definition
#------------------------------------------------------

.SUFFIXES: .c .cc .h. .o

#---
# Source Files
#---

SRC=$(shell ls *.cc)
HED=$(shell ls *.h)
OBJ=$(SRC:.cc=.o)

#------------------------------------------------------
# rules
#------------------------------------------------------

all: $(TARGET)
$(TARGET): $(OBJ)
	$(CC) $(CPPFLAGS) -o $(TARGET) $(OBJ) $(LDFLAGS)

.cc.o:
	$(CC) $(CPPFLAGS) -c $< 

.PHONY: clean dep

dep:
	g++ -MM -MG -std=c++11 $(SRC) >makefile.depend

makefile.depend: 
	g++ -MM -MG -std=c++11 $(SRC) >makefile.depend

clean:
	rm -f $(TARGET) $(OBJ) gmon.*.out gmon.out

tar:
	tar cvzf $(TARGET).tar.gz *.cfg $(SRC) $(HED) makefile

#--------------------------------------------------
-include makefile.depend
