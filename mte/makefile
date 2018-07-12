CXX = g++
C11 = ${CXX} -std=c++11
current_dir = $(notdir $(shell pwd))
parent_dir = $(notdir ${current_dir}/..)
.PHONY : clean
# setting up suffix rules
.SUFFIXES : .cc .h .o
# for creating object files
.cc.o :
	${C11} $*.cc -c -o $*.o

##<<Different commands for makefile>>

## all : runs main code and creates figure
all : shift

## run : runs the numerical calculations only
shift : shiftStochasticity.py
	python shiftStochasticity.py

## clean : remove auto-generated files
clean :
	rm -f *.pyc

help : makefile
	@sed -n 's/^##//p' $<
