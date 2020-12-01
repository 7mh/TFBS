# Note - you should compile with this make file to make sure your files
#  are compiled into the right file names for the check program

CC = g++
CARGS = -Wall -g

all: tf8merFilter tf8merAlign

tf8merFilter: tf8merFilter.cpp mer8.h seq.h
	$(CC) $(CARGS) tf8merFilter.cpp -o tf8merFilter -lm

tf8merAlign: tf8merAlign.cpp mer8.h seq.h
	$(CC) $(CARGS) tf8merAlign.cpp -o tf8merAlign -lm

