CC=g++

MAKE_EXE=-o exe

CFLAGS=-O2 -L/usr/X11R6/lib -lm -lpthread -lX11 -std=c++11 -Dcimg_use_png -lpng  -lz

output:	./includes/voraldo/v.cc ./includes/voraldo/v.h main.cc
	$(CC) $(MAKE_EXE) main.cc ./includes/voraldo/v.cc $(CFLAGS)
