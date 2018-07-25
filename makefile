CC = gcc
CFLAGS = -O2 
DEPS = sirah_radsup.h 

OBJ = inti_radsup.o mimitian_radsup.o carana_radsup.o mt19937ar.o

all: $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -lm -lgsl -lgslcblas  

debug: 
	$(MAKE) CFLAGS="-g"

inti_radsup.o: inti_radsup.c $(DEPS)
	$(CC) -c $(CFLAGS) inti_radsup.c

mimitian_radsup.o: mimitian_radsup.c $(DEPS) 
	$(CC) -c $(CFLAGS) mimitian_radsup.c

carana_radsup.o: carana_radsup.c $(DEPS) 
	$(CC) -c $(CFLAGS) carana_radsup.c

mt19937ar.o: mt19937ar.c $(DEPS) 
	$(CC) -c $(CFLAGS) mt19937ar.c

clean:
	rm -rf $(OBJ)




