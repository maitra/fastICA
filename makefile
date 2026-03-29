CC =	gcc

CFLAGS = -ansi -std=c11 -Wall -pedantic -O3 -msse4 

#OBJ =	cephes_eigens.o order.o mat_vec.o fastICA.o util.o
OBJ = fastICA.o mat_vec.o power_eigens.o optlist.o util.o

run_fastICA:	run_fastICA.c $(OBJ)
	gcc -o run_fastICA run_fastICA.c $(OBJ) $(CFLAGS) -lm -lRmath

testfastICA:	testfastICA.c  $(OBJ) 
	gcc -o testfastICA testfastICA.c $(OBJ) $(CFLAGS) -lm -lRmath

clean:	
	rm $(OBJ) run_fastICA testfastICA

