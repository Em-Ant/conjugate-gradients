#EmAnt Standard Makefile

CC = g++
LDFLAGS= -lm
CFLAGS= -Wall -O3 

OBJS=	Sparse.o\
		SymSparseProb.o\
		main.o
		
EXEC= CG_test_cc


all:	$(OBJS)
	$(CC) -o $(EXEC) $(OBJS) $(LDFLAGS) 
	
clean:
	-@rm -f *.o $(EXEC) *~ *.a

exec: all
	./$(EXEC)

#DEPS gxx -MM 

main.o: main.cpp symsparseprob.h sparse.h
Sparse.o: Sparse.cpp sparse.h
SymSparseProb.o: SymSparseProb.cpp symsparseprob.h sparse.h
