CC = gcc
LDFLAGS= -lm
CFLAGS= -Wall -O3 #-pedantic -ansi

OBJS=	CG.o\
		main.o
		
EXEC= CG_test_c


all:	$(OBJS)
	cc -o $(EXEC) $(OBJS) $(LDFLAGS) 
	
clean:
	-@rm -f *.o $(EXEC) *~ *.a

exec: all
	./$(EXEC)

CG.o:	CG.c CG.h
main.o:	CG.h
