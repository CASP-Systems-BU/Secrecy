CC= gcc -std=c99
CFLAGS= -O3 -Wall -lsodium
MPI= mpicc
SRC=../src
PRIMITIVES= $(SRC)/comm.c $(SRC)/party.c $(SRC)/primitives.c $(SRC)/sharing.c $(SRC)/utils.c
RELATIONAL= $(SRC)/relational.c
RELATIONAL= $(SRC)/relational.c

all:	comorbidity join-distinct

comorbidity:	comorbidity.c $(PRIMITIVES) $(RELATIONAL)
			$(MPI) $(CFLAGS) -o comorbidity comorbidity.c $(PRIMITIVES) $(RELATIONAL)

join-distinct:	join_distinct.c $(PRIMITIVES) $(RELATIONAL)
			$(MPI) $(CFLAGS) -o join-distinct join_distinct.c $(PRIMITIVES) $(RELATIONAL)

clean:
			rm -f comorbidity join-distinct *.o *.exec
