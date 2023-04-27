CC=gcc
LDFLAGS=-lm -lgsl -lgslcblas

SOURCES=xdat-ent-rhoALL_quantum_mod.c xdat-ent-rhoALL_quantum_npt.c xdat-ent-rhoALL_quantum_alloy.c force2dos.c force2band.c
EXECUTABLE=$(SOURCES:.c=)

all : $(SOURCE) $(EXECUTABLE)

.PHONY : all

$(EXECUTABLE) : $(SOURCES)
	$(CC) $(@).c -o $@ $(LDFLAGS) 

clean :
	rm $(EXECUTABLE)
