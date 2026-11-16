CC      = gcc
CFLAGS  = -O3 -march=native -mavx -fopt-info-vec
LDLIBS  = -lm
 
# Para instrumentar com LIKWID:  make LIKWID=1  (LIKWID_HOME padrão: /usr/local)
LIKWID_HOME ?= /usr/local
ifdef LIKWID
CFLAGS += -DLIKWID_PERFMON -I$(LIKWID_HOME)/include
LDLIBS += -L$(LIKWID_HOME)/lib -llikwid
endif
 
OBJS = main.o utils.o funcao.o jacobiana.o solver.o newton.o
 
all: broyden
 
broyden: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDLIBS)
 
main.o:      main.c newton.h utils.h
utils.o:     utils.c utils.h
funcao.o:    funcao.c funcao.h
jacobiana.o: jacobiana.c jacobiana.h
solver.o:    solver.c solver.h jacobiana.h
newton.o:    newton.c newton.h funcao.h jacobiana.h solver.h utils.h
 
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
 
clean:
	rm -f *.o broyden
 
