FC = gfortran
FFLAGS = -Wall -march=native -O5 -fopenmp #compiler flags
LDFLAGS = -fopenmp #link flags

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main #program name

#required objects: 
OBJS =
OBJS += fortplot.o
OBJS += Interactions.o
OBJS += Inits.o
OBJS += Plotroutines.o
OBJS += main_functions.o
OBJS += main.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f95
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
