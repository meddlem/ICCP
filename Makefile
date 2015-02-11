FC = gfortran
FFLAGS = -Wall -march=native -O5 #compiler flags
LDFLAGS = #link flags

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main #program name

#required objects: 
OBJS =
OBJS += maxwell.o
OBJS += LJ.o
OBJS += main.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
