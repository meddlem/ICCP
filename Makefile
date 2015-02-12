FC = gfortran
FFLAGS = -Wall -march=native -O5 #compiler flags
LDFLAGS = #link flags

FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

PROG = main #program name

#required objects: 
OBJS =
OBJS += functions.o
OBJS += Incs.o
OBJS += Inits.o
OBJS += plotpoints.o
OBJS += main.o

all: $(PROG)

main: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) $(PROG) $(OBJS) *.mod
