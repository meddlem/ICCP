FC = gfortran
FFLAGS = -Wall -Wextra -march=native -O5
LDFLAGS =
LIBS = -llapack
FFLAGS += $(shell pkg-config --cflags plplotd-f95)
LIBS += $(shell pkg-config --libs plplotd-f95)

COMPILE = $(FC) $(FFLAGS)
LINK = $(FC) $(LDFLAGS)

OBJS =
OBJS += PLplotDemo.o

all: PLplotDemo

PLplotDemo: $(OBJS)
	$(LINK) -o $@ $^ $(LIBS)

%.o: %.f90
	$(COMPILE) -o $@ -c $<

.PHONY: clean
clean:
	$(RM) myprog $(OBJS) *.mod
