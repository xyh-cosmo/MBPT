EXEC    = MBPT
OP      = ./src/
INC     = ./include/
OBJS    = $(OP)mbpt_main.o\
          $(OP)mbpt_background.o\
          $(OP)mbpt_ic.o\
          $(OP)mbpt_ode.o\
          $(OP)mbpt_gsl_tools.o\
          $(OP)mbpt_power.o

INCL    = $(INC)mbpt.h

CC      = /usr/bin/gcc
OPTIMIZE    = -o3 -Wall
LIBS        = -lm -lgsl -lgslcblas

$(EXEC):$(OBJS)
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)

.PHONY:clean tidy run
clean:
	rm -f $(OBJS)
tidy:
	rm -f $(OBJS) $(EXEC)
run:
	time ./$(EXEC)
auto:
	make tidy; make; make run
