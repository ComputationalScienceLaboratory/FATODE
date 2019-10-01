# Minimun modification of the path of SuperLU source is required.

default: fortran 
# set the original SuperLU direcotry
export SLUPATH = /home/hong/Software/SuperLU_4.2

include $(SLUPATH)/make.inc

HEADER = $(SLUPATH)/SRC

fortran:c_fortran_zgssv.o c_fortran_dgssv.o

c_fortran_zgssv.o: c_fortran_zgssv.c
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $<

c_fortran_dgssv.o: c_fortran_dgssv.c
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $<

clean:
	rm c_fortran*.o
