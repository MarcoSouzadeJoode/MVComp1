FC = gfortran
FFLAGS = -O2 -Wall -Wextra -fbounds-check -g
LIBS = /home/marco/lapack-3.11.0/liblapack.a /home/marco/lapack-3.11.0/librefblas.a /home/marco/lapack-3.11.0/libtmglib.a



EXEC = week7
SRCS = $(EXEC).f90

all: $(EXEC)

$(EXEC): $(SRCS)
	$(FC) -o $(EXEC) $(SRCS) $(LIBS)

	
clean:
	rm -f $(EXEC)