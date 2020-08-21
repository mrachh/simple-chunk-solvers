
EXEC = int3

FC = gfortran
FFLAGS = -O2 -c -w -std=legacy 
FLINK = gfortran -w -o $(EXEC) 
FEND = -L/usr/local/opt/openblas/lib -lopenblas 
#FEND = -lopenblas

ifeq ($(OMP),ON)
FFLAGS = -O2 -c -w --openmp -std=legacy
FLINK = gfortran -w -o $(EXEC) --openmp

endif

SRC = ../src

.PHONY: all clean list

OBJECTS =  dielectric_square.o \
  $(SRC)/prini.o \
  $(SRC)/legeexps.o \
  $(SRC)/pplot2.o \
  $(SRC)/hkrand.o \
  $(SRC)/hank103.o \
  $(SRC)/dlaran.o \
  $(SRC)/quads.o \
  $(SRC)/chunks.o \
  $(SRC)/helm_kernels.o \
  $(SRC)/get_mat_guru.o \
  $(SRC)/get_mat_guru_nreg.o \
  $(SRC)/helm_mats.o \
  $(SRC)/gmres_solvers.o \
  $(SRC)/levrtree2d.o \
  $(SRC)/nearquadrouts.o \


#
# use only the file part of the filename, then manually specify
# the build location
#

%.o : %.f %h
	$(FC) $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)
	./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



