
EXEC = int3

FC = gfortran
FFLAGS = -O2 -c -w -std=legacy 
FLINK = gfortran -w -o $(EXEC) 
FEND = -L/usr/local/opt/openblas/lib -lopenblas 
#FEND = -lopenblas

ifeq ($(OMP),ON)
FFLAGS = -O2 -c -w --openmp
FLINK = gfortran -w -o $(EXEC) --openmp

endif

SRC = ../src

EFOL = ../data

FL = ../../../forleslie


.PHONY: all clean list

OBJECTS =  ext_dir_solver.o \
  $(SRC)/prini.o \
  $(SRC)/legeexps.o \
  $(SRC)/pplot2.o \
  $(SRC)/hkrand.o \
  $(SRC)/hank103.o \
  $(SRC)/dlaran.o \
  $(SRC)/quads.o \
  $(SRC)/chunks.o \
  $(SRC)/helm_kernels.o \
  $(SRC)/get_mat.o \
  $(SRC)/gmres_solvers.o \


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
	mkdir -p $(EFOL)
	$(FLINK) $(OBJECTS) $(FEND)
	mv $(EXEC) ./$(EFOL)/
	cd $(EFOL) && ./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



