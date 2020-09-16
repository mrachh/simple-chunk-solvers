
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

OBJECTS =  dielectric_gen.o \
  $(SRC)/projectonpolycheb.o \
  $(SRC)/qrsolve.o \
  $(SRC)/chebexps.o \
  $(SRC)/prini.o \
  $(SRC)/legeexps.o \
  $(SRC)/pplot2.o \
  $(SRC)/hkrand.o \
  $(SRC)/evalpotvol_chunk.o \
  $(SRC)/getdomains.o \
  $(SRC)/hank103.o \
  $(SRC)/dlaran.o \
  $(SRC)/quads.o \
  $(SRC)/chunks.o \
  $(SRC)/chunkfunc_withk.o \
  $(SRC)/helm_kernels.o \
  $(SRC)/get_mat_guru.o \
  $(SRC)/get_mat_guru_nreg2.o \
  $(SRC)/helm_mats.o \
  $(SRC)/gmres_solvers.o \
  $(SRC)/levrtree2d.o \
  $(SRC)/nearquadrouts.o \
  $(SRC)/cdjseval2d.o \
  $(SRC)/dfft.o \
  $(SRC)/fmmcommon2d.o \
  $(SRC)/get_fmm_thresh.o \
  $(SRC)/h2dcommon.o \
  $(SRC)/h2dterms.o \
  $(SRC)/helm_comb_eval.o \
  $(SRC)/helm_dielec_eval.o \
  $(SRC)/helmrouts2d.o \
  $(SRC)/helmrouts2d_dir.o \
  $(SRC)/hfmm2dmain.o \
  $(SRC)/hfmm2dpart.o \
  $(SRC)/next235.o \
  $(SRC)/wideband2d.o \

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



