EXEC = int2

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

OBJECTS =  ext_dir_daria.o daria_getspecs.o \
  prini.o \
  legeexps.o \
  pplot2.o \
  hkrand.o \
  hank103.o \
  dlaran.o \
  quads.o \
  chunks.o \
  chunkfunc_withk.o \
  getdomains_sm2.o \
  lap_dlp_eval.o \
  helm_kernels.o \
  get_mat_guru.o \
  get_mat_guru_nreg.o \
  helm_mats.o \
  gmres_solvers.o \
  levrtree2d.o \
  nearquadrouts.o \
  cdjseval2d.o \
  dfft.o \
  fmmcommon2d.o \
  get_fmm_thresh.o \
  h2dcommon.o \
  h2dterms.o \
  l2dterms.o \
  helm_comb_eval.o \
  helmrouts2d.o \
  helmrouts2d_dir.o \
  hfmm2dmain.o \
  hfmm2dpart.o \
  next235.o \
  wideband2d.o \
  laprouts2d.o \
  laprouts2d_dir.o \
  lfmm2dmain.o \
  lfmm2dpart.o \

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

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



