EXEC = dirichletsolver

FC = gfortran
FFLAGS = -O2 -w -std=legacy -fPIC 
FLINK = gfortran -w -o $(EXEC) 
FEND = -L/usr/local/opt/openblas/lib -lopenblas 
#FEND = -lopenblas

FOMP = --openmp
FENDOMP = -lgomp


SRC = ../src

.PHONY: all lib clean list

OBJECTS =  prini_new.o \
  legeexps.o \
  pplot2.o \
  hkrand.o \
  hank103.o \
  dlaran.o \
  quads.o \
  chunks.o \
  chunks_fmm2dbie.o \
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
	$(FC) -c $(FFLAGS) $< -o $@

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(EXEC)
	$(FLINK) $(OBJECTS) $(FEND)

lib: $(OBJECTS)
	$(FC) -shared -fPIC $(OBJECTS) -o libdirdatgen.so $(FEND)

example: lib 
	$(FC) $(FFLAGS) $(FOMP) ext_dir_daria_datgen.f -o $(EXEC) $(FEND) -L. -ldirdatgen $(FENDOMP)
	./$(EXEC)

clean:
	rm -f $(OBJECTS)
	rm -f $(EXEC)

list: $(SOURCES)
	$(warning Requires:  $^)



