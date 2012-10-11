# cfitsio source files, libraries, and flags
# modify for your computer
prefix=/usr/local/bin
exec_prefix=${prefix}
libdir=${exec_prefix}/lib
includedir=${prefix}/include

#cfitsio libraries and flags
CFIOlibs = -L${libdir} -lcfitsio
CFIOflags = -I${includedir}

# end cfitsio modifications
       
# directory of dophot source files, structs, and headers
# modify for your directory structure as needed
odir = ./

#dophot libraries and flags
DOPHflags =  -I${odir}structs -I${odir}headers 

CC= gcc
CFLAGS= -g -fno-inline -Wall ${DOPHflags} ${CFIOflags} 
#CFLAGS= -O2 ${DOPHflags} ${CFIOflags} 
LIBS= ${CFIOlibs} 
# optimize with flag: -O2 if things compile and run nicely
# debug flag: -Wall
# Valgrind flag: -g -fno-inline
# run: valgrind --dsymutil=yes --leak-check=full ./dophot

DoPhotCOBJECTS = $(odir)mini_mathlib.o $(odir)cast_arr.o \
	$(odir)write_to_log.o $(odir)diag_mmult.o\
	$(odir)parupd.o $(odir)errupd.o $(odir)covarupd.o \
        $(odir)lu_comp.o \
        $(odir)tagi4.o $(odir)medfil.o \
        $(odir)makenoise.o $(odir)paravg.o \
        $(odir)ellipse.o $(odir)ellipint.o \
        $(odir)stdinpt.o $(odir)suminpt.o \
	$(odir)stdotpt.o $(odir)sumout.o \
        $(odir)shdout.o  $(odir)covarout.o \
        $(odir)badotpt.o \
        $(odir)newfits.o $(odir)readfits.o \
        $(odir)elarea.o $(odir)impaper.o \
        $(odir)toofaint.o $(odir)toobright.o $(odir)offpic.o \
        $(odir)addlims.o $(odir)oblims.o \
        $(odir)skyfun_plane.o $(odir)skyfun_hub.o \
        $(odir)varipar.o \
        $(odir)parinterp.o $(odir)galaxy.o \
        $(odir)transmask.o $(odir)makemask.o \
        $(odir)guess.o \
	$(odir)pgauss.o $(odir)gauss.o \
	$(odir)extpgauss.o \
        $(odir)empiricals.o \
	$(odir)onefit.o $(odir)twofit.o \
        $(odir)fillerup.o \
        $(odir)oblit.o $(odir)cosmic.o \
	$(odir)probgal.o $(odir)chisq.o \
        $(odir)addstar.o stampout.o \
        $(odir)add_analytic_or_empirical_obj.o \
        $(odir)shape.o \
	$(odir)isearch.o \
        $(odir)bestab.o \
        $(odir)improve.o \
	$(odir)warmstart.o \
	$(odir)paramfile.o \
	$(odir)io_dophotc.o \
	$(odir)tuneup.o \
        $(odir)dophot.o 

#	$(odir)sersic.o \

dophot: $(DoPhotCOBJECTS)
	     $(CC) $(DoPhotCOBJECTS) \
	     -o ./working_data/dophot $(CFLAGS) $(LIBS) -lm

test: $(DoPhotCOBJECTS)
	     $(CC) $(DoPhotCOBJECTS) \
	     -o ./verif_data/dophot $(CFLAGS) $(LIBS) -lm

clean: 
	rm -f *.o ./working_data/dophot ./verif_data/dophot
