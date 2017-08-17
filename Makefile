########## User-definable stuff ##########
#
###Compiler and compilation options
COMP_CC = icc
OPTIONS = -Wall -O3 -std=gnu99
#
###Path to libraries and headers
###If two or more of the dependencies reside in the same paths, only
###one instance is necessary.
#cfitsio
FITS_INC = -I/users/damonge/include
FITS_LIB = -L/users/damonge/lib
#healpix
HPIX_INC =
HPIX_LIB =
#
########## End of user-definable ##########

DEFINEFLAGS += -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF 

OPTIONS += -fopenmp
DEFINEFLAGS += -D_HAVE_OMP

OPTIONS += $(DEFINEFLAGS)

INC_ALL = -I./src $(FITS_INC) $(HPIX_INC)
LIB_ALL = $(FITS_LIB) $(HPIX_LIB) -lcfitsio -lchealpix
LIB_ALL += -lm

COMMONO = src/common.o
HPIXO = src/healpix_extra.o
MAIN = src/main.c
OFILES = $(COMMONO) $(HPIXO)

EXEC = FieldXCorr

default : $(EXEC)

%.o : %.c
	$(COMP_CC) $(OPTIONS) $(INC_ALL) -c $< -o $@

$(EXEC) : $(OFILES)
	$(COMP_CC) $(OPTIONS) $(INC_ALL) $(OFILES) $(MAIN) -o $(EXEC) $(LIB_ALL)

clean :
	rm -f src/*.o

cleaner :
	rm -f *~ src/*.o src/*~  $(EXEC)
