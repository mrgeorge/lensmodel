#
#  Makefile for the Adiabatic Contraction code
#    Written by:	Oleg Gnedin
#    Last modified:	August 17, 2005
#
#-------------------  macro definitions  -------------------------------
FC       = gfortran
CC       = gcc
OPTS     = -O
FFLAGS   = -O -Wall
CFLAGS   = -O -Wall
LIB      = -lm
LDFLAGS  = -lcc_dynamic
#---------------------  implicit rules  --------------------------------
%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

OBJ = contra.o spline.o los.o

contra: $(OBJ) makefile
	$(FC) $(OPTS) -o $@ $(OBJ) $(LIB)

$(OBJ): makefile


# ----- contra_aw -----
OBJ2 = contra_aw.o spline.o los.o

contra_aw: $(OBJ2) makefile
	$(FC) $(OPTS) -o $@ $(OBJ2) $(LIB)

$(OBJ2): makefile
