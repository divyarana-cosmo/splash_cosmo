FC = gfortran
FFLAGS = -O

#FC = g77
#FFLAGS = -O -static

#FC = f95
#FFLAGS = -O

#FC = ifort
#FFLAGS = -w95 -O -i_dynamic


all:
	make mandc.x

mandc.x: mandc.f makefile
	$(FC) -o mandc.x $(FFLAGS) mandc.f 

 
