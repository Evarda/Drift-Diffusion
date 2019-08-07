# Start of the makefile
# Defining variables
objects = driftDiffusion.o globalConstants.o generateMesh.o LUDecomp.o
f90compiler = gfortran
debugOp = -fcheck=all -Wall

# Makefile

driftDiffusion:	driftDiffusion.o	$(objects)
	$(f90compiler) -o driftDiffusion $(objects)

globalConstants.mod:	globalConstants.f90
	$(f90compiler) -c -g $(debugOp) globalConstants.f90

globalConstants.o:	globalConstants.f90
	$(f90compiler) -c -g $(debugOp) globalConstants.f90

generateMesh.o:	generateMesh.f90 globalConstants.o
	$(f90compiler) -c -g $(debugOp) generateMesh.f90

LUDecomp.o:	LUDecomp.f90 globalConstants.o
	$(f90compiler) -c -g $(debugOp) LUDecomp.f90

driftDiffusion.o:	driftDiffusion.f90	globalConstants.mod
	$(f90compiler) -c -g $(debugOp) driftDiffusion.f90

# Cleaning everything
clean:
	rm globalConstants.mod
	rm $(objects)
	rm driftDiffusion
# End of the makefile