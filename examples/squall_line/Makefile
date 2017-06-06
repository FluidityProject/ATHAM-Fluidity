include envcheck.mk

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 3d:
	gmsh -2 -optimize -bin -o vbox.msh src/vbox.geo

run: envcheck
	@echo **********Calling flredecomp in parallel with verbose log output enabled:
	mpiexec -n 24 ${BINPREFIX}flredecomp -i 1 -o 24 -v -l squall_line squall_line_parallel
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n 24 ${BINPREFIX}fluidity -v2 -l squall_line_parallel.flml

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *_? *_?? *_??? *stat *msh *Mesh *halo *vtu fluidity* flredecomp* matrixdump* *_parallel.flml *checkpoint*
