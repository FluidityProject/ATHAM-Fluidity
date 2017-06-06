include envcheck.mk

preprocess: envcheck
	@echo **********Generating the mesh using gmsh in 3d:
	gmsh -3 -optimize -bin -o structured_thin.msh src/structured_thin.geo

run: envcheck
	@echo **********Calling flredecomp in parallel with verbose log output enabled:
	mpiexec -n 24 ${BINPREFIX}flredecomp -i 1 -o 24 -v -l precip_extreme precip_extreme_parallel
	@echo **********Calling fluidity in parallel with verbose log output enabled:
	mpiexec -n 24 ${BINPREFIX}fluidity -v2 -l precip_extreme_parallel.flml

input: clean preprocess

clean:
	@echo **********Cleaning the output from previous fluidity runs:
	rm -rf *_? *_?? *_??? *stat *msh *halo *vtu fluidity* flredecomp* matrixdump* *_parallel.flml *checkpoint*
