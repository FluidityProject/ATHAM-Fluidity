PROJECT = single_mountain

input: clean
	gmsh -2 -o vbox.msh src/vbox.geo
        cp -f stored_output/stored_$(PROJECT).stat .

clean:
	rm -f $(PROJECT).stat
	rm -rf $(PROJECT)_?
	rm -rf $(PROJECT)_??
	rm -f *vtu
	rm -f fluidity*
	rm -f *.msh *.halo
	rm -f *_parallel.flml
	rm -f matrixdump matrixdump.info
