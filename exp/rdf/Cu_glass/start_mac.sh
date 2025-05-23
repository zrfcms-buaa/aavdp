#!/bin/bash
../../../bin/AAVDP_mac --rdf ./Cu_glass.lammps -r 8 -n 160 -o ./Cu_glass.rdf
../../../bin/AAVDP_mac --rdf ./Cu_glass.lammps -r 3 -n 60 -o ./Cu_glass_coord.rdf -coord -coord_o ./Cu_glass_coord.lmc