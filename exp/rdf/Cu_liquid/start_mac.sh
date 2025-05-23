#!/bin/bash
../../../bin/AAVDP_mac --rdf ./Cu_liquid.lammps -r 8 -n 160 -o ./Cu_liquid.rdf
../../../bin/AAVDP_mac --rdf ./Cu_liquid.lammps -r 3 -n 60 -o ./Cu_liquid_coord.rdf -coord -coord_o ./Cu_liquid_coord.lmc