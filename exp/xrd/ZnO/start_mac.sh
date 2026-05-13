#!/bin/bash
../../../bin/AAVDP_mac --xrd ./ZnO.vasp -2t 10 80 -o ./ZnO_line.xrd
../../../bin/AAVDP_mac --xrd ./ZnO.vasp -2t 10 80 -o ./ZnO.xrd -pseudo -pseudo_uvw 0.07313 0.0 0.07313 -pseudo_na 0.5 -pseudo_nb 0.0 -pseudo_d2t 0.02