#!/bin/bash
../../../bin/AAVDP_mac --xrd ./ZnO.vasp -2t 10 80 -o ./ZnO_line.xrd
../../../bin/AAVDP_mac --xrd ./ZnO.vasp -2t 10 80 -o ./ZnO.xrd -scherrer -scherrer_m 0.5 -scherrer_d 294 -scherrer_d2t 0.02