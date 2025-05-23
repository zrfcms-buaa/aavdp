#!/bin/bash
../../../bin/AAVDP_mac --dkd ./hcp.vasp -q 1.5 -px 600 -py 600 -o ./hcp.dkd -monte -monte_seed ../../../RandomSeeds.data -monte_o ./hcp.mc
../../../bin/AAVDP_mac --dkd ./hcp.vasp -q 1.5 -px 600 -py 600 -o ./hcp.no_mc.dkd