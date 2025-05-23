#!/bin/bash
../../../bin/AAVDP_mac --ned ./FeCo.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -o ./FeCo.ned
../../../bin/AAVDP_mac --ned ./FeCo_random.lmp -e Fe Co -dw 0.5500 0.5500 -2t 0 150 -l 1.5400 -c 10 10 10 -o ./FeCo_random.ned