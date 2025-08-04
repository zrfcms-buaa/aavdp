#!/bin/bash
../../../bin/AAVDP_mac --ked ../Cu.lmp -e Cu -q 0.8 -z 0 0 1 -o ./Cu.ked -gauss -gauss_dx 0.00375 -gauss_sig 0.0075
../../../bin/AAVDP_mac --ked ../Cu_vacancies.lmp -e Cu -q 0.8 -z 0 0 1 -o ./Cu_vacancies.ked -gauss -gauss_dx 0.00375 -gauss_sig 0.0075