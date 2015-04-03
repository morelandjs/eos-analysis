#!/usr/bin/bash
path=$1
parallel --gnu --progress ./iss2flow.py $path/HotQCD-EOS-{}/\*.gz \> results/iss/flow_{}.dat ::: {0..19}
