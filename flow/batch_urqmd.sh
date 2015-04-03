#!/usr/bin/bash
path=$1
parallel --gnu --progress zcat $path/HotQCD-EOS-{}/\*.gz \| ./flow.e \> results/v3/flow_{}.dat ::: {0..19}
