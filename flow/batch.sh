#!/usr/bin/bash

parallel --gnu --progress zcat /var/phy/project/nukeserv/jsm55/eos-study/HotQCD-errors/2015-03-18/HotQCD-EOS-{}/\*1_block.gz \| ./urqmd2flow.py \| ./flow.e \> results/flow_{}.dat ::: {0..99}
