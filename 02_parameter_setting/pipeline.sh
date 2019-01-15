#!/usr/bin/env bash

#### Training / testing set preparation pipeline
#### Phosphorylation sites -> screened + residual
#### Non-phosphorylated seq -> (phos-similar) + screened + residual

# 100x resampling analysis (with binary feature / pssm analysis)
# Should result in values listed in sample_binary_index_testing.txt
# (Actual values could be marginally different because of random sampling)
python binary_index.py 0
python binary_index.py 1
python binary_index.py 2
python binary_index.py 3
python binary_index.py 4

# Preset parameter calculation
# Should result in parameter sets found in 03~ folder / PHOSforUS folder
python parameter_pickler.py 0
python parameter_pickler.py 1
python parameter_pickler.py 2
python parameter_pickler.py 3
python parameter_pickler.py 4
