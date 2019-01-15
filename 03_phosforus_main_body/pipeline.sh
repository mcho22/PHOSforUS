#!/usr/bin/env bash

# Validation analysis with residual datasets (see folder 01)
# Would result in similar values listed in log_pfos_alt_res.txt
# (Random noise may occur)
python binary_index.py 0
python binary_index.py 1
python binary_index.py 2
python binary_index.py 3
python binary_index.py 4


