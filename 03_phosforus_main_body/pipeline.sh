#!/usr/bin/env bash

# Validation analysis with residual datasets (see folder 01)
# Cross-test using parameter application change included
python pfos_cross.py 0
python pfos_cross.py 1
python pfos_cross.py 2
python pfos_cross.py 3
python pfos_cross.py 4

# Please check printed outputs
