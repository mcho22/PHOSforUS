#!/usr/bin/env bash

#### Training / testing set preparation pipeline
#### Phosphorylation sites -> screened + residual
#### Non-phosphorylated seq -> (phos-similar) + screened + residual

# 30% identity screening of phosphorylation sites
python phos_screen.py 0
python phos_screen.py 1
python phos_screen.py 2
python phos_screen.py 3
python phos_screen.py 4

# Extract residual non-phosphorylated sequences (intermediate - screened = residual)
python phos_residual.py 0
python phos_residual.py 1
python phos_residual.py 2
python phos_residual.py 3
python phos_residual.py 4

# 30% identity screening of non-phosphorylated sequences against phosphorylation sites
python nonphos_pre_screen.py 0
python nonphos_pre_screen.py 1
python nonphos_pre_screen.py 2
python nonphos_pre_screen.py 3
python nonphos_pre_screen.py 4

# 30% identity screening of non-phosphorylated sequences (intermediate -> screened)
python nonphos_post_screen.py 0
python nonphos_post_screen.py 1
python nonphos_post_screen.py 2
python nonphos_post_screen.py 3
python nonphos_post_screen.py 4

# Extract residual non-phosphorylated sequences (intermediate - screened = residual)
python nonphos_residual.py 0
python nonphos_residual.py 1
python nonphos_residual.py 2
python nonphos_residual.py 3
python nonphos_residual.py 4
