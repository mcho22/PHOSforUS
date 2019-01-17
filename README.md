[PHOSforUS predictor usage guide]

[1] PHOSforUS features

PHOSforUS is a fully biophysical parameter-based protein phosphorylation site predictor.
From input sequences, it calculates 18 biophysical features, including:
     - Horizontal features (6): Hydrophobicity, Secondary structure propensity, PII propensity, etc.
     - Vertical features (4): AA molecular weight, residue charge, surface area, etc.
     - eSCAPE thermodynamic features (8): dG, dHap, dHpol, TdS, etc.
 All individual biophysical features are sent to sub-predictors, based on Naive Bayes classifier.
 Resulting scores are sent to downstream metapredictor, based on gradient boosting classifier.
 Prediction results are compiled into single result file per each sequence, which could be found in result_output folder.

[2] System requirement

 PHOSforUS is built on python 2.7 environment.
 Along with base packages, it requires installation of numpy and scikit-learn package.
 As of now, no other preparations are needed.

[3] Sample commands

 Option #1: python phosforus.py
 
 If no other arguments are given, PHOSforUS automatically identifies sequence files in input_seq folder.
 Sequence files should follow fasta format for smooth analysis.

 Option #2: python phosforus.py [filename]
 
 If file name is explicitly given, PHOSforUS runs for that specific file only.

 Option #3: python phosforus.py -m [sequence]
 
 If "-m" argument is given, PHOSforUS takes manually typed sequence argument as its input.

[4] Optional arguments

 "-m": manual sequence input option. Should be followed by sequence input.
 "-v": verbose option. All site-wise prediction results will be printed in command interface.
