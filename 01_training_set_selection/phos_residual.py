#!/usr/bin/env python

import datetime
import numpy
import random
import sys

################################################################

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_"]
p_class = ["S-", "T-", "Y-", "SP", "TP"]

################################################################

################################################################

inter_raw = []
screen_raw = []

with open("whole_seq/phos_human_"+sys.argv[1]+".txt") as base_phos:
    for line in base_phos:
    	inter_raw.append(line)

with open("screened/phos_screened_"+sys.argv[1]+".txt") as base_nonphos:
    for line in base_nonphos:
        screen_raw.append(line)

################################################################

sel_file = open("screened/phos_residual_"+sys.argv[1]+".txt", 'w')
sel_list = []
sel_cut = 0.0

for i_1 in range(len(inter_raw)):
	det_1 = 0
	for i_2 in range(len(screen_raw)):
		if inter_raw[i_1] == screen_raw[i_2]:
			det_1 += 1
			break
	if det_1 < 1:
		sel_list.append(i_1)
		sel_file.write(inter_raw[i_1])


	if i_1 % 100 == 99:
		print p_class[int(sys.argv[1])], i_1, len(sel_list)

print len(inter_raw), len(screen_raw), len(sel_list)
sel_file.close()
