#!/usr/bin/env python

import datetime
import numpy
import random
import sys

################################################################

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_"]
p_class = ["S-", "T-", "Y-", "SP", "TP"]

################################################################

def seq_id(si_seq1, si_seq2):
	si_ids = 0
	for si_1 in range(len(si_seq1)):
		if si_seq1[si_1] == si_seq2[si_1]:
			si_ids += 1
	return si_ids

################################################################

phos_raw = []
phos_list = []

with open("whole_seq/phos_human_"+sys.argv[1]+".txt") as base_phos:
    for line in base_phos:
        ssl = line.rstrip('\n').split('\t')[3]
        phos_list.append(ssl)
        phos_raw.append(line)

################################################################

sel_file = open("screened/phos_screened_"+sys.argv[1]+".txt", 'w')
sel_list = []

for i_1 in range(len(phos_list)):
	seq_1 = phos_list[i_1][4:14] + phos_list[i_1][16:25]
	score_list = [0]
	for i_2 in range(i_1):
		if i_1 != i_2:
			seq_2 = phos_list[i_2][4:14] + phos_list[i_2][16:25]
			score_sub = seq_id(seq_1, seq_2)
			score_list.append(score_sub)
			if score_sub >= 7:
				break
	if numpy.amax(score_list) < 7:
		sel_file.write(phos_raw[i_1])
		sel_list.append(i_1)
	if i_1 % 100 == 99:
		print p_class[int(sys.argv[1])], i_1, len(sel_list)

