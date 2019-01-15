#!/usr/bin/env python

#### Descriptions

#### Generic module import
import datetime
import numpy
import os
import pickle
import random
import sys

from sklearn import metrics
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB

#### Custom module import

from custom_module import input_output_handle, sequence_calc, misc_functions

#### Universal variables

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_"]
p_class = ["S-", "T-", "Y-", "SP", "TP"]

weight_types = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
es_presel_1 = [2, 6, 9, 12]
es_presel_2 = [[0, 2, 4, 2], [8, 10, 12, 4]]

silent_flag = 0

#### Temporary functions

#### Test sequences import

iterat_cycle = 100
const_length = 100
p_tseq = []
n_tseq = []
p_vseq = []
n_vseq = []

with open("training_dataset/screened_p/phos_screened_"+sys.argv[1]+".txt") as p_file:
    for line in p_file:
        sst = line.rstrip('\n').split('\t')[3]
        p_tseq.append(sst)

with open("training_dataset/screened_np/nonphos_screened_"+sys.argv[1]+".txt") as n_file:
    for line in n_file:
        sst = line.rstrip('\n').split('\t')[3]
        n_tseq.append(sst)

with open("testing_dataset/phos_residual_"+sys.argv[1]+".txt") as p_file:
    for line in p_file:
        sst = line.rstrip('\n').split('\t')[3]
        p_vseq.append(sst)

with open("testing_dataset/nonphos_residual_"+sys.argv[1]+".txt") as n_file:
    for line in n_file:
        sst = line.rstrip('\n').split('\t')[3]
        n_vseq.append(sst)

#### File import

index_label, index_list = input_output_handle.input_indices()
esc_label, esc_list = input_output_handle.input_escape()
model_list = input_output_handle.preset_clf_import()

#### Analysis & export

tv_compile = []

####

for i_0 in range(iterat_cycle):
	tv_scores = []
	tar_files = [["Training_construct", ""], ["Validating_construct", ""]]
	p_tsample = random.sample(p_tseq, const_length)
	n_tsample = random.sample(n_tseq, const_length)
	p_vsample = random.sample(p_vseq, const_length)
	n_vsample = random.sample(n_vseq, const_length)
	for i_00 in range(const_length):
		tar_files[0][1] += p_tsample[i_00]
		tar_files[0][1] += n_tsample[i_00]
		tar_files[1][1] += p_vsample[i_00]
		tar_files[1][1] += n_vsample[i_00]

	for i_1 in range(len(tar_files)):
		tar_code = sequence_calc.seq_codifier(tar_files[i_1][1])
		tar_sites = sequence_calc.site_identifier(tar_code)
		tar_values = sequence_calc.index_scorer(tar_code, index_list, weight_types) + sequence_calc.escape_scorer(tar_code, esc_list, es_presel_1, es_presel_2)
		tar_res = sequence_calc.phosforus_scorer(tar_sites, tar_values, model_list)
		tv_scoresub = 0.0
		tv_slist = []
		tv_stats = []
		for i_2 in range(len(tar_res)):
			if tar_res[i_2][0] % 58 == 15:
				tv_slist.append(tar_res[i_2][3])
				tv_stats.append(1)
				if tar_res[i_2][3] > 0:
					tv_scoresub += 0.5/float(const_length)
			elif tar_res[i_2][0] % 58 == 44:
				tv_slist.append(tar_res[i_2][3])
				tv_stats.append(0)
				if tar_res[i_2][3] < 0:
					tv_scoresub += 0.5/float(const_length)
		tv_scores.append(tv_scoresub)
		tv_scores.append(metrics.roc_auc_score(tv_stats, tv_slist))
	print "#### Cycle_"+str(i_0+1)+":", tv_scores[0], tv_scores[2], tv_scores[1], tv_scores[3]
	tv_compile.append(tv_scores)

tv_carray = numpy.asarray(tv_compile)
print numpy.mean(tv_carray[:, 0]), numpy.mean(tv_carray[:, 2]), numpy.mean(tv_carray[:, 1]), numpy.mean(tv_carray[:, 3])

####	
