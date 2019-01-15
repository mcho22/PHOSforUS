#!/usr/bin/env python

#### Descriptions

#### Generic module import
import datetime
import numpy
import os
import pickle
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

#### File import

if len(sys.argv) >= 2 and sys.argv[1] == "-h":
	misc_functions.help_print()
	sys.exit()

if "-v" in sys.argv:
	silent_flag += 1
	sys.argv.remove("-v")

tar_files = input_output_handle.input_fileset(sys.argv)
index_label, index_list = input_output_handle.input_indices()
esc_label, esc_list = input_output_handle.input_escape()
model_list = input_output_handle.preset_clf_import()

#### Analysis & export

for i_1 in range(len(tar_files)):
	start_time = datetime.datetime.now()
	tar_code = sequence_calc.seq_codifier(tar_files[i_1][1])
	tar_sites = sequence_calc.site_identifier(tar_code)
	tar_values = sequence_calc.index_scorer(tar_code, index_list, weight_types) + sequence_calc.escape_scorer(tar_code, esc_list, es_presel_1, es_presel_2)
	tar_res = sequence_calc.phosforus_scorer(tar_sites, tar_values, model_list)

	for i_2 in range(len(tar_res)):
		if silent_flag > 0:
			if i_2 == 0:
				print ""
			print tar_res[i_2][0], '\t', tar_res[i_2][1], '\t', tar_res[i_2][2], '\t', tar_res[i_2][3], '\t', tar_res[i_2][4]

	end_time = datetime.datetime.now()
	input_output_handle.output_formatter(tar_files[i_1], tar_res, i_1, start_time, end_time)
	print "## Entry "+str(i_1+1)+" calculation completed: "+tar_files[i_1][0][:min(25, len(tar_files[i_1][0]))]

####	
