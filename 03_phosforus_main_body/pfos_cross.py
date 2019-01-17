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
non_proline = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "Q", "R", "S", "T", "V", "W", "Y"]
p_class = ["S-", "T-", "Y-", "SP", "TP"]

weight_types = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
es_presel_1 = [2, 6, 9, 12]
es_presel_2 = [[0, 2, 4, 2], [8, 10, 12, 4]]

silent_flag = 0

#### Temporary functions

#### Test sequences import

iterat_cycle = 10
const_length = 100
p_tseq = []
n_tseq = []
p_vseq = []
n_vseq = []

with open("training_dataset/screened_p/phos_screened_"+sys.argv[1]+".txt") as p_file:
    for line in p_file:
        sst = line.rstrip('\n').split('\t')[3][0:29]
        p_tseq.append(sst)

with open("training_dataset/screened_np/nonphos_screened_"+sys.argv[1]+".txt") as n_file:
    for line in n_file:
        sst = line.rstrip('\n').split('\t')[3][0:29]
        n_tseq.append(sst)

with open("testing_dataset/phos_residual_"+sys.argv[1]+".txt") as p_file:
    for line in p_file:
        sst = line.rstrip('\n').split('\t')[3][0:29]
        p_vseq.append(sst)

with open("testing_dataset/nonphos_residual_"+sys.argv[1]+".txt") as n_file:
    for line in n_file:
        sst = line.rstrip('\n').split('\t')[3][0:29]
        n_vseq.append(sst)

#### File import

index_label, index_list = input_output_handle.input_indices()
esc_label, esc_list = input_output_handle.input_escape()
model_list = input_output_handle.preset_clf_import()

#### Analysis & export

tv_compile = []

####

for i_0 in range(iterat_cycle):
    tar_files = [["Training_construct", ""], ["Validating_construct", ""]]
    p_tsample = random.sample(p_tseq, const_length)
    n_tsample = random.sample(n_tseq, const_length)
    p_vsample = random.sample(p_vseq, const_length)
    n_vsample = random.sample(n_vseq, const_length)
    
    for i_0_1 in range(const_length):
        tar_files[0][1] += p_tsample[i_0_1]
        tar_files[0][1] += n_tsample[i_0_1]
        tar_files[1][1] += p_vsample[i_0_1]
        tar_files[1][1] += n_vsample[i_0_1]
    
    #print len(tar_files[0][1]), len(tar_files[1][1])
    #print tar_files[0][1][0:29]
    
    for i_2 in range(len(tar_files)):
        tar_code = sequence_calc.seq_codifier(tar_files[i_2][1])
        tar_sites = sequence_calc.site_identifier(tar_code)
        tv_scores = []
        tar_ds = 0
        for i_1 in range(5):
            tar_vex = []
            for i_2_1 in range(len(tar_sites)):
                if tar_sites[i_2_1][0] % 29 == 28:
                    tar_vex.append([tar_sites[i_2_1][0], i_1])
                    tar_ds += i_1
            #print i_1, len(tar_vex)
            tar_values = sequence_calc.index_scorer(tar_code, index_list, weight_types) + sequence_calc.escape_scorer(tar_code, esc_list, es_presel_1, es_presel_2)
            tar_res = sequence_calc.phosforus_scorer(tar_vex, tar_values, model_list)
            tv_slist = []
            tv_stats = []
            for i_3 in range(len(tar_res)):
                if tar_res[i_3][0] % 58 == 15:
                    tv_slist.append(tar_res[i_3][3])
                    tv_stats.append(1)
                elif tar_res[i_3][0] % 58 == 44:
                    tv_slist.append(tar_res[i_3][3])
                    tv_stats.append(0)
            tv_scores.append(metrics.roc_auc_score(tv_stats, tv_slist))
        
        tv_compile.append(tv_scores)
        if i_0 % 5 == 4 and i_2 == 1:
            tv_csub = numpy.asarray(tv_compile)
            print "#### Cycle_"+str(i_0+1)+":", p_class[int(sys.argv[1])], numpy.mean(tv_csub[-10:, 0]), numpy.mean(tv_csub[-10:, 1]), numpy.mean(tv_csub[-10:, 2]), numpy.mean(tv_csub[-10:, 3]), numpy.mean(tv_csub[-10:, 4])

tv_carray = numpy.asarray(tv_compile)
print "Used sequence:", p_class[int(sys.argv[1])]
print "Calculation result:", numpy.mean(tv_carray[:, int(sys.argv[1])])
print "Applied parameters // S-      T-      Y-      SP      TP"
print numpy.mean(tv_carray[:, 0]), numpy.mean(tv_carray[:, 1]), numpy.mean(tv_carray[:, 2]), numpy.mean(tv_carray[:, 3]), numpy.mean(tv_carray[:, 4])

####	
