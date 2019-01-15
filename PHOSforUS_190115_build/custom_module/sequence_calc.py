#!/usr/bin/env python

#### Generic module import

import numpy

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.naive_bayes import GaussianNB

#### Module variables

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_"]
p_class = ["S-", "T-", "Y-", "SP", "TP"]
es_prelabel = ['dHapN', 'TdSN', 'ddGD', 'dHpD', 'dGN-dGD', 'dHapN-dHapD', 'dHpN-dHpD', 'dHapN-dHpN']

escape_params = [[0.8195, 0.8195, 0.7665, 0.7665, 0.7791, 0.7791, 0.7047, 0.7047, 0.7023, 0.7023, 0.7801, 0.7801, 0.7115, 0.7115, 0.7319, 0.7319], 
[0.7492, 0.7492, 0.7632, 0.7632, 0.7524, 0.7524, 0.7507, 0.7507, 0.7479, 0.7479, 0.6708, 0.6708, 0.8047, 0.8047, 0.6782, 0.6782], 
[4696.0, 4696.0, -5068.0, -5068.0, 6195.0, 6195.0, 1998.0, 1998.0, -4154.0, -4154.0, 447.53, 447.53, -94.37, -94.37, 3872.0, 3872.0]]

#### Module functions

def seq_codifier(sc_seq):
    sc_list = []
    for sc_0 in range(14):
    	sc_list.append(20)
    for sc_1 in range(len(sc_seq)):
        sc_det_1 = 0
        for sc_2 in range(21):
            if sc_seq[sc_1] == aa_list[sc_2]:
                sc_list.append(sc_2)
                sc_det_1 += 1
        if sc_det_1 < 1:
            sc_list.append(20)
            sc_det_1 += 1
    for sc_3 in range(14):
    	sc_list.append(20)
    return sc_list

def site_identifier(sc_code):
	si_list = []
	for si_1 in range(len(sc_code)-1):
		if sc_code[si_1] == 15:
			if sc_code[si_1+1] == 12:
				si_list.append([si_1, 3])
			else:
				si_list.append([si_1, 0])
		elif sc_code[si_1] == 16:
			if sc_code[si_1+1] == 12:
				si_list.append([si_1, 4])
			else:
				si_list.append([si_1, 1])
		elif sc_code[si_1] == 19:
			si_list.append([si_1, 2])
	return si_list

def index_scorer(is_code, is_indices, is_weights):
    is_raw = []
    for is_1 in range(len(is_code)):
        is_sub = []
        for is_2 in range(len(is_indices)):
            is_sub.append(is_indices[is_2][is_code[is_1]])
        is_raw.append(is_sub)
    is_trans = numpy.transpose(is_raw)
    is_vlist = []
    for is_3 in range(len(is_indices)):
        if is_weights[is_3] == 0:
            is_vlist.extend(is_trans[is_3][4:-4])
        else:
            for is_4 in range(len(is_code)-8):
                is_vsub = numpy.vdot(is_trans[is_3][(is_4):(9+is_4)], [1.0,2.0,3.0,4.0,5.0,4.0,3.0,2.0,1.0])
                is_vlist.append(is_vsub*0.04)
    return is_vlist

def escape_scorer(es_code, es_indices, es_sel1, es_sel2):
    es_slen = len(es_code)-8
    es_raw = []
    for es_1 in range(len(es_code)-2):
        es_sub = []
        es_key = es_code[es_1]*441 + es_code[es_1+1]*21 + es_code[es_1+2]
        for es_2 in range(len(es_indices)):
            es_sub.append(es_indices[es_2][es_key])
        es_raw.append(es_sub)
    es_trans = numpy.transpose(es_raw)
    es_vlist = []
    for es_3 in range(len(es_indices)):
        for es_4 in range(es_slen):
            es_vsub = numpy.amin(es_trans[es_3][(es_4):(es_4+7)]) * escape_params[0][es_3] + numpy.amax(es_trans[es_3][(es_4):(es_4+7)]) * escape_params[1][es_3] + escape_params[2][es_3]
            es_vlist.append(es_vsub)
    es_vsel = []
    for es_5 in range(len(es_sel1)):
        es_vsel.extend(es_vlist[((es_sel1[es_5])*es_slen):((es_sel1[es_5]+1)*es_slen)])
    for es_6 in range(len(es_sel2[0])):
        for es_7 in range(es_slen):
            es_vsel.append(es_vlist[((es_sel2[0][es_6])*es_slen)+es_7] - es_vlist[((es_sel2[1][es_6])*es_slen)+es_7])
    return es_vsel

def phosforus_scorer(ps_sites, ps_values, ps_clfs):
	ps_ilen = len(ps_clfs[0])-1
	ps_tlen = len(ps_values) / ps_ilen
	ps_scorelist = []
	ps_complist = []

	for ps_1 in range(len(ps_sites)):
		ps_subscore = []
		for ps_2 in range(ps_ilen):
			ps_tarvalues = [ps_values[(ps_2*ps_tlen+(ps_sites[ps_1][0]-14)):(ps_2*ps_tlen+(ps_sites[ps_1][0]+7))]]
			ps_subscore.append(float(ps_clfs[ps_sites[ps_1][1]][ps_2+1].predict_log_proba(ps_tarvalues)[:, 1]))
		ps_xscore = ps_clfs[ps_sites[ps_1][1]][0].predict_proba([ps_subscore])
		ps_tscore = ps_clfs[ps_sites[ps_1][1]][0].predict_log_proba([ps_subscore])
		ps_scorelist.append((ps_tscore[:, 1] - ps_tscore[:, 0])[0])
		if ps_scorelist[-1] > 0:
			ps_complist.append([ps_sites[ps_1][0]-13, p_class[ps_sites[ps_1][1]], "PHOSPHO", ps_scorelist[-1], ((ps_xscore[:, 1])[0])])
		else:
			ps_complist.append([ps_sites[ps_1][0]-13, p_class[ps_sites[ps_1][1]], "NONPHOS", ps_scorelist[-1], ((ps_xscore[:, 1])[0])])

	return ps_complist

####
