#!/usr/bin/env python

#### Generic module import
import datetime
import numpy
import pickle
import random
import sys

from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score
from sklearn.naive_bayes import GaussianNB
from sklearn import metrics

#### Universal variables / values

aa_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y", "_"]
p_class = ["S-", "T-", "Y-", "SP", "TP"]
pswm_xrange = [9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29]
pswm_preset = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
escape_params = [[0.8195, 0.8195, 0.7665, 0.7665, 0.7791, 0.7791, 0.7047, 0.7047, 0.7023, 0.7023, 0.7801, 0.7801, 0.7115, 0.7115, 0.7319, 0.7319], 
[0.7492, 0.7492, 0.7632, 0.7632, 0.7524, 0.7524, 0.7507, 0.7507, 0.7479, 0.7479, 0.6708, 0.6708, 0.8047, 0.8047, 0.6782, 0.6782], 
[4696.0, 4696.0, -5068.0, -5068.0, 6195.0, 6195.0, 1998.0, 1998.0, -4154.0, -4154.0, 447.53, 447.53, -94.37, -94.37, 3872.0, 3872.0]]

weight_option = [[1.0], [1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0]]
weight_types = [1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
es_presel_1 = [2, 6, 9, 12]
es_presel_2 = [[0, 2, 4, 2], [8, 10, 12, 4]]
es_prelabel = ['dHapN', 'TdSN', 'ddGD', 'dHpD', 'dGN-dGD', 'dHapN-dHapD', 'dHpN-dHpD', 'dHapN-dHpN']

#### (Temporary) functions

def seq_codifier(sc_seq):
    sc_list = []
    for sc_1 in range(len(sc_seq)):
        sc_det_1 = 0
        for sc_2 in range(21):
            if sc_seq[sc_1] == aa_list[sc_2]:
                sc_list.append(sc_2)
                sc_det_1 += 1
        if sc_det_1 < 1:
            sc_list.append(20)
            sc_det_1 += 1
    return sc_list

def td_matrix(tm_dim1, tm_dim2, tm_base):
    td_ret = []
    for tm_1 in range(tm_dim1):
        td_sub = []
        for tm_2 in range(tm_dim2):
            td_sub.append((0.0 + tm_base))
        td_ret.append(td_sub)
    return td_ret

def pswm_scorer(ps_code, ps_pswm, ps_range):
    ps_ret = 0.0
    ps_cut = (29-ps_range)/2
    for ps_1 in range(ps_range):
        ps_ret += ps_pswm[ps_cut+ps_1][ps_code[ps_cut+ps_1]]
    return ps_ret

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

#### Index file import

ind_label = []
ind_list = []
escape_label = ['dGN', 'ddGN', 'dHapN', 'ddHapN', 'dHpN', 'ddHpN', 'TdSN', 'dTdSN', 'dGD', 'ddGD', 'dHapD', 'ddHapD', 'dHpD', 'ddHpD', 'TdSD', 'dTdSD']
escape_raw = []

with open("preset_indices/index_reselect.txt") as ind_file:
    for line in ind_file:
        ssl = line.rstrip('\n').split()
        ind_label.append(ssl[0])
        ind_list.append([float(i) for i in ssl[1:]])

with open("preset_indices/eScape_sorted_filled.csv") as ind_file:
    for line in ind_file:
        ssl = line.rstrip('\n').split()
        escape_raw.append([float(i) for i in ssl[1:17]])

escape_list = numpy.transpose(escape_raw)

#### Sequence file import

p_list = []
n_list = []

with open("training_dataset/phos_screened_"+sys.argv[1]+".txt") as p_file:
    for line in p_file:
        sst = line.rstrip('\n').split('\t')[3]
        p_list.append(sst)

with open("training_dataset/nonphos_screened_"+sys.argv[1]+".txt") as n_file:
    for line in n_file:
        sst = line.rstrip('\n').split('\t')[3]
        n_list.append(sst)

#### P-sample equalizer

#pn_ratio = float(len(n_list))/float(len(p_list))
#pn_residual = len(n_list) - len(p_list)*int(numpy.floor(pn_ratio))

#p_xsample = []

#for i_0_1 in range(int(numpy.floor(pn_ratio))):
    #p_xsample.extend(p_list)

#p_xsample.extend(random.sample(p_list, pn_residual))

#print len(p_list), len(n_list), len(p_xsample)

print len(p_list), len(n_list)
weight_base = float(len(p_list) + len(n_list))

#### Pickling body

t_clist = []
t_stats = []
t_weight = []

for i_1 in range(len(p_list)):
    t_clist.append(seq_codifier(p_list[i_1]))
    t_stats.append(1)
    t_weight.append(float(len(n_list)) / weight_base)

for i_2 in range(len(n_list)):
    t_clist.append(seq_codifier(n_list[i_2]))
    t_stats.append(0)
    t_weight.append(float(len(p_list)) / weight_base)

print "#### Sequence codification completed //", str(datetime.datetime.now())[8:19]

t_vlist = []

for i_3 in range(len(t_clist)):
    t_ival = index_scorer(t_clist[i_3], ind_list, weight_types)
    t_eval = escape_scorer(t_clist[i_3], escape_list, es_presel_1, es_presel_2)
    t_vlist.append(t_ival + t_eval)
    if i_3 % 10000 == 9999:
        print "#### Sequence numericization under process: " + str(i_3+1) +"th sequence //", str(datetime.datetime.now())[8:19]

print "#### Sequence numericization completed //", str(datetime.datetime.now())[8:19]

t_varray = numpy.asarray(t_vlist)
t_tlist = []

for i_4 in range(18):
    if i_4 < 10:
        clf_part = GaussianNB()
        clf_part.fit(t_varray[:, (i_4*21):((i_4+1)*21)], t_stats, sample_weight = t_weight)
        t_tlist.append(numpy.ndarray.tolist(clf_part.predict_log_proba(t_varray[:, (i_4*21):((i_4+1)*21)])[:, 1]))
        clf_filename = "preset_params/class_"+sys.argv[1]+"/param_"+sys.argv[1]+"_0"+str(i_4)+".txt"
        pickle.dump(clf_part, open(clf_filename, 'w'))
        print "#### Sub-classifier #0" + str(i_4)+" fitting completed //", str(datetime.datetime.now())[8:19]

    elif i_4 >= 10:
        clf_part = GaussianNB()
        clf_part.fit(t_varray[:, (i_4*21):((i_4+1)*21)], t_stats, sample_weight = t_weight)
        t_tlist.append(numpy.ndarray.tolist(clf_part.predict_log_proba(t_varray[:, (i_4*21):((i_4+1)*21)])[:, 1]))
        clf_filename = "preset_params/class_"+sys.argv[1]+"/param_"+sys.argv[1]+"_"+str(i_4)+".txt"
        pickle.dump(clf_part, open(clf_filename, 'w'))
        print "#### Sub-classifier #" + str(i_4)+" fitting completed //", str(datetime.datetime.now())[8:19]

t_tarray = numpy.asarray(numpy.transpose(t_tlist))

clf_total = GradientBoostingClassifier()
clf_total.fit(t_tarray, t_stats, sample_weight = t_weight)

clf_filename = "preset_params/class_"+sys.argv[1]+"/metaparam_"+sys.argv[1]+".txt"
pickle.dump(clf_total, open(clf_filename, 'w'))
print "#### Metapredictor fitting completed //", str(datetime.datetime.now())[8:19]

#### Pickling test

t_plist = []

eighteen_label = [ind_label[0], ind_label[1], ind_label[2], ind_label[3], ind_label[4], ind_label[5], 
ind_label[6], ind_label[7], ind_label[8], ind_label[9], es_prelabel[0], es_prelabel[1], 
es_prelabel[2], es_prelabel[3], es_prelabel[4], es_prelabel[5], es_prelabel[6], es_prelabel[7]]

for i_5 in range(18):
    if i_5 < 10:
        clf_filename = "preset_params/class_"+sys.argv[1]+"/param_"+sys.argv[1]+"_0"+str(i_5)+".txt"
        clf_pp = pickle.load(open(clf_filename))
        pp_pred = clf_pp.predict(t_varray[:, (i_5*21):((i_5+1)*21)])
        pp_score = clf_pp.score(t_varray[:, (i_5*21):((i_5+1)*21)], t_stats, sample_weight = t_weight)
        t_plist.append(numpy.ndarray.tolist(clf_pp.predict_log_proba(t_varray[:, (i_5*21):((i_5+1)*21)])[:, 1]))
        print p_class[int(sys.argv[1])], eighteen_label[i_5], pp_score, metrics.precision_score(t_stats, pp_pred), metrics.recall_score(t_stats, pp_pred), metrics.f1_score(t_stats, pp_pred), metrics.matthews_corrcoef(t_stats, pp_pred), metrics.roc_auc_score(t_stats, t_plist[-1])
    elif i_5 >= 10:
        clf_filename = "preset_params/class_"+sys.argv[1]+"/param_"+sys.argv[1]+"_"+str(i_5)+".txt"
        clf_pp = pickle.load(open(clf_filename))
        pp_pred = clf_pp.predict(t_varray[:, (i_5*21):((i_5+1)*21)])
        pp_score = clf_pp.score(t_varray[:, (i_5*21):((i_5+1)*21)], t_stats, sample_weight = t_weight)
        t_plist.append(numpy.ndarray.tolist(clf_pp.predict_log_proba(t_varray[:, (i_5*21):((i_5+1)*21)])[:, 1]))
        print p_class[int(sys.argv[1])], eighteen_label[i_5], pp_score, metrics.precision_score(t_stats, pp_pred), metrics.recall_score(t_stats, pp_pred), metrics.f1_score(t_stats, pp_pred), metrics.matthews_corrcoef(t_stats, pp_pred), metrics.roc_auc_score(t_stats, t_plist[-1])

t_parray = numpy.asarray(numpy.transpose(t_plist))

clf_filename = "preset_params/class_"+sys.argv[1]+"/metaparam_"+sys.argv[1]+".txt"
clf_pt = pickle.load(open(clf_filename))
pt_pred = clf_pt.predict(t_tarray)
pt_grad = clf_pt.predict_log_proba(t_tarray)[:, 1]
pt_score = clf_pt.score(t_tarray, t_stats, sample_weight = t_weight)

print p_class[int(sys.argv[1])], "WHOLE", pt_score, metrics.precision_score(t_stats, pt_pred), metrics.recall_score(t_stats, pt_pred), metrics.f1_score(t_stats, pt_pred), metrics.matthews_corrcoef(t_stats, pt_pred), metrics.roc_auc_score(t_stats, pt_grad)




