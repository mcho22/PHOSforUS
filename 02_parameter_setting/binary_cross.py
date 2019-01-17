#!/usr/bin/env python

#### Cross-usage of binary feature-based (PSSM-based) predictor 
#### Calculates differences between subclasses divided by phosphorylation sites / +1 proline

#### Sample result ()

# S- 0.879895103973 0.815059171367 0.762053070662 0.649820259651 0.674778273368 0.673587581434
# T- 0.918203127221 0.791296970748 0.75605359001  0.669760955024 0.697197941959 0.698064516129
# Y- 0.93332128     0.661973500336 0.666030226179 0.733300553802 0.655922571295 0.645016821076
# SP 0.870917665184 0.751379011438 0.72836445751  0.683022719707 0.750161179892 0.756850522296
# TP 0.917417736906 0.718072565772 0.697125221122 0.658577665701 0.719547969438 0.755019007114

#### Generic module import 
import datetime
import numpy
import random
import sys

from sklearn import metrics
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.model_selection import cross_val_score
from sklearn.naive_bayes import GaussianNB

#### Universal variables / values

iterat_cycle = 100

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

def binary_codifier(bc_seq):
    bc_qlist = []
    for bc_1 in range(len(bc_seq)):
        det_bc = 0
        for bc_2 in range(len(aa_list)):
            if bc_seq[bc_1] == aa_list[bc_2]:
                bc_qlist.append(1)
                det_bc += 1
            else:
                if bc_2 == (len(aa_list)-1) and det_bc < 1:
                    bc_qlist.append(1)
                    det_bc += 1
                else:
                    bc_qlist.append(0)
    return bc_qlist

#

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

p_list = [[], [], [], [], []]
n_list = [[], [], [], [], []]

for i_0 in range(5):
    with open("training_dataset/phos_screened_"+str(i_0)+".txt") as p_file:
        for line in p_file:
            sst = line.rstrip('\n').split('\t')[3][0:29]
            p_list[i_0].append(sst)

    with open("training_dataset/nonphos_screened_"+str(i_0)+".txt") as n_file:
        for line in n_file:
            sst = line.rstrip('\n').split('\t')[3][0:29]
            n_list[i_0].append(sst)

#### Main body

score_compile = [[], []]

for i_1 in range(iterat_cycle):
    t_str = []
    t_sts = []
    for i_2 in range(5):
        p_sample = random.sample(p_list[i_2], min(len(p_list[i_2]), len(n_list[i_2])))
        n_sample = random.sample(n_list[i_2], min(len(p_list[i_2]), len(n_list[i_2])))
        tlen = int(numpy.floor(float(len(p_sample))*0.9))
        tshort = len(p_sample)-tlen
        tfloat = float(tlen)
        
        t_blist = []
        t_stats = []
        
        for i_3 in range(len(p_sample)):
            t_blist.append(binary_codifier(p_sample[i_3]))
            t_blist.append(binary_codifier(n_sample[i_3]))
            t_stats.append(1)
            t_stats.append(0)
        
        t_barray = numpy.asarray(t_blist)
        clf_base = GradientBoostingClassifier()
        clf_base.fit(t_barray[:2*tlen, :], t_stats[:2*tlen])
        #b_trscore = clf_base.score(t_barray[:2*tlen, :], t_stats[:2*tlen])
        #b_tsscore = clf_base.score(t_barray[2*tlen:, :], t_stats[2*tlen:])
        b_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_base.predict_log_proba(t_barray[:2*tlen, :])[:, 1])
        b_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_base.predict_log_proba(t_barray[2*tlen:, :])[:, 1])
        
        t_str.append(b_trscore)
        
        for i_4 in range(5):
            if i_2 == i_4:
                t_sts.append(b_tsscore)
            else:
                p_comp = random.sample(p_list[i_4], tshort)
                n_comp = random.sample(n_list[i_4], tshort)
                
                c_blist = []
                c_stats = []
                
                for i_5 in range(len(p_comp)):
                    c_blist.append(binary_codifier(p_comp[i_5]))
                    c_blist.append(binary_codifier(n_comp[i_5]))
                    c_stats.append(1)
                    c_stats.append(0)
                
                c_barray = numpy.asarray(c_blist)
                #c_score = clf_base.score(c_barray, c_stats)
                c_score = metrics.roc_auc_score(c_stats, clf_base.predict_log_proba(c_barray)[:, 1])
                t_sts.append(c_score)
    
    print i_1+1, "S-", str(t_str[0])[0:5], str(t_sts[0])[0:5], str(t_sts[1])[0:5], str(t_sts[2])[0:5], str(t_sts[3])[0:5], str(t_sts[4])[0:5]
    print i_1+1, "T-", str(t_str[1])[0:5], str(t_sts[5])[0:5], str(t_sts[6])[0:5], str(t_sts[7])[0:5], str(t_sts[8])[0:5], str(t_sts[9])[0:5]
    print i_1+1, "Y-", str(t_str[2])[0:5], str(t_sts[10])[0:5], str(t_sts[11])[0:5], str(t_sts[12])[0:5], str(t_sts[13])[0:5], str(t_sts[14])[0:5]
    print i_1+1, "SP", str(t_str[3])[0:5], str(t_sts[15])[0:5], str(t_sts[16])[0:5], str(t_sts[17])[0:5], str(t_sts[18])[0:5], str(t_sts[19])[0:5]
    print i_1+1, "TP", str(t_str[4])[0:5], str(t_sts[20])[0:5], str(t_sts[21])[0:5], str(t_sts[22])[0:5], str(t_sts[23])[0:5], str(t_sts[24])[0:5]
    
    score_compile[0].append(t_str)
    score_compile[1].append(t_sts)

score_trarray = numpy.asarray(score_compile[0])
score_tsarray = numpy.asarray(score_compile[1])

out_file = open("cross_analysis_result.txt", 'w')

for i_6 in range(5):
    print p_class[i_6], numpy.mean(score_trarray[:, i_6]), numpy.mean(score_tsarray[:, (i_6*5)+0]), numpy.mean(score_tsarray[:, (i_6*5)+1]), numpy.mean(score_tsarray[:, (i_6*5)+2]), numpy.mean(score_tsarray[:, (i_6*5)+3]), numpy.mean(score_tsarray[:, (i_6*5)+4])
    out_file.write(p_class[i_6] + '\t' + str(numpy.mean(score_trarray[:, i_6])) + '\t' + str(numpy.mean(score_tsarray[:, (i_6*5)+0])) + '\t' + str(numpy.mean(score_tsarray[:, (i_6*5)+1])) + '\t' + str(numpy.mean(score_tsarray[:, (i_6*5)+2])) + '\t' + str(numpy.mean(score_tsarray[:, (i_6*5)+3])) + '\t' + str(numpy.mean(score_tsarray[:, (i_6*5)+4])) + '\n')

out_file.close()
