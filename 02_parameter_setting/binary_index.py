#!/usr/bin/env python

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
        sst = line.rstrip('\n').split('\t')[3][0:29]
        p_list.append(sst)

with open("training_dataset/nonphos_screened_"+sys.argv[1]+".txt") as n_file:
    for line in n_file:
        sst = line.rstrip('\n').split('\t')[3][0:29]
        n_list.append(sst)

#### Main body

score_compile = [[], []]

for i_1 in range(iterat_cycle):
    random.shuffle(p_list)
    random.shuffle(n_list)
    
    p_sample = random.sample(p_list, min(len(p_list), len(n_list)))
    n_sample = random.sample(n_list, min(len(p_list), len(n_list)))
    tlen = int(numpy.floor(float(len(p_sample))*0.9))
    tfloat = float(tlen)

    t_blist = []
    t_stats = []

    for i_2 in range(len(p_sample)):
        p_sub = []
        n_sub = []
        for i_3 in range(29):
            det_p = 0
            for i_4 in range(21):
                if p_sample[i_2][i_3] == aa_list[i_4]:
                    p_sub.append(1)
                    det_p += 1
                else:
                    p_sub.append(0)
                    if i_4 == 20 and det_p < 1:
                        p_sub[-1] += 1
                        det_p += 1
            det_n = 0
            for i_5 in range(21):
                if n_sample[i_2][i_3] == aa_list[i_5]:
                    n_sub.append(1)
                    det_n += 1
                else:
                    n_sub.append(0)
                    if i_5 == 20 and det_n < 1:
                        n_sub[-1] += 1
                        det_n += 1
        t_blist.append(p_sub)
        t_blist.append(n_sub)
        t_stats.append(1)
        t_stats.append(0)
    
    score_sub_tr = []
    score_sub_ts = []

    t_barray = numpy.asarray(t_blist)

    tlem = tlen - 100
    clf_b = GradientBoostingClassifier()
    clf_b.fit(t_barray[:2*tlen, :], t_stats[:2*tlen])
    b_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_b.predict_log_proba(t_barray[:2*tlen, :])[:, 1])
    b_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_b.predict_log_proba(t_barray[2*tlen:, :])[:, 1])

    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "BINARY_FEATURE", b_trscore, b_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(b_trscore)
    score_sub_ts.append(b_tsscore)
    
    p_clist = []
    n_clist = []
    for i_12 in range(len(p_sample)):
        p_clist.append(seq_codifier(p_sample[i_12]))
        n_clist.append(seq_codifier(n_sample[i_12]))

    t_vlist = []
    t_elist = []
    t_tlist = []

    for i_6 in range(len(p_sample)):
        t_vlist.append(index_scorer(p_clist[i_6], ind_list, weight_types))
        t_vlist.append(index_scorer(n_clist[i_6], ind_list, weight_types))
        t_elist.append(escape_scorer(p_clist[i_6], escape_list, es_presel_1, es_presel_2))
        t_elist.append(escape_scorer(n_clist[i_6], escape_list, es_presel_1, es_presel_2))
        #t_tlist.append(index_scorer(p_clist[i_6], ind_list, weight_types))
        #t_tlist[-1].extend(escape_scorer(p_clist[i_6], escape_list, es_presel_1, es_presel_2))
        #t_tlist.append(index_scorer(n_clist[i_6], ind_list, weight_types))
        #t_tlist[-1].extend(escape_scorer(n_clist[i_6], escape_list, es_presel_1, es_presel_2))

    t_varray = numpy.asarray(t_vlist)
    t_earray = numpy.asarray(t_elist)
    #t_tarray = numpy.asarray(t_tlist)

    for i_7 in range(10):
        clf_v = GaussianNB()
        clf_v.fit(t_varray[:2*tlen, (i_7*21):((i_7+1)*21)], t_stats[:2*tlen])
        v_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_v.predict_log_proba(t_varray[:2*tlen, (i_7*21):((i_7+1)*21)])[:, 1])
        v_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_v.predict_log_proba(t_varray[2*tlen:, (i_7*21):((i_7+1)*21)])[:, 1])
        print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], ind_label[i_7], v_trscore, v_tsscore
        t_tlist.append(numpy.ndarray.tolist(clf_v.predict_log_proba(t_varray[:, (i_7*21):((i_7+1)*21)])[:, 1]))
        score_sub_tr.append(v_trscore)
        score_sub_ts.append(v_tsscore)

    for i_8 in range(8):
        clf_e = GaussianNB()
        clf_e.fit(t_earray[:2*tlen, (i_8*21):((i_8+1)*21)], t_stats[:2*tlen])
        e_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_e.predict_log_proba(t_earray[:2*tlen, (i_8*21):((i_8+1)*21)])[:, 1])
        e_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_e.predict_log_proba(t_earray[2*tlen:, (i_8*21):((i_8+1)*21)])[:, 1])
        print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], es_prelabel[i_8], e_trscore, e_tsscore
        t_tlist.append(numpy.ndarray.tolist(clf_e.predict_log_proba(t_earray[:, (i_8*21):((i_8+1)*21)])[:, 1]))
        score_sub_tr.append(e_trscore)
        score_sub_ts.append(e_tsscore)
    
    t_tarray = numpy.asarray(numpy.transpose(t_tlist))
    
    ####

    clf_a = GradientBoostingClassifier()
    clf_a.fit(t_varray[:2*tlen, 0:126], t_stats[:2*tlen])
    a_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_a.predict_log_proba(t_varray[:2*tlen, 0:126])[:, 1])
    a_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_a.predict_log_proba(t_varray[2*tlen:, 0:126])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "HOR-GBoost", a_trscore, a_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(a_trscore)
    score_sub_ts.append(a_tsscore)

    clf_b = GradientBoostingClassifier()
    clf_b.fit(t_varray[:2*tlen, 126:210], t_stats[:2*tlen])
    b_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_b.predict_log_proba(t_varray[:2*tlen, 126:210])[:, 1])
    b_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_b.predict_log_proba(t_varray[2*tlen:, 126:210])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "VER-GBoost", b_trscore, b_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(b_trscore)
    score_sub_ts.append(b_tsscore)

    clf_c = GradientBoostingClassifier()
    clf_c.fit(t_earray[:2*tlen, :], t_stats[:2*tlen])
    c_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_c.predict_log_proba(t_earray[:2*tlen, :])[:, 1])
    c_tsscore = metrics.roc_auc_score(t_stats[2*tlen:], clf_c.predict_log_proba(t_earray[2*tlen:, :])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "ESC-GBoost", c_trscore, c_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(c_trscore)
    score_sub_ts.append(c_tsscore)
    
    ####

    clf_x = GradientBoostingClassifier()
    clf_x.fit(t_tarray[:2*tlen, 0:6], t_stats[:2*tlen])
    x_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_x.predict_log_proba(t_tarray[:2*tlen, 0:6])[:, 1])
    x_tsscore = metrics.roc_auc_score(t_stats[2*tlem:], clf_x.predict_log_proba(t_tarray[2*tlem:, 0:6])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "HOR-NESTED", x_trscore, x_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(x_trscore)
    score_sub_ts.append(x_tsscore)

    clf_y = GradientBoostingClassifier()
    clf_y.fit(t_tarray[:2*tlen, 6:10], t_stats[:2*tlen])
    y_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_y.predict_log_proba(t_tarray[:2*tlen, 6:10])[:, 1])
    y_tsscore = metrics.roc_auc_score(t_stats[2*tlem:], clf_y.predict_log_proba(t_tarray[2*tlem:, 6:10])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "VER-NESTED", y_trscore, y_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(y_trscore)
    score_sub_ts.append(y_tsscore)

    clf_z = GradientBoostingClassifier()
    clf_z.fit(t_tarray[:2*tlen, 10:18], t_stats[:2*tlen])
    z_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_z.predict_log_proba(t_tarray[:2*tlen, 10:18])[:, 1])
    z_tsscore = metrics.roc_auc_score(t_stats[2*tlem:], clf_z.predict_log_proba(t_tarray[2*tlem:, 10:18])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "ESC-NESTED", z_trscore, z_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(z_trscore)
    score_sub_ts.append(z_tsscore)
    
    ####

    clf_t = GradientBoostingClassifier(n_estimators = 50)
    clf_t.fit(t_tarray[:2*tlen, :], t_stats[:2*tlen])
    t_trscore = metrics.roc_auc_score(t_stats[:2*tlen], clf_t.predict_log_proba(t_tarray[:2*tlen, :])[:, 1])
    t_tsscore = metrics.roc_auc_score(t_stats[2*tlem:], clf_t.predict_log_proba(t_tarray[2*tlem:, :])[:, 1])
    print "#### Cycle " + str(i_1+1) + ":", p_class[int(sys.argv[1])], "COMPOSITE", t_trscore, t_tsscore, str(datetime.datetime.now())[8:19]
    score_sub_tr.append(t_trscore)
    score_sub_ts.append(t_tsscore)

    score_compile[0].append(score_sub_tr)
    score_compile[1].append(score_sub_ts)

output_label = ["BINARY_FEATURE", ind_label[0], ind_label[1], ind_label[2], ind_label[3], ind_label[4], ind_label[5], 
ind_label[6], ind_label[7], ind_label[8], ind_label[9], es_prelabel[0], es_prelabel[1], 
es_prelabel[2], es_prelabel[3], es_prelabel[4], es_prelabel[5], es_prelabel[6], es_prelabel[7], "HOR-GBoost", "VER-GBoost", "ESC-GBoost", "HOR-NESTED", "VER-NESTED", "ESC-NESTED", "COMPOSITE"]

score_train_array = numpy.asarray(score_compile[0])
score_test_array = numpy.asarray(score_compile[1])

out_file = open("index_analysis_result_"+sys.argv[1]+".txt", 'w')

for i_9 in range(len(output_label)):
    if i_9 == 0:
        print p_class[int(sys.argv[1])], output_label[i_9], "AUROC_Training", "AUROC_Testing", "SD_Training", "SD_Testing"
        out_file.write(p_class[int(sys.argv[1])] + '\t' + output_label[i_9] + '\t' + "AUROC_Training" + '\t' + "AUROC_Testing" + '\t' + "SD_Training" + '\t' + "SD_Testing" + '\n')
    print p_class[int(sys.argv[1])], output_label[i_9], numpy.mean(score_train_array[:, i_9]), numpy.mean(score_test_array[:, i_9]), numpy.std(score_train_array[:, i_9]), numpy.std(score_test_array[:, i_9])
    out_file.write(p_class[int(sys.argv[1])] + '\t' + output_label[i_9] + '\t' + str(numpy.mean(score_train_array[:, i_9])) + '\t' + str(numpy.mean(score_test_array[:, i_9])) + '\t' + str(numpy.std(score_train_array[:, i_9])) + '\t' + str(numpy.std(score_test_array[:, i_9])) + '\n')

out_file.close()
