import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
from os import mkdir

import numpy as np
from scipy.stats import truncnorm
import matplotlib.pyplot as plt
from os import mkdir

try:
    mkdir('./Combinatorics')
except:
    print('Combinatorics arleady exists')

def make_trunc_norm(partons, jets, std):
    array = np.empty([partons,jets,jets])
    for i in range(partons):
        array[i] =  truncnorm.rvs(a=(0 - .5) / std, b=(1 - 0.5) / std, loc=0.5, scale=std, size=[jets,jets])
        for j in range(jets):
            for k in range(jets):
                if k >= j:
                    array[i,j,k] = -999
    return array

def make_bimodal(partons, jets, std):
    array = np.empty([partons,jets,jets])
    for i in range(partons):
        norm = truncnorm.rvs(a=(0) / std, b=(1) / std,
                                 loc=0, scale=std, size=[jets, jets])
        modifier = (np.random.rand(jets, jets) > 0.5).astype(int)
        array[i] = np.where(modifier - norm > 0, modifier - norm, norm - modifier)
        for j in range(jets):
            for k in range(jets):
                if k >= j:
                    array[i,j,k] = -999
    return array

def check_all_jets(p, jets, checked=[]):
    p_max = -999
    for i in range(jets):
        for j in range(jets):
            if p[i,j] == -999:
                continue
            elif p[i,j] > p_max:
                if i not in checked and j not in checked:
                    p_max = p[i,j]
                    i_max = i
                    j_max = j

    return p_max, i_max, j_max

def SpaNet_Combinatorics(array, partons, jets):
    claimed = np.ones([partons,3]) * -999
    while -999 in claimed:
        for p in range(partons):
            if claimed[p][0] == -999:
                val_temp, i_temp, j_temp = check_all_jets(array[p], jets)
                checked = np.full_like(claimed, -999)
                while i_temp in claimed or j_temp in claimed:
                    val_temp, i_temp, j_temp = check_all_jets(array[p], jets, checked)
                    if i_temp in claimed:
                        where_i = np.nonzero(claimed == i_temp)[0][0]
                        if j_temp in claimed:
                            where_j = np.nonzero(claimed == j_temp)[0][0]
                            if val_temp > claimed[where_i][2]:
                                if val_temp > claimed[where_j][2]:
                                    claimed[where_i] = [-999, -999, -999]
                                    claimed[where_j] = [-999, -999, -999]
                                else:
                                    checked[where_j] = claimed.copy()[where_j]
                                    claimed[where_i] = claimed.copy()[where_i]
                            else:
                                checked[where_j] = claimed.copy()[where_j]
                                claimed[where_i] = claimed.copy()[where_i]
                        else:
                            if val_temp > claimed[where_i][2]:
                                claimed[where_i] = [-999, -999, -999]
                            else:
                                checked[where_i] = claimed.copy()[where_i]
                    elif j_temp in claimed:
                        where_j = np.nonzero(claimed == j_temp)[0][0]
                        if val_temp > claimed[where_j][2]:
                            claimed[where_j] = [-999, -999, -999]
                        else:
                            checked[where_j] = claimed.copy()[where_j]
                val_temp, i_temp, j_temp = check_all_jets(array[p], jets, checked)
                claimed[p] = [i_temp, j_temp, val_temp]
            else:
                continue

    return np.prod(claimed[:,2]), np.sum(claimed[:,2])

def PartonIteration_Combinatorics(array, partons, jets):
    sum_array = np.empty([partons])
    prod_array = np.empty([partons])
    for parton in range(partons):
        claimed = np.ones([partons, 3]) * -999
        a_max = np.max(array[parton])
        idx = np.nonzero(array[parton] == a_max)
        array[parton][idx] = 999
        while -999 in claimed:
            for p in range(partons):
                if claimed[p][0] == -999:
                    val_temp, i_temp, j_temp = check_all_jets(array[p], jets)
                    checked = np.full_like(claimed, -999)
                    while i_temp in claimed or j_temp in claimed:
                        val_temp, i_temp, j_temp = check_all_jets(array[p], jets, checked)
                        if i_temp in claimed:
                            where_i = np.nonzero(claimed == i_temp)[0][0]
                            if j_temp in claimed:
                                where_j = np.nonzero(claimed == j_temp)[0][0]
                                if val_temp > claimed[where_i][2]:
                                    if val_temp > claimed[where_j][2]:
                                        claimed[where_i] = [-999, -999, -999]
                                        claimed[where_j] = [-999, -999, -999]
                                    else:
                                        checked[where_j] = claimed.copy()[where_j]
                                        claimed[where_i] = claimed.copy()[where_i]
                                else:
                                    checked[where_j] = claimed.copy()[where_j]
                                    claimed[where_i] = claimed.copy()[where_i]
                            else:
                                if val_temp > claimed[where_i][2]:
                                    claimed[where_i] = [-999, -999, -999]
                                else:
                                    checked[where_i] = claimed.copy()[where_i]
                        elif j_temp in claimed:
                            where_j = np.nonzero(claimed == j_temp)[0][0]
                            if val_temp > claimed[where_j][2]:
                                claimed[where_j] = [-999, -999, -999]
                            else:
                                checked[where_j] = claimed.copy()[where_j]
                    val_temp, i_temp, j_temp = check_all_jets(array[p], jets, checked)
                    claimed[p] = [i_temp, j_temp, val_temp]
                else:
                    continue
        array[parton][idx] = a_max
        claimed[np.nonzero(claimed == 999)[0][0]][2] = a_max
        sum_array[parton] = np.sum(claimed[:,2])
        prod_array[parton] = np.prod(claimed[:,2])
        
    return np.max(prod_array), np.max(sum_array)

def run_models_and_plot(std, parts=[1, 2], js=[2, 4], dist='trunc_norm'):
    sum_perc = []
    prod_perc = []
    for partons, jets in zip(parts, js):
        print('partons: ' + str(partons))
        print('jets: ' + str(jets))
        sn_prods = []
        sn_sums = []
        pi_prods = []
        pi_sums = []
        for i in range(1000):
            if dist=='trunc_norm':
                array = make_trunc_norm(partons, jets, std)
            elif dist=='bimodal':
                array = make_bimodal(partons, jets, std)
            sn_prod, sn_sum = SpaNet_Combinatorics(array, partons, jets)
            pi_prod, pi_sum = PartonIteration_Combinatorics(array, partons, jets)
            sn_prods.append(sn_prod)
            sn_sums.append(sn_sum)
            pi_prods.append(pi_prod)
            pi_sums.append(pi_sum)
        sn_prods = np.array(sn_prods, dtype='float')
        sn_sums = np.array(sn_sums, dtype='float')
        pi_prods = np.array(pi_prods, dtype='float')
        pi_sums = np.array(pi_sums, dtype='float')
        try:
            sum_perc.append((len(pi_sums[pi_sums > sn_sums]) / len(pi_sums)) * 100)
        except:
            sum_perc.append(0)
        try:
            prod_perc.append((len(pi_prods[pi_prods > sn_prods]) / len(pi_prods)) * 100)
        except:
            prod_perc.append(0)

    pstring = []
    pnum =[]
    jstring = []
    jnum = []
    columns = []
    rows = []
    for i in range(len(parts)):
        pstring.append('p=')
        pnum.append(parts[i])
        jstring.append(' j=')
        jnum.append(js[i])
    for elem in zip(pstring, pnum, jstring, jnum):
        columns.append(elem[0] + str(elem[1]) + elem[2] + str(elem[3]))
    rows = (r'$\Pi P.I. > \Pi S.N.$', r'$\sum P.I. > \sum S.N.$')
    cell_text = ['%.1f%%' % x for x in sum_perc]
    cell_text2 = ['%.1f%%' % y for y in prod_perc]
    cell_text_array = [cell_text, cell_text2]
    plt.close()
    fig, ax = plt.subplots()

    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')

    ax.table(cellText=cell_text_array,
             rowLabels=rows,
             colLabels=columns,
             loc='center')
    plt.title('% Parton-Iteration Sln. >  SpaNet Sln. (n=1000)')

    fig.tight_layout()
    if dist=='trunc_norm':
        fig.savefig('./Combinatorics/Trunc_Norm_Table_std=' + str(std) + '_p=' + str(parts) + '_j=' + str(js) + '.jpeg')
    elif dist=='bimodal':
        fig.savefig('./Combinatorics/Bimodal_Table_std=' + str(std) + '_p=' + str(parts) + '_j=' + str(js) + '.jpeg')

    plt.close()
    plt.hist(array[array != -999].flatten(), bins=30)
    plt.title('loc=0.5, std=%.2f' % std)
    
    if dist=='trunc_norm':
        plt.savefig('./Combinatorics/Trunc_Norm_Distr_std=' + str(std) + '_p=' + str(parts) + '_j=' + str(js) + '.jpeg')
        plt.close()
    elif dist=='bimodal':
        plt.savefig('./Combinatorics/Bimodal_Distr_std=' + str(std) + '_p=' + str(parts) + '_j=' + str(js) + '.jpeg')
        plt.close()

# stds = [0.15, 0.2, 0.25, .3, .35]
# for std in stds:
#     run_models_and_plot(std)
