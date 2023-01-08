#%%
import numpy as np
from itertools import combinations

def jandeterms(filestr):
    with open(filestr,'r') as f:
        text = f.read()
        f.close()
    monomials = text.split('\\')
    monomials = [line.split(':')[1] for line in monomials if ':' in line]
    monomials = [line.replace('&','') for line in monomials]
    return monomials

def chiraldic(jelist):
    chirals = {}
    for mon in jelist:
        Jterm,Eterm = mon.split('    ')
        Jterm,Eterm = Jterm.split('-'), Eterm.split('-')
        tempterm = {}
        for term in mon.split(' '):
            if '_' in term:
                if term not in chirals:
                    chirals[term] = []
                if term not in tempterm:
                    tempterm[term] = {}
                tempterm[term]['J'] = [str(Jterm.index(i)) for i in Jterm if term in i]
                tempterm[term]['E'] = [str(Eterm.index(i)) for i in Eterm if term in i]

        for term1 in tempterm:
            rel = ['','']
            for term2 in tempterm:
                if term1!=term2:
                    types = ['J','E']
                    for i,t in enumerate(types):
                        if tempterm[term1][t] != [] and tempterm[term2][t] != []:
                            if tempterm[term1][t] == tempterm[term2][t]:
                                rel[i] += 'inv'+term2+'*'
                            else:
                                rel[i] += term2+'*'
            for r in rel:
                if r !='':
                    chirals[term1] += [r]
    return chirals

def findkmatrix(xdict):
    chiralkeys = list(xdict.keys())

    'SEARCH FOR A RELABELING'
    for relabels in combinations(range(len(xdict)),nodenum+3):
        K = {}
        relabels = [chiralkeys[i] for i in list(relabels)]
        relabels = [i for i in xdict if 'P' in i or '{12}' in i]
        for chiral in chiralkeys:
            if chiral in relabels:
                K[chiral] = str(relabels.index(chiral))
            else:
                for expression in xdict[chiral]:
                    for r in relabels:
                        expression = expression.replace(r,str(relabels.index(r)))
                    expression = expression.replace('*',',')
                    expression = expression.replace('inv','-')
                    if chiral not in K:
                        K[chiral] = expression
                    if expression.count('_')<K[chiral].count('_'):
                        K[chiral] = expression
        Unlabeledp = [0]
        Unlabeled = [K[key].count('_') for key in K]
        while np.sum(Unlabeledp)<np.sum(Unlabeled):
            Unlabeled = [K[key].count('_') for key in K]
            for chiral in K:
                for item in K[chiral].split(','):
                    chiral2 = item.replace('-','')
                    if '_' in item and chiral2 in K:
                        if '_' not in K[chiral2]:
                            K[chiral] = K[chiral].replace(item+',',K[chiral2])
                            K[chiral] = K[chiral].replace(item,K[chiral2])
            Unlabeledp = [K[key].count('_') for key in K]
        if np.sum(Unlabeledp)==0:
            break

    'TURN THE STRING DATA INTO ARRAYS'
    for chiral in K:
        datum = K[chiral].split(',')
        kvec = np.zeros(nodenum+3)
        for d in datum:
            if d!='':
                zero = np.zeros(nodenum+3)
                if '-' in d:
                    zero[abs(int(d))] = -1
                else:
                    zero[abs(int(d))] = 1
                kvec += zero
        K[chiral] = kvec
    return K


jeterms = jandeterms('c4z4.txt')
xdict = chiraldic(jeterms)
k = findkmatrix(xdict)

#%%