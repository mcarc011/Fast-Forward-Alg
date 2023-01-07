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

def findkmatrix(xdict,nodenum):
    chiralkeys = list(xdict.keys())
    c = 0
    for relabels in combinations(range(len(xdict)),nodenum+3):
        K = {}
        c +=1
        relabels = [chiralkeys[i] for i in list(relabels)]
        for chiral in chiralkeys:
            if chiral in relabels:
                temp = np.zeros(len(relabels))
                temp[relabels.index(chiral)] = 1
                K[chiral] = temp
            else:
                found = False
                for expression in xdict[chiral]:
                    for r in relabels:
                        expression = expression.replace(r,str(relabels.index(r)))
                    expression = expression.replace('*',',')
                    expression = expression.replace('inv','-')
                    if '_' not in expression:
                        found = True
                        Kvector = np.zeros(len(relabels))
                        for ind in expression.split(','):
                            if ind!='':
                                zerovector = np.zeros(len(relabels))
                                if '-' in ind:
                                    zerovector[abs(int(ind))] = -1
                                else:
                                    zerovector[abs(int(ind))] = 1
                                Kvector += zerovector
                        break
                if found:
                    K[chiral] = Kvector 
        if len(K) == len(xdict):
            break


jeterms = jandeterms('c4z4.txt')
xdict = chiraldic(jeterms)
nodenum = 4

#%%