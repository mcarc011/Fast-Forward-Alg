#%%
import numpy as np
import sympy
from itertools import combinations
from itertools import product

def jandeterms(filestr):
    with open(filestr,'r') as f:
        text = f.read()
        f.close()
    monomials = text.split('\\')
    monomials = [line.split(':')[1] for line in monomials if ':' in line]
    monomials = [line.replace(' & ','|') for line in monomials]
    monomials = [line.replace('&','') for line in monomials]
    return monomials

def chiraldic(jelist):
    chirals = {}
    for mon in jelist:
        Jterm,Eterm = mon.split('|')[1],mon.split('|')[2]
        Jterm,Eterm,mon = Jterm.replace('|',''), Eterm.replace('|',''), mon.replace('|','')
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


def findtmatrix(kmatrix, nodenum):
    tlist = []
    for t in product(range(2),repeat=nodenum+3):
        store = True
        for k in kmatrix:
            if np.dot(t,k)<0:
                store = False
                break
        if store and t not in tlist:
            tlist += [list(t)]
    tlist = [np.array(i) for i in tlist]
    tvectors = tlist.copy()

    def treducefunc(tmatrix,irepeat):
        tnew = []
        for n,ti in enumerate(tmatrix):
            save = True
            for addedarr in combinations(range(len(tmatrix)),irepeat):
                if n not in addedarr:
                    tn = addedarr[0]
                    for a in addedarr[1:]:
                        tn += tmatrix[a]
                    if all(tn==ti):
                        save = False
                if not save:
                    break
            if save:
                tnew +=[ti]
        return tnew

    i = 2
    while i<len(tlist):
        tlist = treducefunc(tlist,i)
        i+=1
        
    return np.matrix([t0 for t0 in tlist if np.sum(t0)>0]),tvectors
    #return [tmatrix[i] for i in inds]



jeterms = jandeterms('model3.txt')
xdict = chiraldic(jeterms)
k = findkmatrix(xdict,8)
#K = np.matrix([k[n] for n in k])

# #SPPxC
# K = [
#     [1,1,0,0,0,0,0,0,0,0],
#     [0,1,1,0,0,0,0,0,0,0],
#     [0,0,0,1,0,0,-1,0,0,0],
#     [0,0,0,0,1,0,1,0,0,0],
#     [0,0,0,0,0,1,1,0,0,0],
#     [0,0,0,0,0,0,0,1,1,1]
#     ]
# K = np.transpose(K)

T,tvec = findtmatrix(K,3)
print(K*np.transpose(T))

#%%