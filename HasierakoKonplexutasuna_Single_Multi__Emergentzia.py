# -*- coding: utf-8 -*-
"""
@author: Maite Larrarte Mayoz

Emergentzian oinarritutako 
hasierako konplexutasuna 
geruza bakarreko eta anitzeko ASBetarako
"""

import rbn
import time
import matplotlib.pyplot as plt
import numpy as np
from  tqdm import trange
from scipy import stats

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    N=5
    p=0.5
    T=100
    
    number_of_iterations=50
    
    red=rbn.RBN()
    g1=[]
    yerr=[]
    
    rango=np.arange(1,6)
    for K in trange(1,6):
        C=[]
        for x in range(number_of_iterations):
            red.CreateNet(K, N, p)
            State = red.RunNet(2*T)
            C.append(np.mean(rbn.complexity(State[-T:])))
        g1.append(np.mean(C))
        yerr.append(C)
    
# make data:
x = np.arange(1,K+1)
y = g1

# plot
fig, ax = plt.subplots()
plt.title("N="+str(N)) 
plt.ylabel("C")
plt.xlabel("K")
plt.ylim(top=1)
ax.bar(x, y, width=0.8, align='center', yerr=stats.sem(yerr,1),
       ecolor='r', edgecolor="white", linewidth=0.7)
ax.set(xlim=(0,6), xticks=np.arange(1, 6), ylim=(0,1))
plt.savefig("HasierakoKonplexutasuna_N="+str(N)+"_SINGLE.pdf")
plt.show()

# MULTILAYER

N=10
Ng=4
Nc=5
p=0.5
T=100 
number_of_iterations=50 
    
red=rbn.RBN()
g1=[]
yerr=[]
    

for K in trange(1,6):
        C=[]
        for x in range(number_of_iterations):
            red.CreateNetMultilayer(K, N, p, Ng, Nc)
            Cell_State = red.RunNetMultilayer(2*T)
            for i in range(Nc):
                State=np.array(Cell_State[i])
                State=State[-T:]
                Cell_State[i]=State.tolist()
                
            
            C.append(np.mean(rbn.complexityMultilayer(Cell_State)))
            
        g1.append(np.mean(C))
        yerr.append(C)
    

# make data:
x = np.arange(1,K+1)
y = g1

# plot
fig, ax = plt.subplots()
plt.title("Hasierako konplexutasuna geruza anitzeko ASBentzako")
plt.ylabel("C0")
plt.xlabel("K")
plt.ylim(top=1)
ax.bar(x, y, width=0.8, align='center', yerr=stats.sem(yerr,1),
       ecolor='r', edgecolor="white", linewidth=0.7)
ax.set(xlim=(0,6), xticks=np.arange(1, 6), ylim=(0,1))
plt.savefig("HasierakoKonplexutasuna_N="+str(N)+"_NC="+str(Nc)+"Ng="+str(Ng)+"multi.pdf")
plt.show()

