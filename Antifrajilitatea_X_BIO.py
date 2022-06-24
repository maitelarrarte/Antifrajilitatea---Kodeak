# -*- coding: utf-8 -*-
"""
@author: okarim
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
    plt.title("Batez besteko antizaurgarritasuna Sare Boolear biologikotarako O=1")
    plt.ylabel("Antizaurgarritasuna")
    plt.xlabel("X")
    
    colors=['b', 'orange', 'g', 'brown', 'purple', 'pink', 'y' ]
    
    red=rbn.RBN()
    r=100 #100
    
    red.CreateBioNet(1)
    f=[]
  
    for i in trange(r):
 
        f.append(red.antifragile(100, runs=10, O=1, fraction=1)) #100,runs=10
        
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="CD4+ T Zelula", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[0] )
    
    red.CreateBioNet(2)
    f=[]
    for i in trange(r):
        f.append(red.antifragile(100, runs=10, O=1, fraction=1)) #100,runs=10
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Ugaztunak", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[1] )
    
    red.CreateBioNet(3)
    f=[]
    for i in trange(r):
        f.append(red.antifragile(100, runs=10, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Bihotza", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[2] )

    red.CreateBioNet(4)
    f=[]
    for i in trange(r):
        f.append(red.antifragile(100, runs=1, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Mikrobioma", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[3] )

    red.CreateBioNet(5)
    f=[]
    for i in trange(r):
        f.append(red.antifragile(100, runs=10, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Heriotza", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[4] )
    
    red.CreateBioNet(6)
    f=[]
    for i in trange(r):
        f.append(red.antifragile(100, runs=10, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="A. thaliana", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[5] )
    
    red.CreateBioNet(7)
    f=[]
    for i in trange(r):
        f.append(red.antifragile(100, runs=10, O=1, fraction=1))
    plt.errorbar(np.arange(1,red.N+1), np.mean(f,0), label="Tumorea", 
                 yerr=stats.sem(f, 0), ecolor='r', color=colors[6] )
    
    plt.legend()
    plt.savefig("Antifrajilitatea_X_BIO.pdf")
    print("--- %s seconds ---" % (time.time() - start_time))