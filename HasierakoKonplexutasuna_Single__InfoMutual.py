# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 19:39:26 2022

@author: Maite
"""

import rbn
import numpy as np
from math import log
from  tqdm import trange
import time
import matplotlib.pyplot as plt
from scipy import stats

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    N=5
    T=100
    p=0.5
    red=rbn.RBN()
    number_of_iterations=50 
    g1=[]
    yerr=[]
    
    rango=np.arange(1,6)
    for K in trange(1,6):
        I=[]
        for x in range(number_of_iterations):
            red.CreateNet(K, N, p)
            State = red.RunNet(2*T)
            mutual=rbn.mutual_info(State[-T:])/N
            I.append(4*mutual*(1-mutual))
        g1.append(np.mean(I))
        yerr.append(I)
    
    # make data:
    x = np.arange(1,K+1)
    y = g1
    # plot
    fig, ax = plt.subplots()
    plt.title("N="+str(N))
    plt.ylabel("C")
    plt.xlabel("K")
    ax.bar(x, y, width=0.8, align='center', yerr=stats.sem(yerr,1),
           ecolor='r', edgecolor="white", linewidth=0.7)
    ax.set(xlim=(0,6), xticks=np.arange(1, 6), ylim=(0,1))
    plt.show()
    plt.savefig("HasierakoKonplexitatea_N="+str(N)+".pdf")

    print("--- %s seconds ---" % (time.time() - start_time))
    
    State0=red.RunNet(2*T)
    State0=State0[-T:,:]
    
   
    
    print("I=",rbn.mutual_info(State0))
    
    State=red.RunNet(2*T,X=3,O=1)
    State=State[-T:,:]
    print("I perturbauta=",rbn.mutual_info(State))
    
    print("diferentzia=",rbn.mutual_info(State) - rbn.mutual_info(State0) )
    
    
    
    