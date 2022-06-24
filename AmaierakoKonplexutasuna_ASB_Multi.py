# -*- coding: utf-8 -*-
"""
@author: Maite Larrarte Mayoz

Amaierako konplexutasuna
geruza anitzeko ASBetarako
perturbatutako X nodoen funtzioan
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
    
    N=10
    Ng=4
    Nc=5
    p=0.5
    T=100
    
    colors=['b', 'orange', 'g', 'brown', 'purple']
    
    number_of_iterations=50 
    fraction=1
    plt.ylabel("C")
    plt.xlabel("X")
    red=rbn.RBN()

    for K in trange(1, 6):
        C=np.zeros((number_of_iterations,Nc*N))
        i=0
        for q in trange(number_of_iterations):
            red.CreateNetMultilayer(K, N, p, Ng, Nc)
            for XX in range(1,Nc*N):
                Cell_State=red.RunNetMultilayer(T,X=XX,O=1)
                for i in range(Nc):
                    State=np.array(Cell_State[i])
                    State=State[-T:]
                    Cell_State[i]=State
                C[q, XX-1]=rbn.complexityMultilayer(Cell_State)
            i+=1
        g1=np.mean(C, 0)

        plt.errorbar(np.arange(1,(N*Nc)+1), g1, label="K= "+str(K), 
                     yerr=stats.sem(C), ecolor='r', color=colors[K-1])
  
    plt.title("Amaierako konplexutasuna geruza anitzeko ASBentzako")
    plt.legend()
    
    plt.savefig("AmaierakoKonplexitatea_N="+str(N)+"_NC="+str(Nc)+"Ng="+str(Ng)+"multi.pdf")
    
    
    plt.show()
    
    
    print("--- %s seconds ---" % (time.time() - start_time))