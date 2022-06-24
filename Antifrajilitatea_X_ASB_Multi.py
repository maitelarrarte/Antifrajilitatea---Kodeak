# -*- coding: utf-8 -*-
"""

Antifrajilitatearen batez bestekoa X perturbatutako nodo kopuruaren funtzioan
geruza anitzeko ASBetarako
L zelulen arteko link kopurua ausaz finkatuko da aldiune bakoitzean

@author: Maite

L KONTUN HARTU GABE X ARDATZEN X JARRITA, K EZBERDINETAKO ANTIFRAGILITATEA

"""

import numpy as np
import rbn
from  tqdm import trange
import time
import numpy as np
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
 

if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    
    p=0.5
    T=100
    N=10
    Nc=5
    Ng=4
    L=50
    
    colors=['b', 'orange', 'g', 'brown', 'purple', 'pink', 'y' ]
    
    plt.title("Geruza anitzeko EAS")
    plt.ylabel("Antizaurgarritasuna")
    plt.xlabel("X")
        
    red=rbn.RBN()
    p1=np.zeros((5,N*Nc))
    r=10
    
    for K in trange(1, 6):        
        f=np.zeros(( r, N*Nc )) 
        i=0
        for x in trange(r):
            red.CreateNetMultilayer(K, N, p, Ng, Nc)
            f[i]=red.antifragileMultilayer(T, runs=10, O=1) #runs=10  
            i+=1
        g1=np.mean(f,0)
        p1[K-1]=np.sum(np.array(f)<0, axis=0)/r #probabilitatea
        plt.errorbar(np.arange(1,(N*Nc)+1), g1, label="K="+str(K), yerr=stats.sem(f, 0), ecolor='r', color=colors[K-1] )
    
    plt.legend()
    plt.savefig("Antifrajilitatea_RBN_Multi_X_N="+str(N)+"_Nc="+str(Nc)+"Ng="+str(Ng)+"_r="+str(r)+"_O=1.pdf")
       
    plt.show()
    
    
    
    print("--- %s seconds ---" % (time.time() - start_time))