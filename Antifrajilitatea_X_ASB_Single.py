# -*- coding: utf-8 -*-
"""
@author: Maite Larrarte Mayoz

Antifrajilitatearen batez bestekoa geruza bakarreko ASBentzat
pertubatutako X nodoen funtzioan
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
    
    N=50
    p=0.5
    T=100
    
    colors=['b', 'orange', 'g', 'brown', 'purple']
    
    number_of_iterations=50 
    fraction=1

    plt.ylabel("Antifrajilitatea")
    plt.xlabel("X")
    red=rbn.RBN()
    p1=np.zeros((5, N))
   
    for K in (1,6):
        f=np.zeros(( number_of_iterations, N ))
        i=0
        for x in trange(number_of_iterations):
            red.CreateNetScaleFree(K, N, p)
            f[i]=red.antifragile(T, O=1, runs=10, fraction=fraction) #X ardatzen X eukitzeko
            i+=1
        g1=np.mean(f, 0)
        #g1=np.insert(g1, 0,0)
        plt.errorbar(np.arange(1,N+1), g1, label="K= "+str(K), 
                     yerr=stats.sem(f), ecolor='r', color=colors[K-1])
       
  
    plt.title("Geruza bakarreko EAS")
    plt.legend()
    
    plt.savefig("Antifrajilitatea_X_ASB_Single_N="+str(N)+".pdf")
    
    
    plt.show()
#    plt.ylabel("Probability")
#    for K in range(1, 5):
#        plt.plot(np.arange(1,int(N/fraction)+1), p1[K-1], label="K= "+str(K))
#    plt.title("Probility of generating antifragile networks")
#    plt.legend()
    
#    red.CreateBioNet()
    
#    f=red.antifragile(100, runs=500, O=1, fraction=1)
#    plt.plot(f, label="BioNet")
#    
    
    print("--- %s seconds ---" % (time.time() - start_time))