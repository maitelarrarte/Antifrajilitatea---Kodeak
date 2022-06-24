# -*- coding: utf-8 -*-
"""
@author: Maite Larrarte Mayoz

Geruza bakarreko ASBetarako perturbatutako 
sarearen konplexutasuna X-ren funtzioan
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
  
    plt.ylabel("C")
    plt.xlabel("X")
    red=rbn.RBN()
    p1=np.zeros((5, int(N/fraction)))
    
    for K in range(1, 6):
        C=np.zeros(( number_of_iterations, int(N/fraction) ))
        i=0
        for x in trange(number_of_iterations):
            red.CreateNet(int(K), N, p)
            for i in (np.arange(1,N+1)):
                State=red.RunNet(2*T,X=i,O=1)
                C[x,i-1]=rbn.complexity(State[-T:]) #X ardatzen X eukitzeko
            #f[i]=red.antifragile(T, X=40, runs=10, fraction=fraction) #X ardatzen O eukitzeko
            i+=1
        g1=np.mean(C, 0)
       
        plt.errorbar(np.arange(1,int(N/fraction)+1), g1, label="K= "+str(K), 
                     yerr=stats.sem(C), ecolor='r', color=colors[K-1])
    
    plt.title("Amaierako konplexutasuna geruza bakarreko ASBetarako")
    plt.legend()
    plt.savefig("AmaierakoKonplexutasuna_N="+str(N)+"_SINGLE.pdf")
    plt.show()

    
    print("--- %s seconds ---" % (time.time() - start_time))