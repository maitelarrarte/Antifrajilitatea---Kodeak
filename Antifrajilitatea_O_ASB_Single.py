# -*- coding: utf-8 -*-
"""
@author: Maite Larrarte Mayoz

Antifrajilitatearen batez bestekoa geruza bakarreko ASBentzat
O frekuentziaren funtzioan
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
    plt.xlabel("O")
    red=rbn.RBN()
    p1=np.zeros((5, 50))
    rango=np.arange(1,6)
    for K in rango:
        f=np.zeros(( number_of_iterations, 50 ))
        i=0
        for x in trange(number_of_iterations):
            red.CreateNet(int(K), N, p)
            f[i]=red.antifragile(T, X=25, runs=10, fraction=fraction) #X ardatzen O eukitzeko
            i+=1
        g1=np.mean(f, 0)
        p1[K-1]=np.sum(np.array(f) < 0, axis=0)/number_of_iterations #PROBABILITATEA
        plt.errorbar(np.arange(1,50+1), g1, label="K= "+str(K), 
                     yerr=stats.sem(f), ecolor='r', color=colors[K-1])

    plt.title("Geruza bakarreko ASB")
    plt.legend()
    

    plt.savefig("Antifrajilitatea_O_ASB_Single_N="+str(N)+".pdf") 
    
    
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