# -*- coding: utf-8 -*-
"""

@author: Maite Larrarte Mayoz

Geruza bakarreko ASBetarako antifrajilitatearen minimoa, 
erakartzaie kopurua eta erakartzaileen luzera K-ren funtzioan.
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

    plt.xlabel("K")
    red=rbn.RBN()
     
    red=rbn.RBN()
    g1=np.zeros(5)
    A1=np.zeros(5)
    Al1=np.zeros(5)
    O=1
    X=4
    for K in trange(1,6):
        vf = []
        vp = []
        numAtt=[]
        avLen0=[]
        avLen=[]
        #fB=np.zeros( number_of_iterations)
        #f=np.zeros( number_of_iterations)
        Abb=np.zeros(number_of_iterations)
        Al=np.zeros(number_of_iterations)
        i=0
        #hurrengon 80
        for x in range(number_of_iterations):
            red.CreateNet(int(K), N, p)
            A0 = red.Attractors(T,runs=1000)
            print("A0=",len(A0))
            #numAtt.append(len(A0)) #Atraktore kopurua (ziklo kopurua) runs bakoitzeako jasotzeu
            if(len(A0) == 0):
                avLen0=0
            else:
                edos=0
                for j in A0:
                    edos+=len(j)
                avLen0=(edos/len(A0))
            #f[i]=red.antifragile(T=T,X=X, O=1, runs=10) #X ardatzen X eukitzeko
       
            #fB[i]=np.amin(f) # f-ren runs bakoitzeko balio minimoa jasotzen jutea
            #X=np.argmin(f)
            A=red.AttractorsRBN(T=T,X=X, O=1, runs=1000)
           
            Abb[i]=len(A)-len(A0)
            if(len(A) == 0):
                    avLen=100
            else:
                    edos=0
                    for j in A:
                        edos+=len(j)
                    avLen=(edos/len(A))
           
            Al[i]=avLen-avLen0
            i+=1
                     
        #g1[K-1]=np.mean(f)
       
  
        A1[K-1]=np.mean(Abb)
        
        Al1[K-1]=np.mean(Al)
        
        
    
    #plt.errorbar([1,2,3,4,5], g1, label="antizaurgarritasuna", 
    #                 yerr=stats.sem(f), ecolor='r', color=colors[1])
    plt.errorbar([1,2,3,4,5], A1, label="erakartzaile kopuruaren diferentzia", 
                     yerr=stats.sem(A1), ecolor='r', color=colors[2])
    plt.errorbar([1,2,3,4,5], Al1, label="erakartzaile luzeraren diferentzia", 
                     yerr=stats.sem(Al1), ecolor='r', color=colors[3])
        
      
    plt.title("Geruza bakarreko ASB")
    plt.legend()
    
    plt.savefig('AntifrajilitateaErakartzaileak_K_N='+str(N)+'_T='+str(T)+'_noi='+str(number_of_iterations)+'_X='+str(X)+'.pdf')   #X ardatzen X eukitzeko
 
    
    
    plt.show()
    
  