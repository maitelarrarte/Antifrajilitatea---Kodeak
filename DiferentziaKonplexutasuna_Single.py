# -*- coding: utf-8 -*-
"""
@author: Maite Larrarte Mayoz

Konplexutasunaren aldaketa 
geruza bakarreko ASB
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
    
    N=50 
    p=0.5
    T=100 
    
    colors=['b', 'orange', 'g', 'brown', 'purple' ] 
    number_of_iterations = 50 
    fraction=1
  
    plt.ylabel("C - C0")

    plt.xlabel("X")  

    red=rbn.RBN()
    #p1=np.zeros((5, int(N/fraction))) #PROBABILITATEA
    red=rbn.RBN()
    
    for K in (1,6):
        CS0=np.zeros((number_of_iterations, int(N/fraction) )) #hasierako konplexitatea
        CS=np.zeros((number_of_iterations, int(N/fraction) )) #amaierako konplexitatea
        Chbi0=np.zeros((number_of_iterations,int(N/fraction))) 
        Chbi=np.zeros((number_of_iterations,int(N/fraction))) 
        for l in trange(number_of_iterations): #interakzio bakoitzan barrun 10 hasierako egoera aztertu ber dia
            runs=10 #zenbat hasierako egoera aztertuko dian
            Ch0=np.zeros((runs,int(N/fraction))) #hasierako konplexitateantzat
            Ch=np.zeros((runs,int(N/fraction))) #maiaerako konplexitateantzat
            for h in range(runs): #10 hasierako egoera
                red.CreateNet(K, N, p)
                initial=np.random.randint(0, 2, N)
                State0=red.RunNet(2*T,initial) 
                Ch0[h,:]=rbn.complexity(State0)
                for X in range(1, int(N/fraction)+1): #hasierako egoera bakoitzeako X guztitatko konplexitatea aztertukoa
                    #### Perturbatutakoa###
                    State=red.RunNet(2*T,initial, X, O=1) 
                    Ch[h,X-1]=rbn.complexity(State) #X ardatzen X eukitzeko
            Chb0=np.mean(Ch0,0)
            Chb=np.mean(Ch,0) #X bakoitzeako 10 hasierako egoera desberdinen arteko konplexitatean batez bestekoa
            Chbi0[l]=Chb0
            Chbi[l]=Chb #interakzio bakoitzeko X bakoitzeko konplexitatea
        
     
        g1=np.mean(Chbi,0)-np.mean(Chbi0,0)
        print("g1=",g1,np.size(g1))
        
        f= Chbi - Chbi0 #X bakoitzeko interakzio guztin batez bestekoa, hau da ploteau ber dana
        print("g1=",f,np.size(f))
        #g1=np.mean(f, 0)
       # p1[K-1]=np.sum(np.array(f) < 0, axis=0)/number_of_iterations #PROBABILITATEA
        #g1=np.insert(g1, 0,0)
        plt.errorbar(np.arange(1,N+1), g1, label="K= "+str(K), 
                     yerr=stats.sem(f,0), ecolor='r', color=colors[K-1])
      

    plt.title("Konplexutasun diferenetzia geruza bakarreko ASBetarako")
    plt.legend()
    plt.savefig("DiferentziaKonplexitatea_N="+str(N)+"single.pdf")

    
    plt.show()

    
    print("--- %s seconds ---" % (time.time() - start_time))