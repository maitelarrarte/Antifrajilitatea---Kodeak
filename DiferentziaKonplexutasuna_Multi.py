# -*- coding: utf-8 -*-
"""

@author: Maite Larrarte Mayoz

Konplexutasunaren aldaketa
geruza anitzeko ASB
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
    Nc=5
    Ng=4
    p=0.5
    T=100 
    
    colors=['b', 'orange', 'g', 'brown', 'purple']
    
    number_of_iterations = 10 #50
    fraction=1

    plt.ylabel("C - C0")
    plt.xlabel("X")  
    
    red=rbn.RBN()
    for K in trange(1, 6):
                
        CS0=np.zeros((number_of_iterations, N*Nc )) #hasierako konplexitatea
        CS=np.zeros((number_of_iterations, N*Nc )) #amaierako konplexitatea
        Chbi0=np.zeros((number_of_iterations, N*Nc)) 
        Chbi=np.zeros((number_of_iterations, N*Nc)) 
        
        
        for l in trange(number_of_iterations): #interakzio bakoitzan barrun 10 hasierako egoera aztertu ber dia
            
    
            runs=10 #zenbat hasierako egoera aztertuko dian
            Ch0=np.zeros((runs,N*Nc)) #hasierako konplexitateantzat
            Ch=np.zeros((runs,N*Nc)) #maiaerako konplexitateantzat
            
            for h in range(runs): #10 hasierako egoera
                red.CreateNetMultilayer(K, N, p, Ng, Nc)
                Cell_initial=[]
                for j in range(Nc):
                    initial = np.random.randint(0, 2, N)
                    Cell_initial.append(initial.tolist())
                
                for X in range(1, (N*Nc)+1): #hasierako egoera bakoitzeako X guztitatko konplexitatea aztertukoa
                    ####Perturbatu gabekoa###
                    Cell_State0=red.RunNetMultilayer(2*T,Cell_initial=Cell_initial) 
                    for j in range(Nc):
                        State0=np.asarray(Cell_State0[j])
                        State0=State0[-T:,:]
                        Cell_State0[j]=State0.tolist()
                         
                    Ch0[h,X-1]=rbn.complexityMultilayer(Cell_State0) #X ardatzen X eukitzeko
                    #### Perturbatutakoa###
                    Cell_State=red.RunNetMultilayer(2*T,Cell_initial=Cell_initial, X=X, O=1) 
                    for j in range(Nc):
                        State=np.asarray(Cell_State[j])
                        State=State[-T:,:]
                        Cell_State[j]=State.tolist()
                    Ch[h,X-1]=rbn.complexityMultilayer(Cell_State) #X ardatzen X eukitzeko
            Chb0=np.mean(Ch0,0)
            Chb=np.mean(Ch,0) #X bakoitzeako 10 hasierako egoera desberdinen arteko konplexitatean batez bestekoa
            Chbi0[l]=Chb0
            Chbi[l]=Chb #interakzio bakoitzeko X bakoitzeko konplexitatea
        
        g1=np.mean(Chbi,0)-np.mean(Chbi0,0) #X bakoitzeko interakzio guztin batez bestekoa, hau da ploteau ber dana
         
       
        plt.errorbar(np.arange(1,(N*Nc)+1), g1, label="K= "+str(K),
                     yerr=stats.sem(Chbi), ecolor='r', color=colors[K-1])

    plt.title("Konplexutasun diferentzia geruza anitzeko ASBetarako")
    plt.legend()
    
    
    plt.savefig("DifrentziaKonplexitatea_N="+str(N)+"_NC="+str(Nc)+"Ng="+str(Ng)+"multi.pdf")
    plt.show()

  
    
    print("--- %s seconds ---" % (time.time() - start_time))