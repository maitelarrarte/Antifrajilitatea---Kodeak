# -*- coding: utf-8 -*-
"""

@author: Maite Larrarte Mayoz


Perturbatu aurreko eta ondorengo erakartzaile kopuruaren diferentzia
X perturbatutako nodo kopuruaren funtzioan
sare biologiko ezberdinetarako

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
    
    T=100
    colors=['b', 'orange', 'g', 'brown', 'purple', 'pink', 'y' ]
    number_of_iterations=50 #50
    #plt.ylabel("Antizaurgarritasuna")
    plt.xlabel("X")
    red=rbn.RBN()
   
       
    red.CreateBioNet(b=1)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=1)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="CD4+ T", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[0])
    
    
    red.CreateBioNet(b=2)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=2)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="Ugaztunak", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[1])
    
     
    
    red.CreateBioNet(b=3)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=3)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="Bihotza", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[2])
    


    red.CreateBioNet(b=4)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=4)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="Mikrobioma", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[3])



    red.CreateBioNet(b=5)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=5)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="Heriotza", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[4])

    
    
    red.CreateBioNet(b=6)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=6)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="A. thaliana", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[5])
    
    
    
    red.CreateBioNet(b=7)
    Abb=np.zeros((number_of_iterations,red.N))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=7)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(1,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-1]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(1,red.N+1), A1, label="Tumorea", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[6])
    
    
    plt.legend()

    plt.title("Erakartzaile kopuruaren diferentzia")
        
    
    plt.savefig("ErakartzaileDiferentzia_X_T="+str(T)+"number_of_it="+str(number_of_iterations)+"BIO.pdf")   #X ardatzen X eukitzeko

    plt.show()
    