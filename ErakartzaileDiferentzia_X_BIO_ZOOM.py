# -*- coding: utf-8 -*-
"""

@author: Maite Larrarte Mayoz

Perturbatu aurreko eta ondorengo erakartzaile kopuruaren diferentzia
X perturbatutako nodo kopuruaren funtzioan
sare biologiko batzuetarako "zoom" eginda
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
    
    
    red.CreateBioNet(b=5)
    Abb=np.zeros((number_of_iterations,red.N-4))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=5)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(5,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-5]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(5,red.N+1), A1, label="Heriotza", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[4])
    
    
    red.CreateBioNet(b=7)
    Abb=np.zeros((number_of_iterations,red.N-4))
    for x in trange(number_of_iterations):
        red.CreateBioNet(b=7)
        A0 = red.Attractors(T,runs=1000)
        for XX in range(5,red.N+1):
            A=red.AttractorsRBN(T,X=XX, O=1, runs=1000)
            Abb[x,XX-5]=len(A)-len(A0)
                          
    A1=np.mean(Abb,0)
    plt.errorbar(np.arange(5,red.N+1), A1, label="Tumorea", 
                     yerr=stats.sem(Abb), ecolor='r', color=colors[6])
    
    
    plt.legend()

    plt.title("Erakartzaile kopuruaren diferentzia")
        
    
    plt.savefig("ErakartzaileDiferentzia_X_T="+str(T)+"number_of_it="+str(number_of_iterations)+"BIO_ZOOM.pdf")   #X ardatzen X eukitzeko
    #plt.savefig("Figure_3b.eps") #X ardatzen O eukitzeko
    plt.show()