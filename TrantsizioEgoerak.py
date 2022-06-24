# -*- coding: utf-8 -*-
"""

@author: okarim

Trantsizio-egoerak adierazten dituzten diagrama eskematikoa

"""

import rbn
import numpy as np
import time
import matplotlib.pyplot as plt


if __name__ == '__main__':
    __spec__ = "ModuleSpec(name='builtins', loader=<class '_frozen_importlib.BuiltinImporter'>)"
    
    start_time = time.time()
    
    K=2
    N=20
    p=0.5
    T=20
    X=4
    O=1
    
    red=rbn.RBN()
    red.CreateNet(K, N, p)
    
#ATRAKTOREK: 
#    A=red.Attractors(T, runs=10) #runs=1000
#    print("\nAttractores: ")
#    print(len(A)) # A bektore bat da, elementu bakoitzak atraktorei dagokion batez besteko luzera emateu
#    edos=0
#    for i in A:
#        edos+=len(i)
#    print("Longitud promedio de Attractores: ")
#    print(edos/len(A))
#    if(edos!= 0):
#        print(str(len(A)/(edos)*100)+"%")
    
    initial=np.random.randint(0, 2, N)
    
    State=red.RunNet(2*T, initial)
    fig= plt.figure()
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.xlabel('Nodoa')
    plt.ylabel('Denbora')
    plt.title("Pertubatu gabe")
    
    plt.show()
    
    C0=rbn.complexity(State[-T:])
    
    print("C0",C0)
    
    State=red.RunNet(2*T, initial, X=X, O=O)
    plt.imshow(State, cmap='Greys', interpolation='None')
    plt.xlabel('Nodoa')
    plt.ylabel('Denbora')
    plt.title("Perturbatuta")
    #plt.savefig("Figures/Figure_1b2.eps")
    plt.show()
    
    C=rbn.complexity(State[-T:])
    
    print("C",C)
    
    print("Antifrajilitatea",rbn.fragility(C,C0,X,O,N,T))
        
    print("--- %s seconds ---" % (time.time() - start_time))