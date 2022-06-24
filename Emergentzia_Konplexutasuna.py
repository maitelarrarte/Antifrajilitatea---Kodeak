# -*- coding: utf-8 -*-
"""

@author: Maite Larrarte Mayoz

Emergentzia, bere osagarria eta konplexutasunaren arteko erlazioa
"""

import rbn
import numpy as np
from math import log
import matplotlib.pyplot as plt


colors=['b', 'orange', 'g', 'brown', 'purple']
rango=np.arange(0.0,1.001,0.001)
E=np.zeros(np.size(rango))
C=np.zeros(np.size(rango))
EE=np.zeros(np.size(rango))
i=0

plt.xlabel("po(p1)")
for p0 in rango:
    p1=1-p0
    if (p0==0.0)==True:
            p0=1.0
    if (p1==0.0)==True:
            p1=1.0
    E[i]=-(p0*log(p0,2)+p1*log(p1,2))
    EE[i]=1-E[i]
    C[i]=4*E[i]*(1-E[i])
    i+=1
    
plt.errorbar(rango, E, label="E", color=colors[0])
       
plt.errorbar(rango, EE, label="1-E", color=colors[1])
plt.errorbar(rango, C, label="C", color=colors[2])  

plt.legend()
    
plt.savefig("Emergentzia.pdf")
plt.show()