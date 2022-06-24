# -*- coding: utf-8 -*-
"""
@author: okarim, Maite Larrarte Mayoz

Gainontzeko kodeetako funtzioak
"""
import time
import numpy as np
import multiprocessing
from functools import partial
import matplotlib.pyplot as plt
import pandas as pd
from  tqdm import trange
from math import log

class RBN:
    #OKARIM
    def CreateNet(self, K, N, p):
        """
        K = number of connections
        N = number of nodes, indexed 0 .. N-1
        p = probability of one
        """
        self.K=K
        self.N=N
        if(type(self.K) is int):
             self.Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:self.K]
             #self.Bool = np.random.choice([0, 1], size=(N, 2**self.K), p=[1-p, p]) # N random Boolean functions, a list of 2^K ones and zeros.
             self.Bool = np.random.randint(0, 2, size=(N, 2**self.K))  # N random boolean functions, a list of 2^k  ones and zeros.
        else:
            Kv=np.random.poisson(self.K, N)
            print("KV=",Kv)
            Kv=Kv/(np.mean(Kv)/K)
            print("KV=",Kv)
            Kv=np.ceil(Kv)
            print("KV=",Kv)
            Kv=Kv.astype(int)
            Kv[np.where(Kv>N)]=N
            maximo=np.amax(Kv)
            
            self.Con=np.zeros((N+1, maximo),dtype=int)
            self.Bool=np.zeros((N+1, 2**maximo),dtype=int)
            
            for i in range(N):
                self.Con[i+1, 0:Kv[i]] = np.random.choice(N, Kv[i], replace=False)+1
                self.Bool[i+1, 0:2**Kv[i]] = (np.random.choice([0, 1], size=2**Kv[i], p=[1-p, p]))
            print("CON=",self.Con)
        return
   
    #MAITE
    def CreateNetM(self, K, N, p):
        """
        K = number of connections
        N = number of nodes, indexed 0 .. N-1
        p = probability of one
        """
    
        if(type(self.K) is int):
             self.Con = np.apply_along_axis(np.random.permutation, 1, np.tile(range(N), (N,1) ))[:, 0:self.K]
             #self.Bool = np.random.choice([0, 1], size=(N, 2**self.K), p=[1-p, p]) # N random Boolean functions, a list of 2^K ones and zeros.
             self.Bool = np.random.randint(0, 2, size=(N, 2**self.K))  # N random boolean functions, a list of 2^k  ones and zeros.
        else:
            Kv=np.random.poisson(self.K, N)
            Kv=Kv/(np.mean(Kv)/K)
            Kv=np.ceil(Kv)
            Kv=Kv.astype(int)
            Kv[np.where(Kv>N)]=N
            maximo=np.amax(Kv)
            
            self.Con=np.zeros((N+1, maximo),dtype=int)
            self.Bool=np.zeros((N+1, 2**maximo),dtype=int)
            
            for i in range(N):
                self.Con[i+1, 0:Kv[i]] = np.random.choice(N, Kv[i], replace=False)+1
                self.Bool[i+1, 0:2**Kv[i]] = (np.random.choice([0, 1], size=2**Kv[i], p=[1-p, p]))
        return
    #MAITE
    def CreateNetMultilayer(self,K,N,p,Ng,Nc,L=0,M=False):
        """
        K = number of connections
        N = number of nodes, indexed 0 .. N-1
        Nc = zelula kopurua
        Ng = gene komunikatzallen kopurua
        p = probability of one
        L = number of links of an inter-network
        M = True --> Mantendu intrazelulen topologia
        Cell_Bool = Funtzio boolearren taulak jasotzettun lista zelula/sare/RBN bakoitzeako
        Cell_Con = Konekxion taula jasotzeu zelula/sare/RBN bakoitzeako
        """
        self.K=K
        self.N=N
        self.p=p
        self.Ng=Ng
        self.Nc=Nc
        #self.L=L
    
        Cell_C=[]
        Cell_B=[]
        L_bekt=np.zeros(Nc,dtype=int) #zelula bakoitza beste zenbatekin konektatukoan finkatzeu
        
        if (L==0)==True:
            L=np.random.randint(1,Nc*Nc) #L number of links beste RBNkin, gene komunikatzallen konekxio kopurua
           
        L_bekt[L-len(L_bekt)*int(L/len(L_bekt)):]=int(L/len(L_bekt))    
        L_bekt[:L-len(L_bekt)*int(L/len(L_bekt))]=int(L/len(L_bekt))+1
                
         
            
        #### L aleatorio zelula bakoitzeako ####
        #L_bekt=np.random.choice(np.arange(Nc),Nc)
        #while (np.sum(L_bekt)==L)==False:    
        #        if (np.sum(L_bekt)>L)==True: 
        #            indizek=np.array(np.where(L_bekt!=0))
        #            hautazkoa=np.random.choice(indizek[0],1)
        #            L_bekt[hautazkoa[0]]=L_bekt[hautazkoa[0]]-1
        #        elif (np.sum(L_bekt)<L)==True:
        #            indizek=np.array(np.where(L_bekt < self.Nc))
        #            hautazkoa=np.random.choice(indizek[0],1)
        #            L_bekt[hautazkoa]=L_bekt[hautazkoa]+1
          
        if (M==True):
            for i in range(Nc):
                self.AddInterCellConnections(L_bekt[i])
            
        else:
            for i in range(Nc):
                self.CreateNetM(K,N,p)
                A=np.zeros((self.N,self.Nc),dtype=int) #kontuz hemen: K<Nc kastako baliou
                A[:,:K]=self.Con
                self.Con=A
                self.AddInterCellConnections(L_bekt[i])
                Cell_C.append((self.Con).tolist())
                self.Cell_Con=Cell_C
                Cell_B.append((self.Bool).tolist())
                self.Cell_Bool=Cell_B
        
        return   
   
    #OKARIM
    def CreateBioNet(self, b=1): #KASU BAKOITZAN DAGOZKION TAULAK IRAKURRI TA KONEKXION MATRIZEA, FUNTZIO BOOLEARRAN MATRIZEA TA BATEZ BESTEKO K LORTU
        if(b==1):
            self.N=18
            
            #data=pd.read_csv("Con_CD4+TCELL.csv", sep=",", header=-1)
            data=pd.read_csv("Con_CD4+TCELL.csv", sep=",", header=None)
            data=data.fillna(0)
            #x=np.array(data.loc[:, 1:]).astype(int)
            x=np.array(data.iloc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x]) #CON taula matrize bihurtzeu izenburuk 0-tzat hartuz
            
            data=pd.read_csv("Bool_CD4+TCELL.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y]) #BOOL taula matrize bihurtzeu izenburuk 0-tzat hartuz
           
            self.e=6 #Con-en konektatu gabeko nodok
            self.K=np.count_nonzero(self.Con)/(self.N-self.e) #konekxiok zati konektautako nodo kopurua --> K=nodo bakoitzeko konexio kopurun batez bestekoa
    
        elif(b==2):
            self.N=20
            #data=pd.read_csv("Con_Mammalian Cell Cycle.csv", sep=",", header=-1)
            data=pd.read_csv("Con_Mammalian Cell Cycle.csv", sep=",", header=None)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
        
            
            data=pd.read_csv("Bool_Mammalian Cell Cycle.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=1
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
           
            print("K=", self.K)
            
        elif(b==3):
            self.N=15
            #data=pd.read_csv("Con_Cardiac development.csv", sep=",", header=-1)
            data=pd.read_csv("Con_Cardiac development.csv", sep=",", header=None)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Cardiac development.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=1
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==4):
            self.N=12
            #data=pd.read_csv("Con_Metabolic.csv", sep=",", header=-1)
            data=pd.read_csv("Con_Metabolic.csv", sep=",", header=None)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            #print(self.Con)
            data=pd.read_csv("Bool_Metabolic.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=1
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==5):
            self.N=28
            #data=pd.read_csv("Con_Death.csv", sep=",", header=-1)
            data=pd.read_csv("Con_Death.csv", sep=",", header=None)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Death.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=3
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==6):
            self.N=14
            
            #data=pd.read_csv("Con_Arabidopsis.csv", sep=",", header=-1)
            data=pd.read_csv("Con_Arabidopsis.csv", sep=",", header=None)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Arabidopsis.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=0
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
            
        elif(b==7):
            self.N=32
            data=pd.read_csv("Con_Tumour Cell.csv", sep=",", header=None)
            #data=pd.read_csv("Con_Tumour Cell.csv", sep=",", header=-1)
            data=data.fillna(0)
            x=np.array(data.loc[:, 1:]).astype(int)
            self.Con=np.vstack([np.zeros(x.shape[1],dtype=int), x])
            
            data=pd.read_csv("Bool_Tumour Cell.csv", sep=",", header=0)
            data=data.fillna(0)
            y=np.array(data.loc[:, :]).astype(int).transpose()
            self.Bool=np.vstack([np.zeros(y.shape[1],dtype=int), y])
            self.e=2
            self.K=np.count_nonzero(self.Con)/(self.N-self.e)
        
        return 

    #MAITE
    def CreateBioNetMultilayer(self, b, Nc, L=0):
        self.Nc=Nc
        self.Cell_Bool=[]
        self.Cell_Con=[]
  
        L_bekt=np.zeros(Nc)
        if (L==0)==True:
            L=np.random.randint(1,Nc*Nc)
        print("L=",L)    
        L_bekt[L-len(L_bekt)*int(L/len(L_bekt)):]=int(L/L_bekt)    
        L_bekt[:L-len(L_bekt)*int(L/len(L_bekt))]=int(L/len(L_bekt))+1
        print("L_bekt=",L_bekt)
        
        
        #L aleatorio zelula bakoitzeako ###
        
        #L_bekt=np.random.choice(np.arange(Nc),Nc)
        
       
        
        for i in range(Nc):
            self.CreateBioNet(b)
            self.AddInterCellConnections(L_bekt[i])
            self.Cell_Con.append(self.Con) 
            self.Cell_Bool.append(self.Bool)
          
        return   
    #MAITE
    def AddInterCellConnections(self,L): #L number of links
        """ 
        Zelulan arteko interakziok finkatzettu
        """
        if (type(self.K) is int) :
            self.Con[self.N-self.Ng:]=0
            Con_connections=self.Con[self.N-self.Ng:]
            indizek=np.where(Con_connections==0)
            indizen_matrizea=np.array(indizek)
            indizea = np.empty(L,dtype=int)
            
            indizea=indizen_matrizea[:,np.random.choice(np.size(indizen_matrizea,1),L,replace=False)]
            
            for i in range(L):
                ind=indizea[:,i]
                Con_connections[ind[0],ind[1]]=np.random.choice(np.arange(self.N-self.Ng,self.N),1)
            self.Con[self.N-self.Ng:]=Con_connections
        
        else:
            self.Con[self.N-self.e+1:,:]=0
            Con_connections=self.Con[self.N-self.e+1:,:]
            
            indizek=np.where(Con_connections==0)
            indizen_matrizea=np.array(indizek)
            
            indizea = np.empty(L,dtype=int)
            indizea=indizen_matrizea[:,np.random.choice(np.size(indizen_matrizea,1),L,replace=False)]
            
            for i in range(L):
                ind=indizea[:,i]
                Con_connections[ind[0],ind[1]]=np.random.choice(np.arange(self.N-self.e+1,self.N+1),1)
            self.Con[self.N-self.e+1:]=Con_connections

            return 
    #MAITE
    def AddInterCellConnections2(self,L): #L number of links
            """ Zelulan arteko interakziok finkatzettu baño kasu hontan konekxioko nodok beste edozein nodokin konektatzettu"""
        
#        if (type(self.K) is int) :
#            Con_connections=self.Con[self.N-self.Ng:]
#            indizek=np.where(Con_connections==0)
#            indizen_matrizea=np.array(indizek)
#            indizea = np.empty(L,dtype=int)
#            indizea=indizen_matrizea[:,np.random.choice(np.size(indizen_matrizea,1),L)]
            
#            for i in range(L):
#                ind=indizea[:,i]
#                Con_connections[ind[0],ind[1]]=np.random.choice(np.arange(self.N-self.Ng,self.N+1),1)
#            self.Con[self.N-self.Nc:]=Con_connections
        
#        else:
            Con_connections=self.Con[self.N-self.e+1:]
            
            indizek=np.where(Con_connections==0)
            indizen_matrizea=np.array(indizek)
            
            indizea = np.empty(L,dtype=int)
            indizea=indizen_matrizea[:,np.random.choice(np.size(indizen_matrizea,1),L)]
            
            for i in range(L):
                ind=indizea[:,i]
                Con_connections[ind[0],ind[1]]=np.random.choice(np.arange(1,self.N+1),1)
            self.Con[self.N-self.e+1:]=Con_connections

            return 
    #OKARIM
    def RunNet(self, T, initial=[], X=0, O=0, Bio=False):
        """
        Con= matrix of connections
        Bool= lookup table
        T = timesteps
        initial = initial state (random if empty)
        X = how many perturbations
        O = how often the perturbations take place
        """
        Pow = 2**np.arange(np.size(self.Con, 1)) # [ 1 2 4 ... ], for converting inputs to numerical value
        
        if(type(self.K) is int):
            a=0
            State = np.zeros((T+1,self.N),dtype=int)
            if np.array_equal(initial, []):
                State[0] = np.random.randint(0, 2, self.N) 
            else:
                State[0] = initial
        else:
            a=1
            State = np.zeros((T+1,self.N+1),dtype=int) 
            if np.array_equal(initial, []):
                State[0] = np.append([0], np.random.randint(0, 2, self.N)) 
            else:
                State[0] = np.append([0],initial)
            
            self.Bool[np.where(self.Con[:,0]==0),0] = State[0, np.where(self.Con[:,0]==0)] # if node doesn't have conections not change 
           
        #print("state0=",State[0], len(State[0]))
        #print("BOOL=",self.Bool, np.size(self.Bool,0))
        #print("CON=",self.Con, np.size(self.Con,0))
        for t in range(T):  # 0 .. T-1
                self.Bool[np.where(self.Con[:,0]==0),0] = State[t, np.where(self.Con[:,0]==0)] #konekxio gabeko nodoi hasierako egoera finkatzeie funtzio boolearren
                State[t+1] = self.Bool[:, np.sum(Pow * State[t,self.Con],1)].diagonal()
                if ( X and O ) != 0:  #Perturbations 
                    if t%O == 0:
                        State[t+1,  np.random.choice(self.N, size=X, replace=False)+a] = np.random.randint(0, 2, X)
                if(Bio and self.e>0):
                    State[t+1, self.N-self.e+1:]=np.random.randint(0, 2, self.e) #koxekcio gabeko nodoi balio eman
        if(Bio and self.e>0):
            return(State[:,1:-self.e])
        elif(type(self.K) is int):
            return(State)
        else:
            return(State[:,1:])
    
    #MAITE
    def RunNetMultilayer(self, T, b=0, Cell_initial=[], L=0, X=0, O=0):
        """
        Cell_Con= list of matrix of connections (zelula bakoitzeako listako elementu bat)
        Cell_Bool= lookup table (zelula bakoitzana listako elementu bakoitzen)
        T = timesteps
        Nc = Number of cells
        initial = initial state (random if empty)
        X = how many perturbations
        O = how often the perturbations take place
        """
        
        
        Cell_State=[]
        
        if(type(self.K) is int):
            a=0
            for i in range(self.Nc):
                State0 = np.zeros((T+1,self.N),dtype=int)
                Cell_State.append(State0.tolist())

            if np.array_equal(Cell_initial, []):
                for i in range(self.Nc):
                    State=np.asarray(Cell_State[i])
                    State[0] = np.random.randint(0, 2, self.N)
                    Cell_State[i]=State.tolist()
            else:
                for i in range(self.Nc):
                    State=np.asarray(Cell_State[i])
                    initial=np.asarray(Cell_initial[i])
                    State[0] = initial
                    Cell_State[i]=State.tolist()          
            
        else:
            a=1
            #Cell_X=np.zeros(self.Nc)
            #Cell_T_zati_O=np.zeros(self.Nc)
            #O_batezbesteko=np.zeros(self.Nc,dtype=int)
            for i in range(self.Nc):
                State0=np.zeros((T+1,self.N+1),dtype=int)
                Cell_State.append(State0.tolist())

        #Zelula adina State eukikottut
        
            if np.array_equal(Cell_initial, []):
                for i in range(self.Nc):
                    State=np.asarray(Cell_State[i])
                    State[0] = np.append([0], np.random.randint(0, 2, self.N)) 
                    Cell_State[i]=State.tolist()
            else:
                for i in range(self.Nc):
                    State=np.asarray(Cell_State[i])
                    initial=np.asarray(Cell_initial[i])
                    State[0] = np.append([0],initial)
                    Cell_State[i]=State.tolist()  
        
            for j in range(self.Nc):
                Bool=np.asarray(self.Cell_Bool[j])
                Con=np.asarray(self.Cell_Con[j])
                State=np.asarray(Cell_State[j])
                
                Bool[np.where(Con[:,0]==0),0] = State[0, np.where(Con[:,0]==0)] # if node doesn't have conections not change 
                
                self.Cell_Bool[j]=Bool.tolist()
                self.Cell_Con[j]=Con.tolist()
                Cell_State[j]=State.tolist()
                
        for t in range(T):  # 0 .. T-1
                #INTRACELULLAR INTERACTIONS Nodoak 0tik (N-e) arte
                for j in range(self.Nc):
                    Bool=np.asarray(self.Cell_Bool[j])
                    Con=np.asarray(self.Cell_Con[j])
                    State=np.asarray(Cell_State[j])
                    if(type(self.K) is int):
                        Con0=Con[:,:self.K]
                        State0=State 
                    Pow = 2**np.arange(np.size(Con0, 1)) # [ 1 2 4 ... ], for converting inputs to numerical value 
                    #print("t",t)
                    #print("State99=",State[t])
                    #print("Con",Con)
                    
                    #print("State[t,con]=",State[t,Con])
                    #print("Pow",Pow)
                    #print("Pow * State[t,Con]",Pow * State[t,Con])
                    #print("np.sum(Pow * State[t,Con],1)",np.sum(Pow * State[t,Con],1))
                    
                    
                    Bool[np.where(Con0[:,0]==0),0] = State0[t, np.where(Con0[:,0]==0)]
                    
                    State0[t+1] = Bool[:, np.sum(Pow * State0[t,Con0],1)].diagonal()
                
                    if(b>0)==True: #gene komunikatzailek utzi hasieran bezela, aurreago aldaukoia hauek UPDATE RULES DESBERDINAK DAZKATELAKO
                        State0[t+1, self.N-self.e+1:]=State0[t,self.N-self.e+1:]
                            #print(" State[t+1, self.N-self.e+1:]=",State[t+1, self.N-self.e+1:])
                    elif(type(self.K) is int):
                        State0[t+1,self.N-self.Ng:]=State0[t,self.N-self.Ng:]
                    Cell_State[j]=State0.tolist()
                
                #INTERCELLULAR INTERACTIONS    
                for j in range(self.Nc):
            
                    Con=np.asarray(self.Cell_Con[j])
                    State=np.asarray(Cell_State[j])
                    
                    if(type(self.K) is int):
                        Nodok=np.arange(self.Nc)
                        for h in range(self.N-self.Ng,self.N): 
                            activated=0
                            ind_matriz=np.array(np.where(Con[h,:]!=0))
                            if (np.size(ind_matriz) > 0 ) ==True:  
                                ind_bekt=ind_matriz[0]
                                for indizek in ind_bekt:
                                    aztertu_beharreko_nodoa=Con[h,indizek]
                                    aztertu_beharreko_sarea=np.random.choice(Nodok,1)
                                    borratzeko=np.array(np.where(Nodok==aztertu_beharreko_sarea))
                                    borratzeko_bekt=borratzeko[0]
                                    Nodok=np.delete(Nodok,borratzeko_bekt[0])
                                    
                                    #aztertu beharreko sarea hartu
                                    Stateinter=np.asarray(Cell_State[aztertu_beharreko_sarea[0]])
                                    if (Stateinter[t,aztertu_beharreko_nodoa]==1)==True:
                                        activated=1
                            if (activated==1)==True:
                                State[t+1,h]=1 #aurreko linketako bat aktibauta baldin badao, hau aktibatukoa
                            else:
                                State[t+1,h]=0
                    elif (b>0)==True:
                        Nodok=np.arange(self.Nc)
                        for h in range(self.N-self.e+1,self.N+1): 
                            activated=0
                            ind_matriz=np.array(np.where(Con[h,:]!=0))
                            if (np.size(ind_matriz) > 0 ) ==True:  
                                ind_bekt=ind_matriz[0]
                                for indizek in ind_bekt:
                                    aztertu_beharreko_nodoa=Con[h,indizek]
                                    aztertu_beharreko_sarea=np.random.choice(Nodok,1)
                                    borratzeko=np.array(np.where(Nodok==aztertu_beharreko_sarea))
                                    borratzeko_bekt=borratzeko[0]
                                    Nodok=np.delete(Nodok,borratzeko_bekt[0])
                                    
                                    #aztertu beharreko sarea hartu
                                    Stateinter=np.asarray(Cell_State[aztertu_beharreko_sarea[0]])
                                    if (Stateinter[t,aztertu_beharreko_nodoa]==1)==True:
                                        activated=1
                            if (activated==1)==True:
                                State[t+1,h]=1 #aurreko linketako bat aktibauta baldin badao, hau aktibatukoa
                            else:
                                State[t+1,h]=0
                    Cell_State[j]=State.tolist()
                    self.Cell_Con[j]=Con.tolist()
                
                #PERTURBAZIOK
                if ( X and O ) != 0:  #Perturbations 
                        if t%O == 0:
                            X_bekt=np.zeros(self.Nc,dtype=int) #sare bakoitzeko X gordekou
                      
                            X_bekt[X-len(X_bekt)*int(X/len(X_bekt)):]=int(X/len(X_bekt))    
                            X_bekt[:X-len(X_bekt)*int(X/len(X_bekt))]=int(X/len(X_bekt))+1
                           
                            ### X aleatorioko zelula bakoitzai ##
                            #X_bekt=np.random.choice(np.arange(1,self.N+1),self.Nc)
                    
                            #while (np.sum(X_bekt)==X)==False:    
                            #    if (np.sum(X_bekt)>X)==True: 
                            #        indizek=np.array(np.where(X_bekt!=0))
                            #        hautazkoa=np.random.choice(indizek[0],1)
                            #        X_bekt[hautazkoa[0]]=X_bekt[hautazkoa[0]]-1
                            #    elif (np.sum(X_bekt)<X)==True:
                            #        indizek=np.array(np.where(X_bekt < self.N))
                            #        hautazkoa=np.random.choice(indizek[0],1)
                            #        X_bekt[hautazkoa]=X_bekt[hautazkoa]+1
                                    
                            for i in range(self.Nc):
                                #Cell_X[i]+=X_bekt[i]
                                if (X_bekt[i]!=0)==True:
                                    #O_batezbesteko[i]+=1
                                    State=np.asarray(Cell_State[i])
                                    State[t+1,  np.random.choice(self.N, size=X_bekt[i], replace=False)+a] = np.random.randint(0, 2, X_bekt[i])
                                    Cell_State[i]=State.tolist()
                            
                        #ONDORENGO IF HONEKIN EZTAKIT ZE IN
                        #Perturbaziok gene komunikatzalletan:
                        #if(Bio and self.e>0):
                        #    State[t+1, N-self.e+1:]=np.random.randint(0, 2, self.e) #koxekcio gabeko nodoi balio eman   
                
                #iterazio bakoitzen gene komunikatzailen konekxiok aldau
                if(type(self.K) is int):
                    M=True
                    self.CreateNetMultilayer(self.K, self.N, self.p, self.Ng, self.Nc, L, M)
    
                elif (b>0)==True:
                    self.CreateBioNetMultilayer(b, self.Nc, L)
                else:
                    self.CreateNetMultilayerScaleFree(self.K, self.N, self.p, self.Ng, self.Nc, L, M)

            
           # Cell_X=Cell_X/O_batezbesteko
           # Cell_O=O_batezbesteko
           # Cell_T_zati_O=T/O_batezbesteko
                         
        if(b>0)==True:
                    for j in range(self.Nc):
                        State=np.asarray(Cell_State[j])
                        State1=State[1:]
                        Cell_State[j]=State1.tolist()
                    return Cell_State
        elif(type(self.K) is int):
                    return Cell_State
        else:
           for j in range(self.Nc):
                        State=np.asarray(Cell_State[j])
                        State1=State[1:]
                        Cell_State[j]=State1.tolist() 

    

    #OKARIM
    def Attractors(self, T, runs=0):
        """
        List of Attractors of R random initial states
        runs = number of runs (if 0 then List of Attractors of every possible initial state)
        T = timesteps
        """
        attractList=[]
        if (runs == 0)==True :
            for i in range(np.power(2,self.N)): #2^N # i=1,2,3,4,...,2^N
                initial=[x=='1' for x in format(i, '0'+str(self.N)+'b')] #b:binary format , 
                State=self.RunNet(T, initial)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0) #zutabe bakoitzen behin azaltzeian numerok kontatzettu 
                # unique elements: matrize guztin zehar behin bakarrik azaltze dian zenbakik
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        else:
            for i in range(runs):
                State=self.RunNet(2*T)
                unique_elements, counts_elements = np.unique(State[-T:], return_counts=True, axis=0)      
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())        
        return attractList
    
    #MAITE 
    def AttractorsMultilayer(self, T, b, L=0, runs=0):
        """
        List of Attractors of R random initial states
        runs = number of runs (if 0 then List of Attractors of every possible initial state)
        T = timesteps
        """
        attractList=[]
        if runs == 0 :
            for i in range(np.power(2,self.N)): #2^N # i=1,2,3,4,...,2^N
                initial=[x=='1' for x in format(i, '0'+str(self.N)+'b')] #b:binary format , 
                State=self.RunNet(T, initial)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0) #zutabe bakoitzen behin azaltzeian numerok kontatzettu 
                # unique elements: matrize guztin zehar behin bakarrik azaltze dian zenbakik
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        else:
            for i in range(runs):
                Cell_State=self.RunNetMultilayer(T, b, L=L)
                for j in range(self.Nc):
                    State=np.asarray(Cell_State[j])
                    print("State",State)
                    unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
                    A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion    
                    if not(A.tolist() in attractList):  #if A is not in attractList then add it
                        attractList.append(A.tolist())
                    
                
        return attractList
    #MAITE
    def AttractorsRBN(self, T, X=0, O=0, runs=0): #ERREPIKATZEIAN EGOERAK BUELTATZETTU
        """
        List of Attractors of R random initial states
        runs = number of runs (if 0 then List of Attractors of every possible initial state)
        T = timesteps
        """
        attractList=[] 
        
        if (runs==0)==True:
            for i in range(np.power(2,self.N)): #2^N # i=1,2,3,4,...,2^N
                initial=[x=='1' for x in format(i, '0'+str(self.N)+'b')] #b:binary format ,
                State=self.RunNet(T, initial, X, O)
                            
                unique_elements, counts_elements = np.unique(State[-T:], return_counts=True, axis=0) #zutabe bakoitzen behin azaltzeian numerok kontatzettu 
                # unique elements: matrize guztin zehar behin bakarrik azaltze dian zenbakik
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one o
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                        attractList.append(A.tolist())
        else :
            for i in range(runs):
                State=self.RunNet(2*T,X=X,O=O)
              
                unique_elements, counts_elements = np.unique(State[-T:], return_counts=True, axis=0) #zutabe bakoitzen behin azaltzeian numerok kontatzettu 
                # unique elements: matrize guztin zehar behin bakarrik azaltze dian zenbakik
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                        attractList.append(A.tolist())    
        return attractList
    #MAITE
    def AttractorsFlipped(self, T, X, O, runs=0):
        """
        List of Attractors of R random initial states
        runs = number of runs (if 0 then List of Attractors of every possible initial state)
        T = timesteps
        """
        attractList=[]
        if runs == 0 :
            for i in range(np.power(2,self.N)): #2^N # i=1,2,3,4,...,2^N
                initial=[x=='1' for x in format(i, '0'+str(self.N)+'b')] #b:binary format , 
                State=self.RunNet(T, initial, X=X, O=O)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0) #zutabe bakoitzen behin azaltzeian numerok kontatzettu 
                # unique elements: matrize guztin zehar behin bakarrik azaltze dian zenbakik
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
                
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        else:
            for i in range(runs):
                State=self.RunNet(T, X=X, O=O)
                unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)      
                A=unique_elements[np.where(counts_elements > 1)] #States that appear more than one occasion
    
                if not(A.tolist() in attractList):  #if A is not in attractList then add it
                    attractList.append(A.tolist())
                
        return attractList
    
       
    
    #OKARIM   
    def RBNSort(self): 
        """
        Sort the nodes by their overall activity
        """
        SRun = 5     # sorting runs
        ST = 200     # sorting timesteps
        Totals = np.zeros(self.N,dtype=int) # [0 0 0 0 .... 0]
        
        for r in range(SRun):
            State=self.RunNet(ST)  # Perturbatu gabeko sarea T=400 arte
            Totals = Totals + np.sum(State, 0) #Nodo bakoitza 1 egoeran zenbat aldiz eon dan kontatzeu T hoietako ta gañea SRun bakoitzeakore bai
        
        Index = np.argsort(Totals)    # permutation indexes for sorted order: Totalse-eko indizen baliok ematerru zenbakik txikinetik haundinea ordenauta
        
        if(type(self.K) is int):
            self.Bool = self.Bool[Index]         # permute the boolean functions: Index-ak ematettun indizen arabera berordenatzettu Bool matrizeko errenkadak
            self.Con = self.Con[Index]           # permute the connections
            
            InvIndex = np.argsort(Index)         # inverse permutation
            self.Con = InvIndex[self.Con]        # relabel the connections: nodo bakoitzeko konekxiok interkanbiatzettu
        else:
            self.Bool[1:] = self.Bool[Index+1]         # permute the boolean functions
            self.Con[1:] = self.Con[Index+1]           # permute the connections
            InvIndex = np.append([-1], np.argsort(Index)) # inverse permutation
            self.Con[1:] = InvIndex[self.Con[1:]]+1        # relabel the connections
        return
    
    #MAITE
    def AddLink(self):
        """
        Konekxio bat gehituko zaio aleatorioki
        """
        
        if (self.e==0)==True:
            a=np.where(self.Con[1:,:]==0)
        else:
            a=np.where(self.Con[1:-self.e,:]==0)
        
        if (np.size(a)>0)==True:
                ar=np.array(a)
                
                indizea = np.empty(2,dtype=int)
                
                indize=ar[:,np.random.choice(np.size(ar,1),1)]
                b=indize[0]
                c=indize[1]
                indizea[0]=b[0]
                indizea[1]=c[0]        
                
                self.Con[indizea[0],indizea[1]]=np.random.choice(self.N,1)
                #print("CON add",self.Con)
        else :
            a=np.where(self.Con[1:-self.e]==0)
            ar=np.array(a)
            
            indizea = np.empty(2,dtype=int)
            
            indize=ar[:,np.random.choice(np.size(ar,1),1)]
            b=indize[0]+1
            c=indize[1]
            indizea[0]=b[0]
            indizea[1]=c[0]        
            #print("indizea",indizea)
            self.Con[indizea[0],indizea[1]]=np.random.choice(self.N,1)
            #print("CON add",self.Con)
        return
    
    # MAITE
    def ChangeLink(self):
            """
            Nodo bat lehen konektatu gabe zeon batekin konektaukoa beste konekxio bat hausiz
            """
            if(type(self.K) is int):
                a=np.where(self.Con!=0)
                ar=np.array(a)
                
                indizea = np.empty(2,dtype=int)
                
                indize=ar[:,np.random.choice(np.size(ar,1),1)]
                b=indize[0]
                c=indize[1]
                indizea[0]=b[0]
                indizea[1]=c[0]        
                
                self.Con[indizea[0],indizea[1]]=np.random.choice(self.N,1)
                #print("CON add",self.Con)
            else :
                if (self.e==0)==True:
                    a=np.where(self.Con[1:,:]!=0)
                else:
                    a=np.where(self.Con[1:-self.e,:]!=0)
                
                ar=np.array(a)
                
                indizea = np.empty(2,dtype=int)
                
                indize=ar[:,np.random.choice(np.size(ar,1),1)]
                b=indize[0]+1
                c=indize[1]
                indizea[0]=b[0]
                indizea[1]=c[0]        
                #print("indizea",indizea)
                self.Con[indizea[0],indizea[1]]=np.random.choice(self.N,1)
                #print("CON add",self.Con)
            return
   

    #MAITE
    def DeleteLink(self):
            """
            Konekxio bat kenduko zaio aleatorioki
            
            """
            if(type(self.K) is int):
                a=np.where(self.Con[:,:]!=0)

                ar=np.array(a)
                
                indizea = np.empty(2,dtype=int)
                
                indize=ar[:,np.random.choice(np.size(ar,1),1)]
                b=indize[0]
                c=indize[1]
                indizea[0]=b[0]
                indizea[1]=c[0]        
                
                self.Con[indizea[0],indizea[1]]=0
                #print("CON add",self.Con)
            else :
                if (self.e==0)==True:
                    a=np.where(self.Con[1:,:]!=0)
                else:
                    a=np.where(self.Con[1:-self.e,:]!=0)
                ar=np.array(a)
                indizea = np.empty(2,dtype=int)
                
                indize=ar[:,np.random.choice(np.size(ar,1),1)]
        
                b=indize[0]+1
                c=indize[1]
                indizea[0]=b[0]
                indizea[1]=c[0]        
                #print("indizea",indizea)
                self.Con[indizea[0],indizea[1]]=0
                #print("CON add",self.Con)
            return
    # MAITE  
    def antifragileMultilayer(self, T, b=0, L=0, runs=1, X=None, O=None, fraction=1):
        """
        plot antifragility of BN MULTILAYER 
        """
        
        f=np.zeros(int(self.N*self.Nc/fraction))
        fO=np.zeros(50)
       
        for j in range(runs):
            Cell_initial=[]
           
            for i in range(self.Nc):
                initial = np.random.randint(0, 2, self.N) 
                Cell_initial.append(initial.tolist())
    
            Cell_State=self.RunNetMultilayer(2*T,b, Cell_initial, L)
            
            for i in range(self.Nc):
                State=np.asarray(Cell_State[i])
                State = State[-T:]
                Cell_State[i]=State.tolist()
            
            C0=complexityMultilayer(Cell_State)
            
            if(O!=None):
                for XX in range(1,(self.N*self.Nc)+1):
                    f[XX-1]+=self.funcMultilayer(T, b, L, Cell_initial, XX, O, C0)
                 
            elif(X!=None):
                for OO in range(1,50+1):
                    fO[OO-1]+=self.funcMultilayer(T, b, L, Cell_initial, X, OO, C0)
   
        f=f/runs # average fragility by perturbation
  
        fO=fO/runs

        if(O!=None):
            return f
        elif(X!=None):
            return fO
    #MAITE
    def antifragileMultilayerUnif(self, T, b=0, L=0, runs=1, X=None, O=None, fraction=1):
        """
        plot antifragility of BN MULTILAYER 
        """
        
        f=np.zeros(int(self.N*self.Nc/fraction))
        fO=np.zeros(int(T))
       
        for j in range(runs):
            Cell_initial=[]
           
            for i in range(self.Nc):
                initial = np.random.randint(0, 2, self.N) 
                Cell_initial.append(initial.tolist())
    
            Cell_State=self.RunNetMultilayer(2*T, b, Cell_initial, L)
            
            for i in range(self.Nc):
                State=np.asarray(Cell_State[i])
                State = State[-T:]
                Cell_State[i]=State.tolist()
            
            C0=complexityMultilayer(Cell_State)
            
            if(O!=None):
                for XX in range(1,(self.N*self.Nc)+1):
                    print("X",X)
                    f[XX-1]+=self.funcMultilayerUnif(T, b, L, Cell_initial, XX, O, C0)
                 
            elif(X!=None):
                for OO in range(1,T+1):
                    print("O=",OO)
                    fO[OO-1]+=self.funcMultilayerUnif(T, b, L, Cell_initial, X, OO, C0)
        
        f=f/runs # average fragility by perturbation
        fO=fO/runs

        if(O!=None):
            return f
        elif(X!=None):
            return fO
    #OKARIM
    def antifragile(self, T, runs=1, X=None, O=None, fraction=1):
        """
        plot antifragility of RBN
        """
        
        f=np.zeros(int(self.N/fraction))
        fO=np.zeros(int(T/2))
        pool = multiprocessing.Pool()
        
        for j in range(runs):
            initial = np.random.randint(0, 2, self.N)
            State=self.RunNet(2*T, initial)
            C0 = complexity(State[-T:])
            ########Hemen noa #######
            if ( X and O ) != None:
                fXO= self.func(X,T,initial,O,C0)
                
            elif(O!=None):
                f+=pool.map(partial(self.func, T=T, initial=initial, O=O, C0=C0, fraction=fraction), range(1, int(self.N/fraction)+1))
            elif(X!=None):
                fO+=pool.map(partial(self.func2, T=T, initial=initial, X=X, C0=C0, fraction=fraction), range(1, int(T/2)+1))
        f/=runs # average fragility by perturbation
        fO/=runs
        pool.close()
        if ( X and O ) != None:
            return fXO
        elif(O!=None):
            return f
        elif(X!=None):
            return fO
    #OKARIM
    def func2(self, i, T, initial, X, C0, fraction=1):
        f=np.zeros(int(self.N/fraction))
        State=self.RunNet(T*2, initial, X, i)
        C = complexity(State[-T:])
        f=fragility(C, C0, X, i, self.N, T)
        return f
    #OKARIM
    def func(self, X, T, initial, O, C0, fraction=1):
        f=np.zeros(int(self.N/fraction))
        State=self.RunNet(T*2, initial, X, O)
        C = complexity(State[-T:])
        f=fragility(C, C0, X, O, self.N, T)
        return f
    
    ## MAITE
    def funcMultilayer(self, T, b, L, Cell_initial, X, O, C0, fraction=1):
        Cell_State=self.RunNetMultilayer(T*2, b, Cell_initial, L, X, O)
        for i in range(self.Nc):
            State=np.asarray(Cell_State[i])
            State = State[-T:]
            Cell_State[i]=State.tolist()
            
        C=complexityMultilayer(Cell_State)
        #f=np.zeros(int(self.N/fraction))
            
        f=fragility(C, C0, X, O, self.N*self.Nc, T)
        return f
#        if f < 0:
#            return 1
#        else:
#            return 0
    
    #MAITE    
    def funcMultilayerUnif(self, T, b, L, Cell_initial, X, O, C0, fraction=1):
        Cell_State=self.RunNetMultilayerUnif(T*2, b, Cell_initial, L, X, O)
        for i in range(self.Nc):
            State=np.asarray(Cell_State[i])
            State = State[-T:]
            Cell_State[i]=State.tolist()
            
        C=complexityMultilayer(Cell_State)
        #f=np.zeros(int(self.N/fraction))
            
        f=fragility(C, C0, X, O, self.N*self.Nc, T)
        return f
#        if f < 0:
#            return 1
#        else:
#            return 0


#OKARIM
def complexity(state):
    """
    Measuring Complexity Based on Shanon Entropy 
    state = matrix of a RBN states
    """
    p1=np.sum(state, axis=0)/np.size(state, 0)
    p0=1-p1
    np.place(p0, p0==0.0, 1.0)
    np.place(p1, p1==0.0, 1.0)
    #column by column
    E=-(p0*np.log2(p0)+p1*np.log2(p1)) #Shannon Entropy
    E=np.mean(E)
    C=4*E*(1-E) #Complexity
    return C

#MAITE
def complexityMultilayer(Cell_State):
    """
    Measuring Complexity Based on Shanon Entropy 
    Cell_state = RBN state desberdinai dagokien lista
    """
    Cell_E=np.zeros(len(Cell_State))
    for i in range(len(Cell_State)):
        State=np.asarray(Cell_State[i])     
        p1=np.sum(State, axis=0)/np.size(State, 0)
        p0=1-p1
        np.place(p0, p0==0.0, 1.0)
        np.place(p1, p1==0.0, 1.0)
        E=-(p0*np.log2(p0)+p1*np.log2(p1)) #Shannon Entropy
        E=np.mean(E)
        Cell_E[i]=E
        #column by column

    E=np.mean(Cell_E)   
    C=4*E*(1-E) #Complexity
    return C

#OKARIM
def fragility(C, C0, X, O, N, T):
    """
    C0 = initial complexity
    C = final complexity
    X = how many perturbations
    O = how often the perturbations take place
    N = number of nodes, indexed 0 .. N-1
    T = timesteps
    """
    dx =(X*(T/O))/(N*T) # degree of perturbation
    sigma = np.mean(C-C0) # degree of satisfaction
    #return sigma
    return -sigma*dx

#MAITE
def mutual_info(State):
    #marginal probabilities
    
    H_marginal=np.zeros(np.size(State,1))
    unique_elements, counts_elements = np.unique(State, return_counts=True, axis=0)

    for j in range(np.size(State,1)):
        i=0
        prob0=0
        prob1=0
        for bekt in unique_elements:
            if (bekt[j]==0)==True:
                prob0=prob0 + counts_elements[i]
             
            elif (bekt[j]==1)==True:
                prob1=prob1 + counts_elements[i]    
            i+=1
        prob0/=np.size(State,0)
        prob1/=np.size(State,0)
        if (prob0==0.0)==True:
            prob0=1.0
        if (prob1==0.0)==True:
            prob1=1.0
        H_marginal[j]= -(prob0*log(prob0,2)+prob1*log(prob1,2))
    H_marginal=np.sum(H_marginal) 
    #joint probabilities
    unique_elements1, counts_elements1 = np.unique(State, return_counts=True, axis=1)
    H_joint=np.zeros(len(counts_elements1))
    h=0
    for joint in counts_elements1:
        joint/=np.size(State,1)
        H_joint[h]=-(joint*log(joint,2))
        h+=1
    H_joint=np.sum(H_joint)
    I=H_marginal-H_joint
    
    return I
   

    
# MAITE
class return_values:
    def __init__(self, RNE,RE,NRE,NRNE):
        self.RNE=RNE
        self.RE=RE
        self.NRE=NRE
        self.NRNE=NRNE
   
def sailkapena(attractList0,attractList):

    RNE=0
    RE=0
    NRNE=0
    NRE=0
    
    AtraktoreIgualak=0
    
    for ar0 in attractList0:
        for ar in attractList:
            igualak=0
            for i in ar0:
                for j in ar:
                    if np.array_equiv(i,j)==True :
                        igualak+=1
        
        #hemen hartuet hasierako listako lehenengo atraktore
            #hemen bi atraktore desberdin dazkat
            if (igualak==len(ar0))==True and (igualak==len(ar))==True: 
                AtraktoreIgualak+=1
                
    if (AtraktoreIgualak==len(attractList0))==True and (AtraktoreIgualak==len(attractList))==True: #robust and not evolvable
        RNE+=1          
         
    if (AtraktoreIgualak==len(attractList0))==True and (AtraktoreIgualak < len(attractList))==True: #robust and evolvable
        RE+=1
    
    if (AtraktoreIgualak < len(attractList0))==True and (AtraktoreIgualak==len(attractList))==True:  #not robust and not evolvable
        NRNE+=1
    
    if (AtraktoreIgualak!=0)==True and (AtraktoreIgualak<len(attractList0))==True and (AtraktoreIgualak<len(attractList))==True: #not robust and evolvable
        NRE+=1  
                    
    t=return_values(NRE,RE,RNE,NRNE)
    return t
