# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 16:47:39 2023

@author: 15413
"""

import numpy as np
import matplotlib.pyplot as plt

def readgri(fname):
    f = open(fname, 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    V = list([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    B = []; Bname = []
    for i in range(NB):
        s = f.readline().split(); Nb = int(s[0]); Bname.append(s[2])
        Bi = list([[int(s) for s in f.readline().split()] for n in range(Nb)])
        B.append(Bi)
    # read elements
    Ne0 = 0; E = []
    while (Ne0 < Ne):
        s = f.readline().split(); ne = int(s[0])
        Ei = list([[int(s) for s in f.readline().split()] for n in range(ne)])
        E = Ei if (Ne0==0) else np.concatenate((E,Ei), axis=0)
        Ne0 += ne
    f.close()
    Mesh = {'V':V, 'E':E, 'B':B, 'Bname':Bname }
    return Mesh

def split(L, n):
    lenth = int(len(L)/n)
    k, m = divmod(len(L), lenth)
    return list(L[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(lenth))

def importU(fname,p):
    if p == 0:
        N = 1
    elif p == 1:
        N = 3
    elif p == 2:
        N = 6
    elif p == 3:
        N = 10
        
    state = []
    for line in open(fname, 'r'):
        A = float(line)
        state.append(A)    
    U = split(state,4)
    U = split(U,N)
    return U

def mach(U):
    gamma = 1.4
    r = U[0]
    q = np.sqrt(U[1]**2 + U[2]**2)/r
    p = (gamma-1)*(U[3] - 0.5*r*q**2)
    c = np.sqrt(gamma*p/r)
    M = q/c
    return M 

def flinspace(number):
    x = []
    y = [] 
    for i in range(1,number+1):
        x.append( (i-1) * (1/(number-1)))
        y.append( (i-1) * (1/(number-1)))
        
    V = []
    for j in range(number):
        for i in range(number-j):
            V.append([x[i],y[j]])
    Vp = list(zip(*V))
    return Vp


def TriLagrange(x,y,order): #create phi matrix component for singular 1d quad point
    if order == 0:
        phi = np.zeros((1,1))
        phi[0] = 1.0;
    elif order == 1:
        phi = np.zeros((3,1))
        phi[0] = 1-x-y;
        phi[1] =   x  ;
        phi[2] =     y;
    elif order == 2:
        phi = np.zeros((6,1))
        phi[0] = 1.0-3.0*x-3.0*y+2.0*x*x+4.0*x*y+2.0*y*y
        phi[1] = 4.0*x-4.0*x*x-4.0*x*y
        phi[2] = -x+2.0*x*x
        phi[3] = 4.0*y-4.0*x*y-4.0*y*y
        phi[4] = 4.0*x*y
        phi[5] = -y+2.0*y*y
    elif order == 3:
        phi = np.zeros((10,1))
        phi[0] = 1.0-11.0/2.0*x-11.0/2.0*y+9.0*x*x+18.0*x*y+9.0*y*y-9.0/2.0*x*x*x-27.0/2.0*x*x*y-27.0/2.0*x*y*y-9.0/2.0*y*y*y;
        phi[1] = x-9.0/2.0*x*x+9.0/2.0*x*x*x;
        phi[2] = y-9.0/2.0*y*y+9.0/2.0*y*y*y;
        phi[3] = -9.0/2.0*x*y+27.0/2.0*x*x*y;
        phi[4] = -9.0/2.0*x*y+27.0/2.0*x*y*y;
        phi[5] = -9.0/2.0*y+9.0/2.0*x*y+18.0*y*y-27.0/2.0*x*y*y-27.0/2.0*y*y*y;
        phi[6] = 9.0*y-45.0/2.0*x*y-45.0/2.0*y*y+27.0/2.0*x*x*y+27.0*x*y*y+27.0/2.0*y*y*y;
        phi[7] = 9.0*x-45.0/2.0*x*x-45.0/2.0*x*y+27.0/2.0*x*x*x+27.0*x*x*y+27.0/2.0*x*y*y;
        phi[8] = -9.0/2.0*x+18.0*x*x+9.0/2.0*x*y-27.0/2.0*x*x*x-27.0/2.0*x*x*y;
        phi[9] = 27.0*x*y-27.0*x*x*y-27.0*x*y*y;
    else:
        print('basis order not supported')
    return phi





def plot_mach(Mesh,state_fname,p):
    Nodec = Mesh['V']; E2N = Mesh['E']
    
    if p == 0:
        X = []
        Y = []
        for l in range(len(Nodec)):
            X.append(Nodec[l][0])
            Y.append(Nodec[l][1])
        U = importU(state_fname,p)
        elem_mach = []
        for i in range(len(U)):
            elem_mach.append(mach(U[i][0]))
            
        E = np.zeros([len(E2N),3])
        for a1 in range(len(E2N)):
            for a2 in range(3):
                E[a1][a2] = int(E2N[a1][a2])-1
        plt.figure()
        plt.tripcolor(X,Y,triangles=E,facecolors=elem_mach,shading='flat')
        plt.title('Mach Number visualization for Mesh0 p = 0')
        plt.set_cmap('jet')
        cbar=plt.colorbar(orientation='horizontal', pad=0.08, fraction=.08)
        cbar.ax.tick_params(labelsize=8)
        plt.savefig('Mach Number visualization.png', dpi=1000)
        return p
    
    elif p == 1:
        N = 3
    elif p == 2:
        N = 5
    elif p == 3:
        N = 7

    Vp = flinspace(N)
    phi = []
    for i in range(len(Vp[0])):
        phi.append(TriLagrange(Vp[0][i],Vp[1][i],p))

    U = importU(state_fname,p)
    elem_mach = U.copy()       
    for i in range(len(U)):
        for j in range(len(U[0])):
            elem_mach[i][j] = mach(U[i][j])
    Up = []       
    for i in range(len(U)):
        Up.append(np.dot(elem_mach[i],phi))   
        
    M = np.zeros((N,N))
    k = 0
    for j in range(0,N):
        for i in range(0,N-j):
            k = k + 1
            M[i][j] = k

    E = np.zeros(((N-1)*(N-1),3))
    k = 0
    for j in range(0,N-1):
        for i in range(0,N-j-1):
            E[k][:] = [M[i][j],M[i+1][j],M[i][j+1]]
            k = k + 1
            if i < (N-j-2):
                E[k][:] = [M[i+1][j],M[i+1][j+1],M[i][j+1]]
                k = k + 1

    total_nodes = []
    total_E2N = []
    k = 0
    for i in range(len(E2N)):
        p1= Nodec[E2N[i][0]-1]
        p2= Nodec[E2N[i][1]-1]
        p3= Nodec[E2N[i][2]-1]
        elem_nodes = []
        for j in range(len(Vp[0])):
            new_P = [0,0]
            new_P[0] = p1[0]*(1-Vp[0][j]-Vp[1][j]) + p2[0]*Vp[0][j] + p3[0]*Vp[1][j]
            new_P[1] = p1[1]*(1-Vp[0][j]-Vp[1][j]) + p2[1]*Vp[0][j] + p3[1]*Vp[1][j]
            elem_nodes.append(new_P)
            
        M = np.zeros((N,N))
        for j in range(0,N):
            for i in range(0,N-j):
                k = k + 1
                M[i][j] = k
                
        
        E = np.zeros(((N-1)*(N-1),3))
        f = 0
        for j in range(0,N-1):
            for i in range(0,N-j-1):
                E[f][:] = [M[i][j],M[i+1][j],M[i][j+1]]
                f = f + 1
                if i < (N-j-2):
                    E[f][:] = [M[i+1][j],M[i+1][j+1],M[i][j+1]]
                    f = f + 1
        total_nodes += elem_nodes
        total_E2N += list(E)


    total_mach = []
    for i in range(len(Up)):
        total_mach += list(Up[i])

    E = np.zeros([len(total_E2N),3])
    for i in range(len(total_E2N)):
        for j in range(3):
            E[i][j] = int(total_E2N[i][j])-1

    X = []
    Y = []
    for i in range(len(total_nodes)):
        X.append(total_nodes[i][0])
        Y.append(total_nodes[i][1])
    with open('X.txt', "w") as f: 
        for i in range(len(X)):
            a = X[i]
            f.write("{}\n".format(a))
    with open('Y.txt', "w") as f: 
        for i in range(len(Y)):
            a = Y[i]
            f.write("{}\n".format(a))
    with open('mach.txt', "w") as f: 
        for i in range(len(total_mach)):
            a = float(total_mach[i])
            f.write("{}\n".format(a))
    with open('E2N.txt', "w") as f: 
        for i in range(len(total_E2N)):
            a = int(total_E2N[i][0])
            b = int(total_E2N[i][1])
            c = int(total_E2N[i][2])
            f.write("{} {} {}\n".format(a,b,c))
             
    
    
    sube_mach = []
    for i in range(len(E)):
        M = ( float(total_mach[int(E[i][0])])
             +float(total_mach[int(E[i][1])])
             +float(total_mach[int(E[i][2])]) )/3
        sube_mach.append(M)
    plt.figure()
    plt.tripcolor(X,Y,triangles=E,facecolors=sube_mach,shading='flat')
    #plt.axis([-0.05, 0.05, -0.03,0.03])
    plt.title('Mach Number visualization for Mesh0 at p = 3')
    plt.set_cmap('jet')
    cbar=plt.colorbar(orientation='horizontal', pad=0.08, fraction=.08)
    cbar.ax.tick_params(labelsize=8)
    plt.savefig('Mach Number visualization.png', dpi=1000)
    plt.show()

    return Up



Mesh = readgri('b1.gri')

Up = plot_mach(Mesh,'U1P3.txt',3)


























































