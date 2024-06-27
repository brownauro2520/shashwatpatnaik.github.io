from funcs import *
import sys
import numpy as np

mesh = str(int(sys.argv[1]))
[E2N,C,B] = readgri('meshes/c'+mesh+'.gri')
[I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices


if sys.argv[2] == "mach":
    sol = 'sols/U'+mesh+'.txt' #change this
    E = np.loadtxt(sol)
else:
    res = '../cpp/FV/sols/res/R'+mesh+'.txt' #change this
    R = np.loadtxt(res)
    adj = '../cpp/FV/sols/adjoints/adjlift'+mesh+'.txt' #change this
    adj = np.loadtxt(adj)
    Area = np.ones((len(E2N))) #change this for residual based view
    #adj = np.ones((len(R)))
    E = R*adj

plotmach(E,E2N,C,sys.argv[2],Area) ##coarse space

# mesh = str(int(sys.argv[1])+1)
# [E2N,C,B] = readgri('../meshes/c'+sys.argv[1]+'.gri')
# [I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices

# res = '../cpp/FV/sols/res/R'+mesh+'.txt' #change this
# R = np.loadtxt(res)
# #adj = '../cpp/FV/sols/adjoints/adj'+mesh+'.txt' #change this
# #adj = np.loadtxt(adj)
# #Area = np.ones((len(E2N))) #change this for residual based view
# adj = np.ones((len(R)))
# E = abs(R*adj)
# EE = np.zeros((int(len(E)/4)))
# c = len(EE)
# for i in range(int(len(EE)/4)): #loop over coarse space elements
#     EE[4*i:4*(i+1)] = E[4*i:4*(i+1)] + E[c:c+4] + E[c+4:c+8] + E[c+8:c+12]

# plotmach(EE,E2N,C,sys.argv[2],Area) ##fine space