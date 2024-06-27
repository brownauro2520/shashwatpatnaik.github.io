import math
import numpy as npd
import sys
from funcs import *
from smoothing import smoothing
main = np.loadtxt('meshes/main.txt')
flap = np.loadtxt('meshes/flap.txt')
slat = np.loadtxt('meshes/slat.txt')
airfoil = np.concatenate((main,flap,slat))

print('Adapting Mesh')

inp = sys.argv[1]
inp2 = str(int(inp)+1)

gri = 'meshes/c'+ inp +'.gri'
gri2 = 'meshes/c'+ inp2 +'.gri'
sol = 'sols/U' + inp + '.txt'

U = np.loadtxt(sol) #loads solution
[E2N,C,B] = readgri(gri) #read gri
mach_plot(U,E2N,C,str(sys.argv[1]))
[I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices

snapnodes = {}; flags = {}

RHh = np.loadtxt('sols/res/RHh.txt')
adj = np.loadtxt('sols/adjoints/adj.txt')
epsilon = RHh*adj

method = 'res'
comp = 'comb'

flags = flagerror(flags,epsilon,Area,method,comp,int(inp))
print(len(flags))
[E2N,B,C,U,snapnodes] = localref(E2N,C,B,U,flags,snapnodes,0) #localref

snap(C,B,airfoil);  #snapping
[I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices
smoothing(E2N,C,bcords,3,0.7)

makegri(E2N,C,B,gri2) #makes new gri of refined mesh
np.savetxt('sols/U'+inp2+'.txt',U)

print('Mesh Verification: ' + str(meshver(E2N,I2E,B2E,In,Bn,Il,Bl))) #mesh verification