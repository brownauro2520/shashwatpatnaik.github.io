import math
import numpy as npd
import sys
from funcs import *
from smoothing import smoothing
main = np.loadtxt('meshes/main.txt')
flap = np.loadtxt('meshes/flap.txt')
slat = np.loadtxt('meshes/slat.txt')
airfoil = np.concatenate((main,flap,slat))

inp = sys.argv[1]
inp2 = str(int(inp)+1)

gri = 'meshes/c'+ inp +'.gri'
gri2 = 'meshes/c.gri'
sol = 'sols/U' + inp + '.txt'

U = np.loadtxt(sol) #loads solution
[E2N,C,B] = readgri(gri) #read gri
snapnodes = {}
[E2N,B,C,U,snapnodes] = localref(E2N,C,B,U,'all',snapnodes,0) #localref

snap(C,B,airfoil);  #snapping
[I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices
smoothing(E2N,C,bcords,0,0.9)

makegri(E2N,C,B,gri2) #makes new gri of refined mesh

[I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices
print('Mesh Verification: ' + str(meshver(E2N,I2E,B2E,In,Bn,Il,Bl))) #mesh verification

########## save matrices for input in c++
np.savetxt('sols/inj/Uinj.txt',U)
np.savetxt('matrices/fine/E2N.txt',E2N)
np.savetxt('matrices/fine/I2E.txt',I2E)
np.savetxt('matrices/fine/B2E.txt',B2E)
np.savetxt('matrices/fine/In.txt',In)
np.savetxt('matrices/fine/Bn.txt',Bn)
np.savetxt('matrices/fine/Il.txt',Il)
np.savetxt('matrices/fine/Bl.txt',Bl)
np.savetxt('matrices/fine/A.txt',Area)
np.savetxt('matrices/fine/C.txt',C)