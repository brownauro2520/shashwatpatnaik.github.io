import math
import numpy as np
import matplotlib.pyplot as plt
from funcs import *
from smoothing import smoothing, calcneigh
import sys

nr = str(sys.argv[1])
space = sys.argv[2]

[E2N,C,B] = readgri('meshes/c'+nr+'.gri') #read gri
[I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices

########## save matrices for input in c++
np.savetxt('matrices/'+space+'/E2N.txt',E2N)
np.savetxt('matrices/'+space+'/I2E.txt',I2E)
np.savetxt('matrices/'+space+'/B2E.txt',B2E)
np.savetxt('matrices/'+space+'/In.txt',In)
np.savetxt('matrices/'+space+'/Bn.txt',Bn)
np.savetxt('matrices/'+space+'/Il.txt',Il)
np.savetxt('matrices/'+space+'/Bl.txt',Bl)
np.savetxt('matrices/'+space+'/A.txt',Area)
np.savetxt('matrices/'+space+'/C.txt',C)