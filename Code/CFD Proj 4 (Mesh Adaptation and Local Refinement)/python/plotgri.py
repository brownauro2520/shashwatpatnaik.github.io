import numpy as np
import matplotlib.pyplot as plt
import sys

#-----------------------------------------------------------
def readgri(fname):
    f = open(fname, 'r')
    Nn, Ne, dim = [int(s) for s in f.readline().split()]
    # read vertices
    V = np.array([[float(s) for s in f.readline().split()] for n in range(Nn)])
    # read boundaries
    NB = int(f.readline())
    B = []; Bname = []
    for i in range(NB):
        s = f.readline().split(); Nb = int(s[0]); Bname.append(s[2])
        Bi = np.array([[int(s)-1 for s in f.readline().split()] for n in range(Nb)])
        B.append(Bi)
    # read elements
    Ne0 = 0; E = []
    while (Ne0 < Ne):
        s = f.readline().split(); ne = int(s[0])
        Ei = np.array([[int(s)-1 for s in f.readline().split()] for n in range(ne)])
        E = Ei if (Ne0==0) else np.concatenate((E,Ei), axis=0)
        Ne0 += ne
    f.close()
    Mesh = {'V':V, 'E':E, 'B':B, 'Bname':Bname }
    return Mesh

#-----------------------------------------------------------
def plotmesh(Mesh, mesh):
    V = Mesh['V']; E = Mesh['E']; 
    f = plt.figure(figsize=(12,12))
    #plt.tripcolor(V[:,0], V[:,1], triangles=E)
    plt.triplot(V[:,0], V[:,1], E, 'k-')
    plt.axis('off')#plt.tick_params(axis='both', labelsize=12)
    f.tight_layout()
    plt.axis('equal')
    plt.xlim(-0.5,1.5)
    plt.ylim(-0.2,0.1)
    plt.xlabel('X - Coordinate'); plt.ylabel('Y - Coordinate'); plt.title('Mesh Plot')
    plt.savefig('pictures/mesh' + mesh + '.pdf')
    plt.close(f)

#-----------------------------------------------------------
def main():
    Mesh = readgri("meshes/c"+sys.argv[1]+".gri")
    plotmesh(Mesh, sys.argv[1]);


if __name__ == "__main__":
    main()

