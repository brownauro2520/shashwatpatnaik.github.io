import math
import numpy as np
import matplotlib.pyplot as plt
from spline import snap2splines
from dgfuncs import *

def readgri(fname):
    f = open(fname, 'r') #open file
    lbl = [] #placeholder for all lines in .gri
    for line in f:
        lbl.append(line.strip()) #store all lines

    nn = int(lbl[0].split()[0]) #number of nodes
    ne = int(lbl[0].split()[1]) #number of elements

    C = [] #list to store coordinates
    E2N = [] #list to store element-to-nodes
    B = {} #dictionary to store boundary faces

    for i in range(1, nn+1): #loop over the beginning of .gri
        c1 = float(lbl[i].split()[0])
        c2 = float(lbl[i].split()[1])
        C.append([c1,c2]) #store all coordinates

    for i in range(nn + 2, len(lbl) - ne - 1): #loop over boundary groups and create map
        if len(lbl[i].split()) == 3: #we discovered new boundary group
            name = lbl[i].split()[2] #save name
        else:
            n1 = lbl[i].split()[0]
            n2 = lbl[i].split()[1]
            B[n1 + ' ' + n2] = int(name) #create dictionary that stores the node sequence according to bgroup

    for i in range(len(lbl) - ne, len(lbl)): #loop over last part of .gri
        e1 = int(lbl[i].split()[0])
        e2 = int(lbl[i].split()[1])
        e3 = int(lbl[i].split()[2])
        E2N.append([e1, e2, e3]) #store all elements
    return E2N, C, B

def makegri(E2N,C,B,fname):
    f = open(fname,'w')
    f.write(str(len(C)) + ' ' + str(len(E2N)) + ' 2\n') #write the first line, nr nodes, nr elems and dim
    for c in C:
        f.write(str(c[0]) + ' ' + str(c[1]) + '\n') #write all the coordinates
    lenc = {} #dictionary to count boundary groups
    for b in B:
        if B[b] in lenc:
            lenc[B[b]] += ' ' + b #append face to existing bgroup
        else:
            lenc[B[b]] = b #create bgroup for face

    f.write(str(len(lenc)) + '\n') #print len of bgroup

    for l in lenc: #for each face in bgroup
        v = lenc[l].split() #create list containing boundary groups
        f.write(str(int(len(v)/2)) + ' 2 ' + str(l) + '\n') #write the header for each bgroup
        for i in range(len(v)//2): #for each pair of coords
            f.write(str(v[2*i]) + ' ' + str(v[2*i + 1]) + '\n') #print face

    f.write(str(len(E2N)) + ' 2 TriLagrange\n') #write header for E2N matrices
    for e in E2N: #for each element
        f.write(str(e[0]) + ' ' + str(e[1]) + ' ' + str(e[2]) + '\n') #print out each node

def makematrices(E2N,B,C):
    F = {} #dictionary to store all possible faces
    N = []
    L = []
    for i, el in enumerate(E2N): #build a map to store all faces
        for j in range(3):
            key = str(el[j]) + ' ' + str(el[(j+1)%3])
            F[key] = [i + 1, j + 1] #stores the element number and local face number

    I2E, B2E, In, Bn, Il, Bl, Area, N, L, bcords = [], [], [], [], [], [], [], [], [], {}
    for i, el in enumerate(E2N): #loop through all elements, i+1 is element number, j+1 is face number
        l = []
        for j in range(3):
            key = str(el[j]) + ' ' + str(el[(j+1)%3]) #eg. '6 10'
            key_rev = str(el[(j+1)%3]) + ' ' + str(el[j]) #eg. '10 6'

            dx = C[el[(j+1)%3]-1][0]-C[el[j]-1][0] #change in x
            dy = C[el[(j+1)%3]-1][1]-C[el[j]-1][1] #change in y
            l.append(math.sqrt(dx**2+dy**2)+1E-16) #length of face
            L.append(l[-1])
            N.append([dy/l[-1],-dx/l[-1]])
            if key in B: #its a boundary face
                nodes = key.split()
                bcords[int(nodes[0])] = 0; bcords[int(nodes[1])] = 0
                B2E.append([i + 1, j + 1, B[key]]) #[element, local face nr., bgroup]
                Bn.append([dy/l[-1],-dx/l[-1]]) #normalized boundary face vector
                Bl.append(l[-1])
            else: #it's an interior face
                elem1 = F[key] #save the element and local face number
                elem2 = F[key_rev] #same thing but for the reverse face combiation
                if elem1[0] < elem2[0]: #check which has lowest el. nr and also avoid duplicates
                    I2E.append(elem1 + elem2)
                    In.append([dy/l[-1],-dx/l[-1]]) #normalized interior face vector
                    Il.append(l[-1])

        p = sum(l[-3:])/2 #perimeter of element
        Area.append(math.sqrt(p * (p - l[-1]) * (p - l[-2]) * (p - l[-3]))) #Heron's formula for area of triangle
    return I2E, B2E, In, Bn, Il, Bl, Area, L, N, bcords

def localref(E2N,C,B,U,flags,snapnodes,order): #order
    #----------------
    order = 0
    npp = ((order+1)*(order+2))/2
    split = {} #dictionary to keep track of which edges are already split
    #------------------------------------------------------------------------------------
    #all edges flagged
    nn = len(E2N) #number of elements
    nc = len(C) #number of nodes
    if flags == 'all':
        flags = {i: 0 for i in range(1, nn+1)}
    for f in flags: #for each flagged element
        newnodes = [] #place for new/existing nodes
        for i in range(3): #for each node on flagged element
            key = str(E2N[f-1][i]) + ' ' + str(E2N[f-1][(i+1)%3]) #eg. '6 10'
            keyrev = str(E2N[f-1][(i+1)%3]) + ' ' + str(E2N[f-1][i]) #eg. '6 10'
            if not key in split: #edge split for the first time
                x1 = C[E2N[f-1][i]-1][0] #x coordinate of first node
                y1 = C[E2N[f-1][i]-1][1] #y coordinate of first node

                x2 = C[E2N[f-1][(i+1)%3]-1][0] #x coordinate of second node
                y2 = C[E2N[f-1][(i+1)%3]-1][1] #y coordinate of second node
                C.append([(x1+x2)/2,(y1+y2)/2]) #create new coordinates
                nc += 1 #one mode node is created
                split[key] = nc #this edge has been split and created node nc
                split[keyrev] = nc #save the reverse combination too

                ##re-assign boundary groups for boundary edges
                if key in B: #check if edge is boundary
                    if B[key] != -1:
                        snapnodes[nc] = [(x1+x2)/2,(y1+y2)/2]
                    n1 = key.split()[0] #node1
                    n2 = key.split()[1] #node2
                    B[n1 + ' ' + str(nc)] = B[key] #part 1 after split
                    B[str(nc) + ' ' + n2] = B[key] #part 2 after split
                    del B[key] #delete original boundary edge
            newnodes.append(split[key])
        E2N.append([newnodes[0],E2N[f-1][1],newnodes[1]]) #create new elements
        E2N.append([newnodes[1],E2N[f-1][2],newnodes[2]]) #create new elements
        E2N.append([newnodes[2],E2N[f-1][0],newnodes[0]]) #create new elements
        E2N[f-1] = newnodes #replace element to retain ordering
        for i in range(3): #transfer the state to new nodes
            #U.append(U[f-1])
            U = np.concatenate((U,U[int((f-1)*4*npp):int((f)*4*npp)]))
    #------------------------------------------------------------------------------------
    #single/double side refinement
    splitedges = {} #what faces are also split
    for i, el in enumerate(E2N): #for each element
        if not i + 1 in flags: #an element that's not initially flagged
            splitedges[i+1] = []
            for j in range(3): #for each node/edge
                key = str(el[j]) + ' ' + str(el[(j+1)%3])
                if key in split:
                    splitedges[i+1].append(j+1)

    for el in splitedges:
        e = E2N[el-1] #nodes in that element to be split
        face = splitedges[el] #local face number
        if not len(face) == 0:
            key1 = str(e[face[0]-1]) + ' ' + str(e[face[0]%3]) #key for face 1
            nc1 = split[key1] # new/existing node that's on the split edge
        if len(face) == 1: #single edge refine
            E2N.append([e[(face[0]+1)%3],e[face[0]-1],nc1])
            E2N[el-1] = [e[(face[0]+1)%3],nc1,e[face[0]%3]]
            U = np.concatenate((U,U[int((el-1)*4*npp):int((el)*4*npp)]))
            #U.append(U[el-1])
        
        if len(face) == 2: #double edge refine
            key2 = str(e[face[1]-1]) + ' ' + str(e[face[1]%3]) #key for face 2
            nn = [split[key1],split[key2]] #vec of nodes to split
            l = 0 #list to store lenghts of sides
            for j in range(3): #for each node
                dx = C[e[(j+1)%3]-1][0]-C[e[j]-1][0] #change in x
                dy = C[e[(j+1)%3]-1][1]-C[e[j]-1][1] #change in y
                if math.sqrt(dx**2+dy**2) > l and j + 1 in face: #check that largest angle isnt illegal
                    l = math.sqrt(dx**2+dy**2)
                    f = j + 1 #localface with node across largest angle
                    key = str(e[j]) + ' ' + str(e[(j+1)%3]) #key of edge with node across largest angle
            e1 = [e[(f+1)%3],e[f-1],split[key]]
            e2 = [e[(f+1)%3],split[key],e[f%3]]

            nn.remove(split[key])
            #loop through all edges of e1 e2, find out if one split[key] for that face equals nn[0]. if so e becomes that element
            for i in [e1,e2]:
                for j in range(3):
                    key = str(i[j]) + ' ' + str(i[(j+1)%3])
                    if key in split:
                        e = i
                        if e == e2:
                            E2N.append(e1)
                        else:
                            E2N.append(e2)


            for j in range(3):
                key = str(e[j]) + ' ' + str(e[(j+1)%3]) #key of edge
                if key in split: #if key has been split, remaining edge is identified
                    f = j + 1 #local face number of new element
                    E2N.append([e[(f+1)%3],split[key],e[f%3]])
                    E2N[el-1] = [e[(f+1)%3],e[f-1],split[key]]
                    #U.append(U[el-1])
                    for i in range(2):
                        U = np.concatenate((U,U[int((el-1)*4*npp):int((el)*4*npp)]))
        if len(face) == 3: #triple refine
            f = el
            newnodes = [] #place for new/existing nodes
            for i in range(3): #for each node on flagged element
                key = str(E2N[f-1][i]) + ' ' + str(E2N[f-1][(i+1)%3]) #eg. '6 10'
                keyrev = str(E2N[f-1][(i+1)%3]) + ' ' + str(E2N[f-1][i]) #eg. '6 10'
                if not key in split: #edge split for the first time
                    x1 = C[E2N[f-1][i]-1][0] #x coordinate of first node
                    y1 = C[E2N[f-1][i]-1][1] #y coordinate of first node

                    x2 = C[E2N[f-1][(i+1)%3]-1][0] #x coordinate of second node
                    y2 = C[E2N[f-1][(i+1)%3]-1][1] #y coordinate of second node
                    C.append([(x1+x2)/2,(y1+y2)/2]) #create new coordinates
                    nc += 1 #one mode node is created
                    split[key] = nc #this edge has been split and created node nc
                    split[keyrev] = nc #save the reverse combination too

                    ##re-assign boundary groups for boundary edges
                    if key in B: #check if edge is boundary
                        if B[key] != -1:
                            snapnodes[nc] = [(x1+x2)/2,(y1+y2)/2]
                        n1 = key.split()[0] #node1
                        n2 = key.split()[1] #node2
                        B[n1 + ' ' + str(nc)] = B[key] #part 1 after split
                        B[str(nc) + ' ' + n2] = B[key] #part 2 after split
                        del B[key] #delete original boundary edge
                newnodes.append(split[key])
            E2N.append([newnodes[0],E2N[f-1][1],newnodes[1]]) #create new elements
            E2N.append([newnodes[1],E2N[f-1][2],newnodes[2]]) #create new elements
            E2N.append([newnodes[2],E2N[f-1][0],newnodes[0]]) #create new elements
            E2N[f-1] = newnodes #replace element to retain ordering
            for i in range(3): #transfer the state to new nodes
                #U.append(U[f-1])
                U = np.concatenate((U,U[int((f-1)*4*npp):int((f)*4*npp)]))
    return E2N, B, C, U, snapnodes

def closest_coordinate(point, file_path): #finds the closest point
    closest_point = None
    closest_distance = float("inf")
    with open(file_path, "r") as f:
        for line in f:
            x, y = map(float, line.strip().split())
            distance = math.sqrt((x - point[0])**2 + (y - point[1])**2)
            if distance < closest_distance:
                closest_distance = distance
                closest_point = (x, y)
    return closest_point

def localr(E2N,cord,x1,y1,r):
    flags = {}
    #finding the centroid of the trangle
    centeroid = []
    for i in range(len(E2N)):
        c = []
        for side in range(3):
            node = E2N[i][side]
            c.append(cord[node-1])
        centroid_x = (c[0][0] + c[1][0] + c[2][0]) / 3
        centroid_y = (c[0][1] + c[1][1] + c[2][1]) / 3  
        centeroid.append([i, centroid_x, centroid_y])
    #flagging the centroids
    for i in range(len(E2N)):
        elm, x2, y2 = centeroid[i]
        #distance formula
        r_prime = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
        if r_prime <= r:
            flags[elm+1] = 0
    return flags

def snap(C,B,airfoil):
    for bb in B:
        if B[bb] != -1:
            for s in bb.split():
                C[int(s)-1] = snap2splines(C[int(s)-1],airfoil)
    snapnodes = {}

def meshver(E2N,I2E,B2E,In,Bn,Il,Bl): #mesh verification
    N = In + Bn
    L = Il + Bl
    IB =  I2E + B2E
    errorE = [[0,0]]*len(E2N)
    for i in range(len(N)):
        e1 = IB[i][0]-1
        e2 = IB[i][2]-1
        if IB[i][2] < 0: #boundary face
            errorE[e1][0] += N[i][0]*L[i]
            errorE[e1][1] += N[i][1]*L[i]
        else: #interior face
            errorE[e1][0] += N[i][0]*L[i]
            errorE[e1][1] += N[i][1]*L[i]
            errorE[e2][0] -= N[i][0]*L[i]
            errorE[e2][1] -= N[i][1]*L[i]

    maxval = abs(max(errorE[0]))
    for i in range(1,len(E2N)):
        for j in range(2):
            if abs(errorE[i][j]) > maxval:
                maxval = abs(errorE[i][j])
    return maxval

def flagedges(solution,gri,r):
    U = np.loadtxt(solution)
    [E2N,C,B] = readgri(gri) #read gri
    [I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #create I2E, B2E, In, Bn, Area matrices

    M = []; p = []
    for i in range(len(U)): #iterate over all cells
        gamma = 1.4
        v = np.array([U[i,1]/U[i,0], U[i,2]/U[i,0]]) #interior velocity
        p.append((gamma-1)*(U[i,3] - 0.5*U[i,0]*(v[0]**2 + v[1]**2))) #pressure
        c = np.sqrt(p[i]*gamma/U[i,0]) #speed of sound
        M.append(np.sqrt(v[0]**2+v[1]**2)/c) #Mach number

    IB = I2E + B2E
    L = Il + Bl
    N = In + Bn
    edgeerror = []
    for i, edge in enumerate(IB):
        if len(edge) == 3: #then we are in B2E
            if edge[2] == -1: #farfield
                e = 0
            else: #touching some airfoil element
                n = np.array(N[i]) #normal at boundary
                v_plus = [U[edge[0],1]/U[edge[0],0], U[edge[0],2]/U[edge[0],0]] #interior velocity
                v_b = v_plus - np.dot(v_plus,n)*n #tangential velocity at the boundary
                rhoE = U[edge[0],3] #total energy 
                pb = (gamma-1)*(rhoE - 0.5*U[edge[0],0]*(v_b[0]**2+v_b[1]**2)) #boundary pressure
                c = np.sqrt(gamma*pb/U[edge[0],0]) #speed of sound
                MT = math.sqrt((v_plus[0]*N[i][0])**2 + (v_plus[1]*N[i][1])**2)/c
                e = abs(MT)*L[i]**r
        else: #interior faces
            e = abs(M[edge[0]-1]-M[edge[2]-1])*L[i]**r #error for interior faces
        edgeerror.append(e)


    slist = sorted(enumerate(edgeerror), key=lambda x: x[1], reverse = True)
    slist = slist[0:round(len(slist)*0.03)]
    indexesint = {f: 0 for f in [i[0] for i in slist]} #flagged edges
    flags = {} #flag edges
    for i, edge in enumerate(IB):
        if i in indexesint:
            if edge[2] < 0:
                flags[edge[0]] = 0
            else:
                flags[edge[0]] = 0
                flags[edge[2]] = 0
    return flags

def nm(N, U):
    gamma = 1.4 #normal at boundary
    v_plus = [U[1]/U[0], U[2]/U[0]] #interior velocity
    v_b = v_plus - np.dot(v_plus,N)*N #tangential velocity at the boundary
    rhoE = U[3] #total energy 
    pb = (gamma-1)*(rhoE - 0.5*U[0]*(v_b[0]**2+v_b[1]**2)) #boundary pressure, scalar
    c = np.sqrt(gamma*pb/U[0]) #speed of sound
    MT = math.sqrt((v_plus[0]*N[0])**2 + (v_plus[1]*N[1])**2)/c
    return MT

def m(U):
    gamma = 1.4
    v = np.array([U[1]/U[0],U[2]/U[0]])
    p = (gamma-1)*(U[3] - 0.5*U[0]*(v[0]**2 + v[1]**2))
    c = np.sqrt(p*gamma/U[0])
    M = np.sqrt(v[0]**2 + v[1]**2)/c
    return M

def curveflag(U,C2E,E2NC,CC,gri,r):
    [E2N,C,B] = readgri(gri) #read gri
    [I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #create I2E, B2E, In, Bn, Area matrices

    IB = I2E + B2E
    edgeerror = []
    for i, edge in enumerate(IB):
        xq, wq = quad1d()
        el = E2NC[int(C2E[edge[0]-1]-1)] #[1 2 3 4 5 6]
        lf = edge[1] #localface
        norms, dsig = lagnorm(el,lf,CC) #get normals at quad points
        
        if len(edge) == 4: #interior edge

            UelL = U[24*(edge[0]-1):24*edge[0]]
            UelR = U[24*(edge[2]-1):24*edge[2]]

            U1L = []; U2L = []; U3L = []; U4L = []
            U1R = []; U2R = []; U3R = []; U4R = []
            for i in range(6):
                U1L.append(UelL[4*i]) #rho
                U2L.append(UelL[4*i + 1]) #rho U
                U3L.append(UelL[4*i + 2]) #rho V
                U4L.append(UelL[4*i + 3]) #rho E
                U1R.append(UelR[4*i]) #rho
                U2R.append(UelR[4*i + 1]) #rho U
                U3R.append(UelR[4*i + 2]) #rho V
                U4R.append(UelR[4*i + 3]) #rho E

            U1QL = quadstate(U1L,edge[1],2,6)
            U2QL = quadstate(U2L,edge[1],2,6)
            U3QL = quadstate(U3L,edge[1],2,6)
            U4QL = quadstate(U4L,edge[1],2,6)
            
            U1QR = quadstate(U1L,edge[3],2,6)
            U2QR = quadstate(U2L,edge[3],2,6)
            U3QR = quadstate(U3L,edge[3],2,6)
            U4QR = quadstate(U4L,edge[3],2,6)

            gamma = 1.4; mLQ = []; mRQ = []
            for i in range(7):
                vL = np.array([U2QL[i]/U1QL[i], U3QL[i]/U1QL[i]])
                vR = np.array([U2QR[i]/U1QR[i], U3QR[i]/U1QR[i]])
                pL = (gamma-1)*(U4QL[i] - 0.5*U1QL[i]*(vL[0]**2 + vL[1]**2))
                pR = (gamma-1)*(U4QR[i] - 0.5*U1QR[i]*(vR[0]**2 + vR[1]**2))
                cL = np.sqrt(pL*gamma/U1QL[i])
                cR = np.sqrt(pR*gamma/U1QR[i])
                mLQ.append(np.sqrt(vL[0]**2 + vL[1]**2)/cL)
                mRQ.append(np.sqrt(vR[0]**2 + vR[1]**2)/cR)

            e = 0
            for i in range(7):
                e += wq[i]*abs(mLQ[i]-mRQ[6-i])*dsig[i] #include norm of stan
            edgeerror.append(L[i]**r*e)
            #print(edgeerror[-1])
            
        elif edge[2] != -1: #boundary edge
            

            Uel = U[24*(edge[0]-1):24*edge[0]] #from B2E edge = [el, lf, bgroup]
            U1 = []; U2 = []; U3 = []; U4 = []
            for i in range(6):
                U1.append(Uel[4*i]) #rho
                U2.append(Uel[4*i + 1]) #rho U
                U3.append(Uel[4*i + 2]) #rho V
                U4.append(Uel[4*i + 3]) #rho E

            U1Q = quadstate(U1,edge[1],2,6)
            U2Q = quadstate(U2,edge[1],2,6)
            U3Q = quadstate(U3,edge[1],2,6)
            U4Q = quadstate(U4,edge[1],2,6)
            
            MT = 0
            for i in range(7):
                MT += wq[i]*nm(norms[i],[U1Q[i],U2Q[i],U3Q[i],U4Q[i]])*dsig[i] #include norm of stan
            edgeerror.append(L[i]**r*abs(MT))
            #print(edgeerror[-1])

    slist = sorted(enumerate(edgeerror), key=lambda x: x[1], reverse = True)
    slist = slist[0:round(len(slist)*0.03)]
    indexesint = {f: 0 for f in [i[0] for i in slist]} #flagged edges
    flags = {} #flag edges
    for i, edge in enumerate(IB):
        if i in indexesint:
            if edge[2] < 0:
                flags[edge[0]] = 0
            else:
                flags[edge[0]] = 0
                flags[edge[2]] = 0
    for i, edge in enumerate(I2E):
        if edge[0] in flags:
            flags[edge[2]] = 0
        if edge[2] in flags:
            flags[edge[0]] = 0
    return flags

def mach_plot(sol,E2N,C,fname):
    U = sol.reshape(-1, 4)
    M = []
    p = []
    node0 = np.array(C)
    elem0 = np.array(E2N)
    for i in range(len(U)): #iterate over all cells
        gamma = 1.4
        v = np.array([U[i,1]/U[i,0], U[i,2]/U[i,0]]) #flow velocity
        p.append((gamma-1)*(U[i,3] - 0.5*U[i,0]*(v[0]**2 + v[1]**2))) #pressure
        c = np.sqrt(p[i]*gamma/U[i,0]) #speed of sound
        M.append(np.sqrt(v[0]**2+v[1]**2)/c) #Mach number

    plt.tripcolor(node0[:,0], node0[:,1], triangles=elem0-1, facecolors=M, shading='flat',vmin = 0, vmax = 1)
    #plt.triplot(node0[:,0], node0[:,1], elem0-1,'k-',linewidth=1)
    plt.set_cmap('jet')
    plt.axis('equal')
    plt.xlim(-1,2.5)
    plt.ylim(-1,1.5)
    plt.colorbar()
    plt.xlabel('X - Coordinate'); plt.ylabel('Y - Coordinate'); plt.title('Contour Plot of Mach Number')
    plt.savefig('pictures/mach' + fname + '.pdf')
    

def loadtext(file):
    f = open(file,'r')
    vec = []
    for line in f:
        v = line.split()
        vv = []
        for l in v:
            vv.append(float(l))
        vec.append(vv)
    return vec

def savetext(name,vec):
    with open(name, 'w') as file:
        for row in vec:
            row_string = ' '.join(str(element) for element in row)
            file.write(row_string + '\n')

def flagerror(flags,epsilon,Area,met,comp,inp):
    if met != 'res':
        Area = np.ones((len(Area)))
    eps = abs(epsilon)
    EE = np.zeros((int(len(epsilon)/4)))
    c = len(EE)
    for i in range(int(len(EE)/4)): #loop over coarse space elements
        EE[4*i:4*(i+1)] = eps[4*i:4*(i+1)] + eps[c:c+4] + eps[c+4:c+8] + eps[c+8:c+12]
        c += 12

    error = []
    if comp == "rho":
        error = abs(EE[0::4])#/np.array(Area)
    elif comp == "xmom":
        error = abs(EE[1::4])#/np.array(Area)
    elif comp == "ymom":
        error = abs(EE[2::4])#/np.array(Area)
    elif comp == "rhoe":
        error = abs(EE[3::4])#/np.array(Area)
    if comp == "comb":
        for i in range(len(Area)):
            error.append((abs(EE[4*i])+abs(EE[4*i+1])+abs(EE[4*i+2])+abs(EE[4*i+3])))#Area[i]) #error indicator for all state components of each element

    slist = sorted(enumerate(error), key=lambda x: x[1], reverse = True) #element indeces sorted in descending order of error
    slist = slist[0:round(len(slist)*0.08)]
    flags = {f + 1: 0 for f in [i[0] for i in slist]} #flagged edges
    return flags