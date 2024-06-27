import numpy as np

def TriLagrange(x,y,order): #create phi matrix component for singular 1d quad point
    if order == 2:
        phi = np.zeros((6,1))
        phi[0] = 1.0-3.0*x-3.0*y+2.0*x*x+4.0*x*y+2.0*y*y
        phi[1] = 4.0*x-4.0*x*x-4.0*x*y
        phi[2] = -x+2.0*x*x
        phi[3] = 4.0*y-4.0*x*y-4.0*y*y
        phi[4] = 4.0*x*y
        phi[5] = -y+2.0*y*y
    else:
        print('basis order not supported')
    return phi

def makephimat(xref):
    phi = np.zeros((len(xref),6))
    for i, xr in enumerate(xref):
        phi[i,:] = np.matrix.transpose(TriLagrange(xr[0],xr[1],2))
    return phi

def GradTriLagrange(x,y,order): #create gradient of phi matrix
    if order == 2:
        phix = np.zeros((6,1))
        phix[0] =  -3.0+4.0*x+4.0*y
        phix[1] =  4.0-8.0*x-4.0*y
        phix[2] =  -1.0+4.0*x
        phix[3] =  -4.0*y
        phix[4] =  4.0*y
        phix[5] =  0.0

        phiy = np.zeros((6,1))
        phiy[0] =  -3.0+4.0*x+4.0*y
        phiy[1] =  -4.0*x
        phiy[2] =  0.0
        phiy[3] =  4.0-4.0*x-8.0*y
        phiy[4] =  4.0*x
        phiy[5] =  -1.0+4.0*y
    else:
        print('basis gradient order not supported')
    return phix, phiy

def ElemJacobian(xref, Q, xyQ): #xyQ = coordinates in global space
    [phix, phiy] = GradTriLagrange(xref[0],xref[1],Q) #x, y coordinate of point in ref elem, Q = order of geom.
    J = np.zeros((2,2))
    for i in range(len(xyQ)):
        J[0,0] += xyQ[i,0]*phix[i]
        J[0,1] += xyQ[i,0]*phiy[i]
        J[1,0] += xyQ[i,1]*phix[i]
        J[1,1] += xyQ[i,1]*phiy[i]
    return J

def quad1d():
    # Order 13 Gauss-Legendre points
    xq = [0.025446043828621, 0.129234407200303, 0.297077424311301, 0.500000000000000, \
    0.702922575688699, 0.870765592799697, 0.974553956171379]
    wq = [0.064742483084435, 0.139852695744638, 0.190915025252560, 0.208979591836735, \
    0.190915025252560, 0.139852695744638, 0.064742483084435]
    return xq, wq

def lagnorm(el,lf,CC): #el is [n1, n2, n3, n4, n5, n6], lf = localface, CC is the coordinates of all points
    xyQ = np.zeros((6,2)); Q = 2
    for i in range(6):
        for j in range(2):
            xyQ[i,j] = CC[int(el[i])-1][j] #make global space coords
    if lf == 1:
        xi_sigma = 1; eta_sigma = 0
    if lf == 2:
        xi_sigma = -1; eta_sigma = 1
    if lf == 3:
        xi_sigma = 0; eta_sigma = -1
    nvec = []
    [xq,wq] = quad1d()
    xref = RefEdge2RefElem(lf,xq)
    
    for xr in xref:
        J = ElemJacobian(xr, Q, xyQ)
        stan = J[:,0]*xi_sigma + J[:,1]*eta_sigma
        nvec.append([stan[1], -stan[0]])
    return np.array(nvec)

def RefEdge2RefElem(edge, xedge): #generates xref from xq from quad1d()
    sigma = np.zeros((len(xedge), 1))
    sigma[:,0] = xedge
    Z = np.zeros((len(sigma),1))
    O = np.ones((len(sigma),1))
    if edge == 1:
        xref = np.concatenate((sigma, Z), axis = 1)
    elif edge == 2:
        xref = np.concatenate((O-sigma, sigma), axis = 1)
    elif edge == 3:
        xref = np.concatenate((Z, O-sigma), axis = 1)
    else:
        print('edge out of bounds')
    return xref

def quadstate(state,lf): #state and localface
    [xq, wq] = quad1d()
    xref = RefEdge2RefElem(lf, xq)
    phi = makephimat(xref)
    Q = np.dot(phi,state)
    return Q

def inheritstatep0p1(U0):
    U1 = np.zeros((len(U0)*3))
    for i in range(int(len(U0)/4)): #for each elem in U0
        for j in range(3):
            U1[12*i + 4*j:12*(i+1) - (2-j)*4] = U0[4*i:4*(i+1)]
    return U1

def inheritstatep1p2(U1):
    U2 = np.zeros((len(U1)*2))
    for i in range(int(len(U1)/12)): #for elem node in U1
        U = U1[12*i:12*(i+1)]
        u1 = U[0:4]
        u3 = U[4:8]
        u6 = U[8:12]
        u2 = np.mean([u1, u3], axis=0)
        u4 = np.mean([u1, u6], axis=0)
        u5 = np.mean([u3, u6], axis=0)
        U2[24*i:24*(i+1)] = np.concatenate([u1,u2,u3,u4,u5,u6],axis=0)
    return U2

def inheritstatep2p3(U2):
    U3 = np.zeros((int(len(U2)*(5/3))))
    for i in range(int(len(U2)/24)): #for elem node in U1
        U = U2[24*i:24*(i+1)]
        u1 = U[0:4]
        u4 = U[8:12]
        u10 = U[20:24]
        mid1 = U[4:8]
        mid2 = U[16:20]
        mid3 = U[12:16]
        u2 = np.mean([u1, mid1], axis=0)
        u3 = np.mean([mid1, u4], axis=0)
        u5 = np.mean([mid3, u1], axis=0)
        u7 = np.mean([mid2, u4], axis=0)
        u8 = np.mean([mid3, u10], axis=0)
        u9 = np.mean([mid2, u10], axis=0)
        u6 = np.mean([u2, u5, u7, u9 ], axis=0)
        U3[40*i:40*(i+1)] = np.concatenate([u1,u2,u3,u4,u5,u6,u7,u8,u9,u10],axis=0)
    return U3