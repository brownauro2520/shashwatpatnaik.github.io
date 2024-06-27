import numpy as np
import sys
from funcs import *

def main():
    gri = 'meshes/c.gri'
    sol = 'sols/U.txt'

    U = np.loadtxt(sol)
    [E2N,C,B] = readgri(gri) #read gri
    [I2E,B2E,In,Bn,Il,Bl,Area,L,N,bcords] = makematrices(E2N,B,C) #make matrices

    outputjac(E2N,B2E,Bn,Bl,U,C)

def outputjac(E2N,B2E,Bn,Bl,U,C):
    rhoinf = 1
    uinf = 0.5
    gamma = 1.4
    DLDU = np.zeros((len(E2N)*4,1))
    DDDU = np.zeros((len(E2N)*4,1))
    DMDU = np.zeros((len(E2N)*4,1))
    DSDU = np.zeros((len(E2N)*4,1))
    M = 0
    Fy = 0; Fx = 0
    for i, B in enumerate(B2E): #[element nr, localface, bgroup]
        if B[2] != -1:
            E = B[0] - 1 #index
            els = E2N[E]
            if B[1] == 1:
                xmid = (C[els[0]-1][0] + C[els[1]-1][0])/2
                ymid = (C[els[0]-1][1] + C[els[1]-1][1])/2
            elif B[1] == 2:
                xmid = (C[els[1]-1][0] + C[els[2]-1][0])/2
                ymid = (C[els[1]-1][1] + C[els[2]-1][1])/2
            else:
                xmid = (C[els[2]-1][0] + C[els[0]-1][0])/2
                ymid = (C[els[2]-1][1] + C[els[0]-1][1])/2

            
            u1 = U[4*E]; u2 = U[4*E+1]; u3 = U[4*E+2]; u4 = U[4*E+3]
            nx = Bn[i][0]; ny = Bn[i][1]
            vbx = u2/u1 - (nx*u2/u1 + ny*u3/u1)*nx
            vby = u3/u1 - (nx*u2/u1 + ny*u3/u1)*ny
            vbnorm = np.sqrt(vbx**2 + vby**2)
            dp1 = (gamma-1)*((Bl[i]/(0.5*rhoinf*uinf**2))*(-u1*np.sqrt((((nx**2-1)*u2 + nx*ny*u3)/(u1**2))**2 + ((u2*nx*ny + (ny**2-1)*u3)/(u1**2))**2))*vbnorm - 0.5*vbnorm**2)
            dp2 = (gamma-1)*(Bl[i]/(0.5*rhoinf*uinf**2))*(-u1*np.sqrt(((1-nx**2)/(u1))**2 + ((-nx*ny)/(u1))**2))*vbnorm
            dp3 = (gamma-1)*(Bl[i]/(0.5*rhoinf*uinf**2))*(-u1*np.sqrt((-nx*ny/u1)**2 + ((1-ny**2)/(u1))**2))*vbnorm
            dp4 = (gamma-1)*(Bl[i]/(0.5*rhoinf*uinf**2))*(1)
            DLDU[4*E] = ny*dp1
            DLDU[4*E+1] = ny*dp2
            DLDU[4*E+2] = ny*dp3
            DLDU[4*E+3] = ny*dp4

            DDDU[4*E] = nx*dp1
            DDDU[4*E+1] = nx*dp2
            DDDU[4*E+2] = nx*dp3
            DDDU[4*E+3] = nx*dp4

            DMDU[4*E] = dp1*(-nx*ymid-ny*xmid)
            DMDU[4*E+1] = dp2*(-nx*ymid-ny*xmid)
            DMDU[4*E+2] = nx*(-nx*ymid-ny*xmid)
            DMDU[4*E+3] = nx*(-nx*ymid-ny*xmid)
            vb = [u2/u1, u3/u1] - np.dot([u2/u1, u3/u1],Bn[i])*np.array(Bn[i])
            fy = ny*Bl[i]*(gamma-1)*(u4-0.5*u1*(vb[0]**2 + vb[1]**2))/(0.5*rhoinf*uinf**2)
            fx = nx*Bl[i]*(gamma-1)*(u4-0.5*u1*(vb[0]**2 + vb[1]**2))/(0.5*rhoinf*uinf**2)
            Fy += fy
            Fx += fx
            M += fx*ymid -fy*(xmid-0.25)
    
    for i in range(0,len(E2N)):
        u1 = U[4*i]; u2 = U[(4*i)+1]; u3 = U[(4*i)+2]; u4 = U[(4*i)+3]
        DSDU[4*i] = (u1*(gamma/u1 - (0.5000*u2**2 + 0.5000*u3**2)/(u1**2*(u4 - (0.5000*u2**2 + 0.5000*u3**2)/u1))))/(gamma - 1) - (math.log((u4 - (0.5000*u2**2 + 0.5000*u3**2)/u1)*(gamma - 1)) - gamma*math.log(u1))/(gamma - 1)
        DSDU[(4*i)+1] = -(2*u1*u2)/((gamma - 1)*(u2**2 + u3**2 - 2*u1*u4))
        DSDU[(4*i)+2] = -(2*u1*u3)/((gamma - 1)*(u2**2 + u3**2 - 2*u1*u4))
        DSDU[(4*i)+3] = (2*u1**2)/((gamma - 1)*(u2**2 + u3**2 - 2*u1*u4))

    alpha = -8*np.pi/180
    rot = np.array([[np.cos(alpha),-np.sin(alpha)],[np.sin(alpha),np.cos(alpha)]])
    L = np.dot(rot,np.array([Fx,Fy]))[1]
    D = np.dot(rot,np.array([Fx,Fy]))[0]
    print(L,D,M)
    np.savetxt('sols/outjacs/DLDU.txt',DLDU)
    np.savetxt('sols/outjacs/DDDU.txt',DDDU)
    np.savetxt('sols/outjacs/DMDU.txt',DMDU)
    np.savetxt('sols/outjacs/DSDU.txt',DSDU)

if __name__ == "__main__":
    main()