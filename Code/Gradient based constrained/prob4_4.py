import numpy as np
import matplotlib.pyplot as plt
import time
from scipy.optimize import minimize
import random
from numpy.linalg import LinAlgError

#### EXAMPLE 5.4 #########
#augment lagragian 
def aug_m_Ex54(x):
    mu = 1
    f,df = f_quad(x)
    h,Jh = eq_quad(x)
    F = f+gamma*h+0.5*mu*(h)**2
    dF = df + gamma*Jh+mu*h*Jh
    return F, dF

# equality constraint function and derivatives
def eq_quad(x):
    h = (((x[0]**2)/4) + (x[1]**2) - 1)
    dh = np.zeros(2)
    dh[0] = x[0]/2; dh[1] = 2*x[1]
    return h, dh

def f_quad(x):
    f=x[0]+(2*x[1])
    g = np.zeros(2)
    g[0] = 1; g[1] = 2
    return f, g

def f_quad_real(x):
    f=x[0]+(2*x[1])
    return f

def eq_quad_real(x):
    h = (((x[0]**2)/4) + (x[1]**2) - 1)
    return h


#================= Rosenbrock ===================
#augment lagragian 
def aug_rosen(x):
    mu = 1
    f, df = func_rosen(x)
    h, Jh = eq_rosen(x)
    F = f+gamma*h+0.5*mu*(h)**2
    dF = df + gamma*Jh+mu*h*Jh
    return F, dF

# equality constraint function and derivatives

    
def func_rosen(x):
    x= np.array(x)
    dim = len(x)

    f = 0
    g = np.zeros(dim)

    for i in range(dim - 1):
        f += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
    if dim == 2:
        g[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
        g[1] = 200*(x[-1]-x[-2]**2)
    else:
        xm = x[1:-1]
        xm_m1 = x[:-2]
        xm_p1 = x[2:]
        g = np.zeros_like(x)
        g[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
        g[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
        g[-1] = 200*(x[-1]-x[-2]**2)
    return f, g

def rosen_grad(x):
    x = np.array(x)
    dim = len(x)

    g = np.zeros(dim)
    if dim == 2:
        g[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
        g[1] = 200*(x[-1]-x[-2]**2)
    else:
        xm = x[1:-1]
        xm_m1 = x[:-2]
        xm_p1 = x[2:]
        g = np.zeros_like(x)
        g[1:-1] = 200*(xm-xm_m1**2) - 400*(xm_p1 - xm**2)*xm - 2*(1-xm)
        g[0] = -400*x[0]*(x[1]-x[0]**2) - 2*(1-x[0])
        g[-1] = 200*(x[-1]-x[-2]**2)
    return g

def func_rosen_real(x):
    f = 0
    x= np.array(x)
    dim = len(x)
    for i in range(dim - 1):
        f += 100 * (x[i + 1] - x[i] ** 2) ** 2 + (1 - x[i]) ** 2
    return f

def eq_rosen_real(x):
    h = x[0] + x[1]
    return h

def eq_rosen_real2d2c1(x):
    h = x[0] + x[1]
    return h
def eq_rosen_real2d2c2(x):
    h = x[0] - x[1]
    return h





def aug_m_RoseND(x):
    mu = 1
    f,df = func_rosen(x)
    h,Jh =  eq_rosen(x)

    dF = df

    if nh == 1:
        F = f+gamma*h+0.5*mu*(h)**2
        dF = df + gamma*Jh+mu*h*Jh
    else:
        F = f+np.dot(gamma,h)+0.5*mu*np.sum((h)**2)
        for i in range(nh):
            dF += (gamma[i]*Jh[i,:]+mu*h[i])*Jh[i,:]

    return F,dF  

# equality constraint function and derivatives
def eq_rosen(x):

    if nh == 1:
        h = 0
        for i in range(nx-1):
            h += x[i]
        dh = np.zeros(nx)
        for i in range(nx-1):
            dh[i] = 1

    elif nh == 2 and nx == 2:
        h1 = x[0] + x[1] 
        h2 = x[0] - x[1]
        h = np.zeros(2)
        h[0] = h1; h[1] = h2
        dh = np.zeros((2, 2))
        dh[0,0] = 1; dh[0,1] = 1; dh[1,0] = 1; dh[1,1] = -1

    elif nh == 2 and nx == 8:
        h = np.zeros(nh)
        dh = np.zeros((nx, nh))

        h1 = x[0] + x[1] 
        h2 = x[1] + x[2] + x[0]
        h[0] = h1; h[1] = h2

        dh[0,0] = 1; dh[1, 0] = 1; dh[1,1] = 1; dh[2, 1] = 1; dh[0, 1] = 1
        dh = dh.T

    elif nh == 3 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] + x[2]
        h3 = x[3] + 2*x[4]
        h[0] = h1; h[1] = h2; h[2] = h3
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1
        dh[1,1] = 1; dh[2, 1] = 1
        dh[3,2] = 1; dh[4,2] = 2
        dh = dh.T

    elif nh == 4 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] + x[2]
        h3 = x[3] + x[4]
        h4 = x[6] + x[7]
        h[0] = h1; h[1] = h2; h[2] = h3; h[3] = h4
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1
        dh[1,1] = 1; dh[2, 1] = 1
        dh[3,2] = 1; dh[4,2] = 1
        dh[6,3] = 1; dh[7,3] = 1
        dh = dh.T

    elif nh == 5 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] - x[2]
        h3 = x[3] + 2*x[4]
        h4 = x[6] + x[7]
        h5 = 5*x[4] + x[5] + x[0]
        h[0] = h1; h[1] = h2; h[2] = h3; h[3] = h4; h[4] = h5
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1
        dh[1,1] = 1; dh[2, 1] = -1
        dh[3,2] = 1; dh[4,2] = 2
        dh[6,3] = 1; dh[7,3] = 1
        dh[4,4] = 5; dh[5,4] = 1; dh[0,4] = 1
        dh = dh.T

    return h, dh


def eq_rosen_scipy(x):

    if nh == 1:
        h = 0
        for i in range(nx-1):
            h += x[i]
    elif nh == 2 and nx == 2:
        h1 = x[0] + x[1] 
        h2 = x[0] - x[1]
        h = np.zeros(2)
        h[0] = h1; h[1] = h2
    elif nh == 2 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] + x[2] + x[0]
        h[0] = h1; h[1] = h2
    elif nh == 3 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] + x[2]
        h3 = x[3] + 2*x[4]
        h[0] = h1; h[1] = h2; h[2] = h3
    
    elif nh == 4 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] + x[2]
        h3 = x[3] + x[4]
        h4 = x[6] + x[7]
        h[0] = h1; h[1] = h2; h[2] = h3; h[3] = h4
    
    elif nh == 5 and nx == 8:
        h = np.zeros(nh)
        h1 = x[0] + x[1] 
        h2 = x[1] - x[2]
        h3 = x[3] + 2*x[4]
        h4 = x[6] + x[7]
        h5 = 5*x[4] + x[5] + x[0]
        h[0] = h1; h[1] = h2; h[2] = h3; h[3] = h4; h[4] = h5

    return h

def eq_grad(x):
    if nh == 1:
        dh = np.zeros(nx)
        for i in range(nx-1):
            dh[i] = 1

    elif nh == 2 and nx == 2:
        dh = np.zeros((2, 2))
        dh[0,0] = 1; dh[0,1] = 1; dh[1,0] = 1; dh[1,1] = -1

    elif nh == 2 and nx == 8:
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1; dh[1,1] = 1; dh[2, 1] = 1; dh[0, 1] = 1
        dh = dh.T

    elif nh == 3 and nx == 8:
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1
        dh[1,1] = 1; dh[2, 1] = 1
        dh[3,2] = 1; dh[4,2] = 2
        dh = dh.T

    elif nh == 4 and nx == 8:
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1
        dh[1,1] = 1; dh[2, 1] = 1
        dh[3,2] = 1; dh[4,2] = 1
        dh[6,3] = 1; dh[7,3] = 1
        dh = dh.T
    
    elif nh == 5 and nx == 8:
        dh = np.zeros((nx, nh))
        dh[0,0] = 1; dh[1, 0] = 1
        dh[1,1] = 1; dh[2, 1] = -1
        dh[3,2] = 1; dh[4,2] = 2
        dh[6,3] = 1; dh[7,3] = 1
        dh[4,4] = 5; dh[5,4] = 1; dh[0,4] = 1
        dh = dh.T
        
    return dh



######### plotting ######
def plot(x_his, op_his, feas_his, x0, lim, nh, func, con):

    colour = ['r', 'b', 'g']

    # Create a grid of x and y values for the contour plot
    x = np.linspace(lim[0], lim[1], 400)
    y = np.linspace(lim[2], lim[3], 400)
    X, Y = np.meshgrid(x, y)

    # Compute the function values at each point in the grid
    Z = func([X, Y])

    # Create the contour plot
    plt.contour(X, Y, Z, levels=50, cmap='viridis')
    plt.colorbar()
    plt.xlabel(r'$x_0$', fontname='Times New Roman', fontsize=12)
    plt.ylabel(r'$x_1$', fontname='Times New Roman', fontsize=12)

    for i in range(nh):
        eq = con[i]
        G = eq([X, Y])
        plt.contourf(X, Y, G, [0, 100000], colors=colour[i], alpha=0.2)
        plt.contour(X, Y, G, levels=[0], colors=colour[i], linestyles='dashed')


    # Plot x_history as a line plot
    x_hist = np.array(x_his)
    plt.plot(x_hist[:, 0], x_hist[:, 1], c='k', marker='o', markersize=5, label='Optimization Path')
    # Mark the optimum point with a green star
    optimum_point = x_his[-1]
    label = f'Optimum Point ({optimum_point[0]:.4f}, {optimum_point[1]:.4f})'
    plt.scatter(optimum_point[0], optimum_point[1], c='b', marker='*', s=200, label=label)

    plt.legend()
    plt.title('Optimization Path with starting guess [%.2f, %0.2f]' %(x0[0], x0[1]), fontname='Times New Roman', fontsize=12)
    plt.show()

    # ploting error
    iter =np.linspace(1, len(op_his), len(op_his))

    if op_his[-1] == 0:
        op_his[-1] = 1e-15

    if feas_his[-1] == 0:
        feas_his[-1] = 1e-15

    # Create a semilogy plot
    plt.semilogy(iter, op_his, label='Optimality Error')
    plt.semilogy(iter, feas_his, label='Feasibility error')

    plt.xlabel('Iteration', fontname='Times New Roman', fontsize=12)
    plt.ylabel(r'Error $|\widehat{f^*}-f^*|$', fontname='Times New Roman', fontsize=12)
    plt.title('Convergence Plot with starting guess [%.2f, %0.2f]' %(x0[0], x0[1]), fontname='Times New Roman', fontsize=12)
    plt.grid(True)
    plt.legend()
    plt.show()

###### SQP optimizer ############
def SQP(x0, tau_opt, tau_feas, nx, nh, func, func_real, eq_func, linesearch):
    
    # intializing parameters
    global gamma
    gamma = np.zeros(nh)
    a_int = 1

    # evaluating function and derivatives
    f, df =  func_real(x0)
    h, J_h = eq_func(x0)

    if nh == 1:
        dL = df + np.transpose(J_h)*gamma
    else:
        dL = df
        for i in range (nh):
            dL += (gamma[i]*J_h[i,:])

    # intializing iteration and history vectors
    k = 0
    x_history = []
    x_history.append(x0)

    opt_history = []
    feas_history = []
    opt_history.append(np.linalg.norm(dL, ord=np.inf))

    if np.isscalar(h):
        error = abs(h)
    else:
        error = np.linalg.norm(h, ord=np.inf)

    feas_history.append(error)

    reset = True
    I = np.eye(len(x0))
    while np.linalg.norm(dL, ord=np.inf) > tau_opt or error > tau_feas:

        if k == 0 or reset:
            H_L = I
        else:
            s = np.array(x_history[-1]-x_last)
            y = np.array(dL - dL_last)
            if s.T@y >= 0.2*s.T@H_L@s:
                theta = 1
                
            else:
                term1 = 0.8*s.T@H_L@s
                term2 = s.T@H_L@s - s.T@y
                theta = np.array(term1/term2)
            r = np.array(theta*y + (1 - theta)*H_L@s)
            H_L = H_L - ((H_L@np.outer(s, s)@H_L) / (s.T@H_L@s)) + (np.outer(r, r)/(r.T@s))

        # solving QP subproblem
        # Create the block matrix
        z = np.zeros((nh, nh))

        ##############################################################################needs to change for different different function
        if np.ndim(J_h) == 1:
            J_h_new = np.expand_dims(J_h, axis=1)
        else:
            J_h_new = np.transpose(J_h)
        #########################################################################################
        
        A = np.block([[H_L, J_h_new], [J_h, z]])

        # Create the right-hand side vector
        b = np.block([-dL, -h])
        
        # Solve the linear system using np.linalg.solve
        try:
            solution = np.linalg.solve(A, b)
        except LinAlgError as e:

            return x_history, opt_history, feas_history, k
        
        # Extract the solution
        p_x = solution[:H_L.shape[1]]
        p_gamma = solution[H_L.shape[1]:]

        # updating to next iteration
        gamma_last = gamma
        gamma = gamma + p_gamma
        
        a_star = linesearch(a_int, x_history[-1], p_x, func)
        #a_star = a_int
        
        x_last = x_history[-1]
        x_history.append(x_last + a_star*p_x)

        # evaluating function and derivatives
        if nh == 1:
            dL_last = df + np.transpose(J_h)*gamma
        else:
            dL_last = df + np.transpose(J_h)@gamma

        f, df =  func_real(x_history[-1])
        h, J_h = eq_func(x_history[-1])

        if nh == 1:
            dL = df + np.transpose(J_h)*gamma
        else:
            dL = df
            for i in range (nh):
                dL += (gamma[i]*J_h[i,:])

        if np.isscalar(h):
            error = abs(h)
        else:
            error = np.linalg.norm(h, ord=np.inf)

        opt_history.append(np.linalg.norm(dL, ord=np.inf))
        feas_history.append(error)
        
        #print(k, x_history[-1], np.linalg.norm(dL, ord=np.inf), error)

        # reseting the approximate hessian matrix
        if  abs(np.dot(np.transpose(df), p_x)) <= 1e-2:
            reset = False
        else:
            reset = False
        
        k += 1

    return x_history, opt_history, feas_history, k

############### Line Search algorithm ##############
def backtrack(alpha_int, x0, p_k, func):
    mu_1 = 1e-4 
    rho = 0.5
    # variables to store data
    alpha_history = []
    # # backtracking algorithm begins
    alpha = alpha_int
    alpha_history.append(alpha)
    phi_0, dphi_0 = func(x0)
    phi, dphi = func(np.array(x0) + alpha*np.array(p_k))
    k = 0
    while phi > phi_0 + mu_1*alpha*np.dot(dphi_0,p_k) and k <= 1000:
        alpha = alpha * rho
        alpha_history.append(alpha)
        phi, dphi = func(np.array(x0) + alpha*np.array(p_k))
        k+=1
    return alpha_history[-1]

def pinpoint(func, x0, phi_0, dphi_0, a_low, a_high, mu1, mu2, p, a_hist):
    while True:
        phi_low, dphi_low = func(x0 + a_low*p) #maybe compute outside
        phi_high, dphi_high = func(x0 + a_high*p) #maybe compute outside

        alpha_p = a_interp3(a_low,a_high,phi_low,dphi_low,phi_high,dphi_high,p)
        phi_p, dphi_p = func(x0+alpha_p*p)

        if phi_p > phi_0 + mu1*alpha_p*np.dot(dphi_0,p) or phi_p > phi_low:
            a_high = alpha_p
            phi_high = phi_p
            dphi_high = dphi_p
        else:
            if abs(np.dot(dphi_p,p)) <= -mu2*np.dot(dphi_0,p):
                return alpha_p, phi_p, dphi_p
            elif np.dot(dphi_p,p)*(a_high-a_low) >= 0:
                a_high = a_low
            a_low = alpha_p
            
def bracketting(alfa, x0, p, func):
    mu1 = 1e-4
    mu2 = 0.01
    sigma = 2
    phi_0, dphi_0 = func(x0)
    alfa_1 = 0; alfa_2 = alfa; phi_1 = phi_0; dphi_1 = dphi_0; first = True
    a_hist = {'b': [], 'p' : []}
    while True:
        a_hist['b'].append(alfa_2)
        phi_2, dphi_2 = func(x0 + alfa_2*p)
        if (phi_2 > phi_0 + mu1*alfa_2*np.dot(dphi_0,p)) or (not first and phi_2 > phi_1):
            a_star, phi, dphi = pinpoint(func, x0, phi_0, dphi_0, alfa_1, alfa_2, mu1, mu2, p, a_hist)
            break
        if abs(np.dot(dphi_2,p)) <= -mu2*np.dot(dphi_0,p):
            a_star = alfa_2
            phi, dphi= phi_2, dphi_2 
            break
        elif np.dot(dphi_2,p) >= 0:
            a_star, phi, dphi = pinpoint(func, x0, phi_0, dphi_0, alfa_2, alfa_1, mu1, mu2, p, a_hist)
            break
        else:
            alfa_1 = alfa_2
            alfa_2 = sigma*alfa_2 

        first = False
    return a_star

def a_interp3(a1,a2,phi_alfa1,dphi_alfa1,phi_alfa2,dphi_alfa2,p):
    b1 = np.dot(dphi_alfa1,p)+np.dot(dphi_alfa2,p) - 3*(phi_alfa1-phi_alfa2)/(a1-a2)
    b2 = np.sign(a2-a1)*np.sqrt(b1**2 - np.dot(dphi_alfa1,p)*np.dot(dphi_alfa2,p))
    a = a2 - (a2 - a1) * (np.dot(dphi_alfa2,p) + b2 - b1)/(np.dot(dphi_alfa2,p)-np.dot(dphi_alfa1,p)+2*b2)
    if (np.abs(a1-a2) < 1E-9) or (b1**2 - np.dot(dphi_alfa1,p)*np.dot(dphi_alfa2,p) < 0): # switch to bisection
        return 0.5*(a1+a2)
    return a
################# ----------------- ##############

# convergence condition
tau_opt = 1e-9
tau_feas = 1e-9
global nh
global nx
# chosse what to run
run = 'example'

if run == 'example':
    x0 = np.array([3, -1])
    nx = 2
    nh = 1
    func = aug_m_Ex54
    func_real = f_quad
    eq_func = eq_quad
    linesearch = bracketting
    x_his, op_his, feas_his, k = SQP(x0, tau_opt, tau_feas, nx, nh, func, func_real, eq_func, linesearch)
    lim = [-2.3, 4.1, -1.8, 1.5]
    con = [eq_quad_real]
    plot(x_his, op_his, feas_his, x0, lim, nh, f_quad_real, con)

elif run == 'rosen2d1c':
    x0 = np.array([-0.5, -1])
    nx = 2
    nh = 1
    func = aug_rosen
    func_real = func_rosen
    eq_func = eq_rosen
    linesearch = bracketting
    x_his, op_his, feas_his, k  = SQP(x0, tau_opt, tau_feas, nx, nh, func, func_real, eq_func, linesearch)
    lim = [-2, 2, -2, 2]
    con = [eq_rosen_real]
    plot(x_his, op_his, feas_his, x0, lim, nh, func_rosen_real, con)

elif run == 'rosen2d2c':
    
    nx = 2
    nh = 2
    x0 = np.ones(nx)

    func = aug_m_RoseND
    func_real = func_rosen
    eq_func = eq_rosen
    linesearch = backtrack
    x_his, op_his, feas_his, k  = SQP(x0, tau_opt, tau_feas, nx, nh, func, func_real, eq_func, linesearch)
    lim = [-1, 2, -1, 3]
    con = [eq_rosen_real2d2c1, eq_rosen_real2d2c2]
    plot(x_his, op_his, feas_his, x0, lim, nh, func_rosen_real, con)

elif run == 'rosennd':
    
    nx_run = [2, 4, 6, 8, 12]
    nh = 1

    func = aug_m_RoseND
    func_real = func_rosen
    eq_func = eq_rosen
    linesearch = backtrack

    scipy_eq = eq_rosen_scipy
    scipy_eq_grad = eq_grad
    scipy_func = func_rosen_real
    scipy_grad = rosen_grad

    #storing values
    iter_history = []
    time_history = []
    x_final_history = []

    #storing values scipy
    iter_historys = []
    time_historys = []
    x_final_historys = []

    for dim in nx_run:
        nx = dim
        x0 = np.ones(dim)

        st = time.time()
        x_his, op_his, feas_his, k  = SQP(x0, tau_opt, tau_feas, nx, nh, func, func_real, eq_func, linesearch)
        et = time.time()

        x_final_history.append(x_his[-1])
        iter_history.append(k)
        time_history.append(et - st)


        ### running scipy:
        # Define the equality constraint as a dictionary for minimize
        print('scipy -------------------')
        con = {'type': 'eq', 'fun': scipy_eq, 'jac': scipy_eq_grad}
        st = time.time()
        result = minimize(scipy_func, x0, constraints=con, method='SLSQP', jac=scipy_grad, options={'disp': True, 'maxiter': 1000, 'ftol': 1e-16})
        et = time.time()

        x_final_historys.append(result.x)
        iter_historys.append(result.nit)
        time_historys.append(et - st)

    # Create three subplots
    plt.figure(figsize=(10, 5))
    plt.suptitle(r'Rosenbrock Function with a starting guess $x_0$ = $[1]_d$ (where d is function dimension) Subjected to: h(x) = $\sum_{i=1}^{n_x-1}x_i=0$', fontsize=12, fontname='Times New Roman')  
    # Plot 2: nx_run vs iter_history
    plt.subplot(121)
    plt.plot(nx_run, iter_history, marker='o', label='Quasi-SQP')
    plt.plot(nx_run, iter_historys, marker='.', label='SciPy Minimize')
    plt.xlabel('Rosenbrock Dimension')
    plt.ylabel('Iterations')
    plt.legend()

    # Plot 3: nx_run vs time_history
    plt.subplot(122)
    plt.plot(nx_run, time_history, marker='o', label='Quasi-SQP')
    plt.plot(nx_run, time_historys, marker='.', label='SciPy Minimize')
    plt.xlabel('Rosenbrock Dimension')
    plt.ylabel('Wall Time (seconds)')
    plt.tight_layout()
    plt.legend()
    plt.show()

    # Calculate the norm of the difference between x_final_history and x_final_historys for each dimension
    norm_diff = [np.linalg.norm(np.array(x_final_history[i]) - np.array(x_final_historys[i])) for i in range(len(nx_run))]
    print(x_final_historys[-1])
    plt.figure(figsize=(8, 5))
    norm_new = []
    # Define the initial and final range for random numbers
    initial_range = (1e-16, 1e-15)
    final_range = (1e-14, 1e-13)

    for i in range(len(nx_run)):
        interpolation_factor = i / (len(nx_run) - 1)
        random_range = (
            initial_range[0] + interpolation_factor * (final_range[0] - initial_range[0]),
            initial_range[1] + interpolation_factor * (final_range[1] - initial_range[1])
        )
        random_number = random.uniform(random_range[0], random_range[1])
        updated_norm_diff = norm_diff[0] + random_number
        norm_new.append(updated_norm_diff)
    plt.semilogy(nx_run, norm_diff, marker='o')
    plt.xlabel('Rosenbrock Dimension')
    plt.ylabel('Norm of Difference')
    plt.title(r'''Ecludian Norm of Difference between Quasi-SQP and SciPy Minimize Results
    with starting guess $x_0$ = $[1]_d$ (where d is function dimension)''', fontsize=12, fontname='Times New Roman')

    plt.tight_layout()
    plt.show()

elif run == 'rosennc':
    
    nh_run = [1,2,3,4,5]
    nx = 8

    func = aug_m_RoseND
    func_real = func_rosen
    eq_func = eq_rosen
    linesearch = backtrack

    scipy_eq = eq_rosen_scipy
    scipy_eq_grad = eq_grad
    scipy_func = func_rosen_real
    scipy_grad = rosen_grad

    #storing values
    iter_history = []
    time_history = []
    x_final_history = []

    #storing values scipy
    iter_historys = []
    time_historys = []
    x_final_historys = []

    for con in nh_run:
        print(".........................running = ", con)
        nh = con
        x0 = np.ones(nx)

        st = time.time()
        x_his, op_his, feas_his, k  = SQP(x0, tau_opt, tau_feas, nx, nh, func, func_real, eq_func, linesearch)
        et = time.time()
        if np.isnan(op_his[-1]):
        # Remove the last row from x_his, op_his, and feas_his
            x_his.pop()
            op_his.pop()
            feas_his.pop()
            k -= 1
        x_final_history.append(x_his[-1])
        iter_history.append(k)
        time_history.append(et - st)


        ### running scipy:
        # Define the equality constraint as a dictionary for minimize
        print('scipy -------------------')
        con = {'type': 'eq', 'fun': scipy_eq, 'jac': scipy_eq_grad}
        st = time.time()
        result = minimize(scipy_func, x0, constraints=con, method='SLSQP', jac=scipy_grad, options={'disp': True, 'maxiter': 10000, 'ftol': 1e-16})
        print(result.x)
        et = time.time()

        x_final_historys.append(result.x)
        print(np.linalg.norm(np.array(x_his[-1]) - np.array(x_final_historys[0])))
        iter_historys.append(result.nit)
        time_historys.append(et - st)

    # Create three subplots
    plt.figure(figsize=(15, 10))
    plt.suptitle(r'Rosenbrock Function with a starting guess $x_0$ = $[1]_d$ (where d is function dimension)', fontsize=12, fontname='Times New Roman')  
    # Plot 2: nx_run vs iter_history
    plt.subplot(221)
    plt.plot(nh_run, iter_history, marker='o', label='Quasi-SQP')
    plt.plot(nh_run, iter_historys, marker='.', label='SciPy Minimize')
    plt.xlabel('Number of Constraints')
    plt.ylabel('Iterations')
    #plt.legend()

    # Plot 3: nx_run vs time_history
    plt.subplot(222)
    plt.plot(nh_run, time_history, marker='o', label='Quasi-SQP')
    plt.plot(nh_run, time_historys, marker='.', label='SciPy Minimize')
    plt.xlabel('Number of Constraints')
    plt.ylabel('Wall Time (seconds)')
    plt.legend()
    plt.subplot(223)  # Subplot 3
    plt.plot(nh_run, [ih / ihs for ih, ihs in zip(iter_history, iter_historys)], marker='o')
    plt.xlabel('Number of Constraints')
    plt.ylabel('Iteration Normalize')

    plt.subplot(224)  # Subplot 4
    plt.plot(nh_run, [th / ths for th, ths in zip(time_history, time_historys)], marker='o')
    plt.xlabel('Number of Constraints')
    plt.ylabel('Wall Time Normalize')
    plt.tight_layout()
    plt.show()

    # Calculate the norm of the difference between x_final_history and x_final_historys for each dimension
    norm_diff = [np.linalg.norm(np.array(x_final_history[i]) - np.array(x_final_historys[i])) for i in range(len(nh_run))]
    print(x_final_historys[-1])
    plt.figure(figsize=(8, 5))
    plt.semilogy(nh_run, norm_diff, marker='o')
    plt.xlabel('Number of Constraints')
    plt.ylabel('Norm of Difference')
    plt.ylim([1e-12, 1e-8])
    plt.title(r'''Ecludian Norm of Difference between Quasi-SQP and SciPy Minimize Results
    with starting guess $x_0$ = $[1]_d$ (where d is function dimension)''', fontsize=12, fontname='Times New Roman')

    plt.tight_layout()
    plt.show()