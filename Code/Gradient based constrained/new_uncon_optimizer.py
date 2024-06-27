"""
This is a template for Assignment 3: unconstrained optimization

You can (and should) call other functions or import functions from other files,
but make sure you do not change the function signature (i.e., function name `uncon_optimizer`, inputs, and outputs) in this file.
The autograder will import `uncon_optimizer` from this file. If you change the function signature, the autograder will fail.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def uncon_optimizer(func, x0, epsilon_g, options=None):
    output = {}

    if options is None:
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.28; sigma = 2; reset_point = 1e-8
        linesearch = backtrack
        namelinesearch="bactrack"
        searchdir = BFGS
        namesearchdir="quasi"

    # Direction Search Algorithm 
    if options == 'quasi_bt':
        mu_1 = 1e-2; mu_2 = 0.01; rho = 0.4; sigma = 2; reset_point = 1e-8
        searchdir = BFGS
        namesearchdir="quasi"
        linesearch = backtrack
        namelinesearch="backtracking"

    elif options == 'quasi_bp':
        mu_1 = 1e-4; mu_2 = 1e-2; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir = BFGS  
        namesearchdir="quasi"
        linesearch = bracketting
        namelinesearch="bracketing"
        
    xopt, fopt, k, x_his, error_his = searchdir(x0, epsilon_g, mu_1, mu_2, rho, sigma, reset_point, func, linesearch)
    output = {'alias':"browncurry2520", 'iteration':k,'x history': x_his, 'error history': error_his, 'namesearchdir':namesearchdir,'namelinesearch':namelinesearch}
    return xopt, fopt, output

# backtracking algorithm 
""" def backtrack(alpha_int, x0, mu_1, mu_2, p_k, rho, sigma, phi_0, dphi_0, func, inf):
    # variables to store data
    alpha_history = []
    # # backtracking algorithm begins
    alpha = alpha_int
    alpha_history.append(alpha)
    phi, dphi = func(np.array(x0) + alpha*np.array(p_k))
    k = 0
    while phi > phi_0 + mu_1*alpha*np.dot(dphi_0,p_k) and k <= 1000:
        alpha = alpha * rho
        alpha_history.append(alpha)
        phi, dphi = func(np.array(x0) + alpha*np.array(p_k))
        k+=1
    return alpha_history[-1], phi, dphi """

def backtrack(a0,x0,mu1,mu2,p,rho,sigma,phi0,dphi0,func,inf):
    alpha = a0 #set initial step size
    phi = 0; first = True
    while phi > phi0 + mu1*alpha*np.dot(dphi0,p) or first: #Function value is above sufficient decrease line, condition from algorithm in mdobook
        alpha *= rho #reduce step size, backtrack
        phi, dphi = func(x0 + alpha*p)
        first = False
    return alpha, phi, dphi

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
            
def bracketting(alfa, x0, mu1, mu2, p, rho, sigma, phi_0, dphi_0, func, inf):
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
    return a_star, phi, dphi 

def a_interp3(a1,a2,phi_alfa1,dphi_alfa1,phi_alfa2,dphi_alfa2,p):
    b1 = np.dot(dphi_alfa1,p)+np.dot(dphi_alfa2,p) - 3*(phi_alfa1-phi_alfa2)/(a1-a2)
    b2 = np.sign(a2-a1)*np.sqrt(b1**2 - np.dot(dphi_alfa1,p)*np.dot(dphi_alfa2,p))
    a = a2 - (a2 - a1) * (np.dot(dphi_alfa2,p) + b2 - b1)/(np.dot(dphi_alfa2,p)-np.dot(dphi_alfa1,p)+2*b2)
    if (np.abs(a1-a2) < 1E-9) or (b1**2 - np.dot(dphi_alfa1,p)*np.dot(dphi_alfa2,p) < 0): # switch to bisection
        return 0.5*(a1+a2)
    return a


# BFGS or Quasi-Newton Direction
def BFGS(x0, epsilon_g, mu_1, mu_2, rho, sigma, reset_point, func, line_search):
    #maximum iteration
    maxk = 300000
    k = 0
    reset = True
    alpha_int = 1

    # intializing
    x_history = []
    x_history.append(x0)

    f, gradf = func(x_history[-1])
    phi, dphi = f, gradf

    #gradf_history = []
    #gradf_history.append(gradf)
    normgradf_history = []
    
    normgradf = np.linalg.norm(gradf)
    normgradf_history.append(normgradf)
    I = np.eye(len(x0))
    p_k_history = []

    while normgradf_history[-1] > epsilon_g or k > maxk: 

        #print('iteration = ', k, ' error = ', normgradf_history[-1])

        if reset:
            V = I * (1 / np.linalg.norm(dphi))
        else:
            s = np.array([x_history[-1]-x_last])
            y = np.array([dphi - dphi_last])
            sigma_BFGS = 1 / ( 0.001+y@s.T)
            term1 = (I-(sigma_BFGS*((s.T@(y)))))
            term2 = (I-(sigma_BFGS*((y.T@(s)))))
            hess_inter = term1@V@term2
            V = hess_inter + (sigma_BFGS*((s.T@(s)))) # BFGS Update

        #last iteration values
        x_last, dphi_last = x_history[-1], dphi

        p_k = np.dot(-V, dphi)
        p_k_history.append(p_k)

        alpha, phi, dphi = line_search(alpha_int, x_history[-1], mu_1, mu_2, p_k_history[-1], rho, sigma, phi, dphi, func, 1)

        x = x_history[-1] + alpha * p_k_history[-1]
        x_history.append(x)

        normgradf = np.linalg.norm(dphi) #update
        normgradf_history.append(normgradf)

        # reseting the approximate hessian matrix
        if  abs(np.dot(np.transpose(dphi), p_k_history[-1])) <= reset_point:
            reset = True
            #print('Updating')
        else:
            reset = False

        k += 1
    
    xopt = x_history[-1]
    fopt, garbage = func(xopt)

    return xopt, fopt, k, x_history, normgradf_history


