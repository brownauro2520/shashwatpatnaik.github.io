import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy.optimize import minimize

def ones_vector(length, value):
    # Initialize an empty list
    x = []
    # Use a loop to add 1 to the list 'length' times
    for i in range(length):
        x.append(1*value)

    return x

class FunctionWrapperWithCounter():
    """
    Attach a function call counter to a user-defined function

    To use:
    --------------------------------------------------
    # instantiate
    func = FunctionWrapperWithCounter(user_defined_function)

    # call function
    func(x)

    # get function calls
    fcalls = func.get_fcalls()
    --------------------------------------------------
    """

    def __init__(self, func):
        self._func = func
        self._fcalls = 0   # initialize function call counter
        self._max_fcalls = 1E7   # maximum iterations

    def __call__(self, x):
        """ calls the user-provided function and count the function calls"""
        self._fcalls += 1

        # kill and exit if the function calls reached the limit
        if self._fcalls >= self._max_fcalls:
            raise RuntimeError('You reached the function call limit of ' + str(self._max_fcalls) + '. Exit.')
        
        return self._func(x)

    def get_fcalls(self):
        """ get function calls"""
        return self._fcalls

#================= slant ===================
def func_slanted_quad(x):
    beta = 1.5

    f = x[0]**2 + x[1]**2 - beta * x[0] * x[1]

    g = np.zeros(2)
    g[0] = 2 * x[0] - beta * x[1]
    g[1] = 2 * x[1] - beta * x[0]
    return f, g

#================= Rosenbrock ===================
def func_rosen(x):
    x= np.array(x)
    """
    Rosenbrock function

    Parameters
    ----------
    x : ndarray, shape (dim,)
        design variables
    dim : int, optional
        Dimension of the Rosenbrock function, default is 2

    Returns
    -------
    f : float
        function value
    g : ndarray, shape (dim,)
        objective gradient
    """
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

#================= Bean ===================
def func_bean(x):
    """
    Bean function from Appendix D of the book

    Parameters
    ----------
    x : ndarray, shape (2,)
        design variables
    
    Returns
    -------
    f : float
        function value
    g : ndarray, shape (2,)
        objective gradient
    """

    x1 = x[0]
    x2 = x[1]

    f = (1. - x1)**2 + (1. - x2)**2 + 0.5 * (2 * x2 - x1**2)**2

    g = np.zeros(2)
    g[0] = -2. * (1 - x1) - 2. * x1 * (2 * x2 - x1**2)
    g[1] = -2. * (1 - x2) + 2. * (2 * x2 - x1**2)

    return f, g

def uncon_optimizer(func, x0, epsilon_g, options=None):
    output = {}

    if options is None:
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        linesearch = backtrack
        namelinesearch="bactrack"
        searchdir = BFGS
        namesearchdir="quasi"

    # Direction Search Algorithm 
    if options == 's_bt':
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir = stepest
        namesearchdir="steepdes"
        linesearch = backtrack
        namelinesearch="backtracking"

    elif options == 'c_bt':
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir=conjgrad  
        namesearchdir="conjgrad"
        linesearch = backtrack
        namelinesearch="backtracking"

    elif options == 'quasi_bt':
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir = BFGS
        namesearchdir="quasi"
        linesearch = backtrack
        namelinesearch="backtracking"
    
    elif options == 's_bp':
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir = stepest
        namesearchdir="steepdes"
        linesearch = bracketting
        namelinesearch="bracketing"
    
    elif options == 'c_bp':
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir=conjgrad  
        namesearchdir="conjgrad"
        linesearch = bracketting
        namelinesearch="bracketing"

    elif options == 'quasi_bp':
        mu_1 = 1e-4; mu_2 = 0.9; rho = 0.3; sigma = 2; reset_point = 1e-8
        searchdir = BFGS  
        namesearchdir="quasi"
        linesearch = bracketting
        namelinesearch="bracketing"
        
    elif options == 'scipy_BFGS':
        # Use BFGS method to find the minimum
        error_history = []
        initial_error = func(x0)[0]
        error_history.append(initial_error)
        
        # Define the callback function to store errors
        def callback(x):
            error = func(x)[0]  # Get the error at the current iteration
            error_history.append(error)
        
        result = minimize(lambda x: func(x)[0], x0, method='BFGS', jac=lambda x: func(x)[1], callback=callback)
        iterations = result.nit
        xopt = result.x
        fopt = result.fun



    if options == 'scipy_BFGS':
        output = {'alias':"browncurry2520", 'iteration':iterations, 'error history': error_history}
        return xopt, fopt, output
    else:
        xopt, fopt, k, x_his, error_his = searchdir(x0, epsilon_g, mu_1, mu_2, rho, sigma, reset_point, func, linesearch)
        output = {'alias':"browncurry2520", 'iteration':k,'x history': x_his, 'error history': error_his, 'namesearchdir':namesearchdir,'namelinesearch':namelinesearch}
        return xopt, fopt, output

# backtracking algorithm 
def backtrack(alpha_int, x0, mu_1, mu_2, p_k, rho, sigma, phi_0, dphi_0, func, inf):
    # variables to store data
    alpha_history = []
    # # backtracking algorithm begins
    alpha = alpha_int
    alpha_history.append(alpha)
    phi, dphi = func(np.array(x0) + alpha*np.array(p_k))
    k = 0
    while phi > phi_0 + mu_1*alpha*np.dot(dphi_0,p_k):
        #print(phi_0, k)
        alpha = alpha * rho
        alpha_history.append(alpha)
        phi, dphi = func(np.array(x0) + alpha*np.array(p_k))
        k+=1
    return alpha_history[-1], phi, dphi

#quadratic interpolation
def alpha_interpolation2(alpha_1, alpha_2, phi_1, phi_2, dphi_1, dphi_2):
    numerator = (2*alpha_1*(phi_2 - phi_2)) + (dphi_1*(alpha_1**2 - alpha_2**2))
    denominator = 2*(phi_2 - phi_1 + dphi_1*(alpha_1 - alpha_2))
    alpha_star_inter = numerator / denominator
    return alpha_star_inter

#cubic interpolation
def alpha_interpolation3(alpha_1, alpha_2, phi_1, dphi_1, phi_2, dphi_2):
    beta_1 = dphi_1 + dphi_2 - 3 * ((phi_1 - phi_2) / (alpha_1 - alpha_2))
    beta_2 = np.sign(alpha_2 - alpha_1) * np.sqrt(beta_1**2 - (dphi_1*dphi_2))
    alpha_star_inter = alpha_2 - (alpha_2 - alpha_1) * ((dphi_2 + beta_2 - beta_1) / (dphi_2 - dphi_1 + 2*beta_2))
    if np.abs(alpha_1-alpha_2) < 1E-9: # switch to bisection
        return 0.5*(alpha_1+alpha_2)
    return alpha_star_inter

# Pinpointing Function
def pinpoint(a_low, a_high, phi_0, dphi_0,  x0, p_k, mu1, mu2, func):#
    while True:
        phi_low, dphi_low = func(np.array(x0) + a_low*np.array(p_k))
        phi_high, dphi_high = func(np.array(x0) + a_high*np.array(p_k))
        alpha_p = alpha_interpolation3(a_low, a_high, phi_low, np.dot(dphi_low, p_k), phi_high, np.dot(dphi_high, p_k)) 
        phi_p, dphi_p = func(np.array(x0) + alpha_p*np.array(p_k))


        if phi_p > phi_0 + mu1*alpha_p*np.dot(dphi_0, p_k) or phi_p > phi_low:
            a_high = alpha_p
            phi_high = phi_p
            dphi_high = dphi_p
        else:
            if abs(np.dot(dphi_p, p_k)) <= -mu2*np.dot(dphi_0, p_k):
                return alpha_p, phi_p, dphi_p
            elif np.dot(dphi_p, p_k)*(a_high-a_low) >= 0:
                a_high = a_low
            a_low = alpha_p

# Bracketting Function
def bracketting(alpha, x0, mu_1, mu_2, p_k, rho, sigma, phi_0, dphi_0, func, inf):
    
    # variables to store data
    alpha_1 = 0
    alpha_2 = alpha
    phi_1 = phi_0
    dphi_1 = dphi_0
    
    first = True
    
    while True:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
        phi_2, dphi_2 = func(np.array(x0) + alpha_2*np.array(p_k))
        if (phi_2 > phi_0 + mu_1*alpha_2*np.dot(dphi_0, p_k)) or (not first and phi_2 > phi_1):
            a_star, phi_star, dphi_star = pinpoint(alpha_1, alpha_2, phi_0, dphi_0, x0, p_k, mu_1, mu_2, func)
            break
        
        if abs(np.dot(dphi_2, p_k)) <= -mu_2*np.dot(dphi_0, p_k):
            a_star = alpha_2
            phi_star, dphi_star = phi_2, dphi_2 
            break
        elif np.dot(dphi_2, p_k) >= 0:
            a_star,  phi_star, dphi_star = pinpoint(alpha_1, alpha_2, phi_0, dphi_0, x0, p_k, mu_1, mu_2, func)
            break
        else:
            alpha_1 = alpha_2
            alpha_2 = sigma*alpha_2

        first = False
    
    return a_star, phi_star, dphi_star

# BFGS or Quasi-Newton Direction
def BFGS(x0, epsilon_g, mu_1, mu_2, rho, sigma, reset_point, func, line_search):
    #maximum iteration
    maxk = 30000
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
            sigma_BFGS = 1 / ( y@s.T)
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

# stepest descent Direction
def stepest(x0, epsilon_g, mu_1, mu_2, rho, sigma, reset_point, func, line_search):
    #maximum iteration
    maxk = 30000
    k = 0
    alpha_int = 1.2
    # intializing
    x_history = []
    x_history.append(x0)

    f, gradf = func(x_history[-1])
    phi, dphi = f, gradf
    normgradf = np.linalg.norm(gradf)
    normgradf_history = []
    normgradf_history.append(normgradf)
    p_k_history = []
    alpha_history = []
    alpha_history.append(alpha_int)

    while normgradf_history[-1] > epsilon_g and k <= maxk: 

        #print('iteration stepest = ', k, ' error = ', normgradf_history[-1])
        p_k = - np.array(dphi) / np.array(normgradf_history[-1])
        p_k_history.append(p_k)

        #guessing alpha from previous step
        if k == 0:
            alpha_int = 1.2
            dphi_last = dphi
        else:
            #alpha_int = 1.2
            alpha_int = alpha_history[-1] * (np.dot(np.array(dphi_last).T, np.array(p_k_history[-2])) / np.dot(np.array(dphi).T, np.array(p_k_history[-1])))
            dphi_last = dphi
        alpha_star, phi, dphi = line_search(alpha_int, x_history[-1], mu_1, mu_2, p_k_history[-1], rho, sigma, phi, dphi, func, 1)
        alpha_history.append(alpha_star)

        x = x_history[-1] + alpha_star * p_k_history[-1]
        x_history.append(x)

        normgradf = np.linalg.norm(dphi) #update
        normgradf_history.append(normgradf)

        k += 1

    # computing optimal function value
    xopt=x_history[-1]
    fopt=func(x_history[-1])[0]

    return xopt, fopt, k, x_history, normgradf_history

# conjugate gradient
def conjgrad(x0, epsilon_g, mu_1, mu_2, rho, sigma, reset_point, func, line_search): #input x0 = starting point, t=convergence rate, and func (as a pure function)
     #maximum iteration
    maxk = 30000
    k = 0
    alpha_int = 1

    # intializing
    x_history = []
    x_history.append(x0)

    f, gradf = func(x_history[-1])
    phi, dphi = f, gradf
    normgradf_history = []
    normgradf = np.linalg.norm(gradf)
    normgradf_history.append(normgradf)
    p_k_history = []

    while normgradf > epsilon_g and k <= maxk:

        if k==0:
            p_k = - np.array(dphi) / np.array(normgradf_history[-1])
            p_k_history.append(p_k)
            dphi_last = dphi
            
        else:
            B = np.dot(np.array(dphi).T, np.array(dphi)) / np.dot(np.array(dphi_last).T, np.array(dphi_last))
            p_k = - np.array(dphi) / np.array(normgradf_history[-1]) + (np.dot(B, p_k_history[-1]))
            p_k_history.append(p_k)
            dphi_last = dphi
        
        alpha_star, phi, dphi = line_search(alpha_int, x_history[-1], mu_1, mu_2, p_k_history[-1], rho, sigma, phi, dphi, func, 1)

        x = x_history[-1] + alpha_star * p_k_history[-1]
        x_history.append(x)

        normgradf = np.linalg.norm(dphi) #update
        normgradf_history.append(normgradf)
        
        k += 1

    # computing optimal function value
    xopt=x_history[-1]
    fopt=func(x_history[-1])[0]

    return xopt, fopt, k, x_history, normgradf_history


if __name__ == "__main__":
        
    epsilon_g=10**-6 #should be as close to zero as possible
    maxj=10 #limit the lopp of pinpoint
    maxk=300000 #limit the loop of search direction

    x0 = ones_vector(400, 0) #starting guess


    # choose function
    func1 = func_slanted_quad
    func_name = 'Slanted Quadratic Function'
    func1 = func_rosen
    func_name = '2D Rosenbrock Function'

    func = FunctionWrapperWithCounter(func1)

    # choose optimization method
    options = 'quasi_bp'
    
    xopt, fopt, outputs = uncon_optimizer(func, x0, epsilon_g, options) 
    #print('Scipy Minimize - BFGS')
    print(outputs["namesearchdir"], outputs["namelinesearch"])
    print('Optimization terminated successfully.')
    print('         alias is: ',outputs['alias'])
    #print('         Optimal value:',xopt)
    print('         Current function value:',fopt)
    print('         Iterations:',outputs['iteration'])
    print('         Function evaluations:', func.get_fcalls())
       