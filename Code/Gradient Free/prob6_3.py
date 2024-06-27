import numpy as np
import matplotlib.pyplot as plt
from gradfree import NDM
from funcgrad import *
from scipy.optimize import minimize


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
    


## tolerance paramters
tau_f=10**-6
l=1

### My gradient free
run = 1

dim_run = [2, 4, 6, 8, 10, 12, 16, 32, 40]

# storing history
iter_hist_ndm = []; iter_hist_sndm = []; iter_hist_sbfa = []; iter_hist_sbfd = []
func_hist_ndm = []; func_hist_sndm = []; func_hist_sbfa = []; func_hist_sbfd = []
f_hist_ndm = []; f_hist_sndm = []; f_hist_sbfa = []; f_hist_sbfd = []
x_normerror_ndm = []; x_normerror_sndm = []; x_normerror_sbfa = []; x_normerror_sbfd = []


for dim in dim_run:
    if dim <13:
        tau_x=10**-6
    else:
        tau_x=10**-3

    #print(' ---------------------------------- Rosenbrock Dimenson ---------------------------------- ', dim)
    #correct answer
    x_c = np.ones(dim)

    #initial guess
    x0= np.ones(dim)*2
    func1=func_rosen
    func = FunctionWrapperWithCounter(func1)


    X_opt, F_opt, error_x, error_f, traingle_hist = NDM(x0, l, func, tau_x, tau_f)
    """ print("Nelder-Mead Gradient free Optimizer")
    print("Optimal parameters:", X_opt)
    print("Optimal function value:",F_opt) """
    iterations = list(range(len(error_x)))
    """ print("Number of iteration = ", iterations[-1])
    print('Function evaluations:', func.get_fcalls()) """
    iter_hist_ndm.append(iterations[-1]); func_hist_ndm.append(func.get_fcalls()); f_hist_ndm.append(F_opt)
    norm = np.linalg.norm(X_opt - x_c)
    x_normerror_ndm.append(norm)
    # Set the tolerance
    tolerance = 1e-6

    # Use scipy minimize with Nedler method
    """ print("SCIPY --- Nelder-Mead Gradient free Optimizer") """
    result = minimize(func_rosen_real, x0, method='Nelder-Mead', tol=tolerance, options={'disp': False, 'maxiter': 30000})
    # Print the result
    """ print("Optimal parameters:", result.x)
    print("Optimal function value:", result.fun)
    print("Number of iterations:", result.nit)
    print("Number of function evaluations:", result.nfev) """
    iter_hist_sndm.append(result.nit); func_hist_sndm.append(result.nfev); f_hist_sndm.append(result.fun)
    norm = np.linalg.norm(result.x - x_c)
    x_normerror_sndm.append(norm)

    # Use scipy minimize with BFGS method - analytical derivatives
    """ print("SCIPY --- BFGS Analytical Derivatives") """
    result = minimize(func_rosen_real, x0, method='BFGS', jac=rosen_grad, tol=tolerance, options={'disp': False})
    # Print the result
    """ print("Optimal parameters:", result.x)
    print("Optimal function value:", result.fun)
    print("Number of iterations:", result.nit)
    print("Number of function evaluations:", result.nfev) """
    iter_hist_sbfa.append(result.nit); func_hist_sbfa.append(result.nfev); f_hist_sbfa.append(result.fun)
    norm = np.linalg.norm(result.x - x_c)
    x_normerror_sbfa.append(norm)

    # Use scipy minimize with BFGS method - finite difference
    """ print("SCIPY --- BFGS Finite Difference Derivatives") """
    result = minimize(func_rosen_real, x0, method='BFGS', options={'disp': False})
    # Print the result
    """ print("Optimal parameters:", result.x)
    print("Optimal function value:", result.fun)
    print("Number of iterations:", result.nit)
    print("Number of function evaluations:", result.nfev) """
    iter_hist_sbfd.append(result.nit); func_hist_sbfd.append(result.nfev); f_hist_sbfd.append(result.fun)
    norm = np.linalg.norm(result.x - x_c)
    x_normerror_sbfd.append(norm)
    

# Plotting iteration histories
plt.plot(dim_run, iter_hist_ndm, label='Nelder-Mead algorithm')
plt.plot(dim_run, iter_hist_sndm, label='SciPy - Nelder-Mead')
plt.plot(dim_run, iter_hist_sbfa, label='SciPy - BFGS Analytical Derivatives')
plt.plot(dim_run, iter_hist_sbfd, label='SciPy - BFGS Finite Difference Derivatives')

# Adding labels and legend
plt.xlabel('Dimension')
plt.yscale('log')
plt.ylabel('Number of Iterations')
plt.title('Comparison of Optimization Algorithms')
plt.legend(fontsize=8)

# Show the plot
plt.show()

# Plotting function calls histories
plt.plot(dim_run, func_hist_ndm, label='Nelder-Mead algorithm')
plt.plot(dim_run, func_hist_sndm, label='SciPy - Nelder-Mead')
plt.plot(dim_run, func_hist_sbfa, label='SciPy - BFGS Analytical Derivatives')
plt.plot(dim_run, func_hist_sbfd, label='SciPy - BFGS Finite Difference Derivatives')

# Adding labels and legend
plt.xlabel('Dimension')
plt.yscale('log')
plt.ylabel('Number of Function Calls')
plt.title('Comparison of Optimization Algorithms')
plt.legend(fontsize=8)

# Show the plot
plt.show()

# Plotting optimal function values
plt.plot(dim_run, f_hist_ndm, label='Nelder-Mead algorithm')
plt.plot(dim_run, f_hist_sndm, label='SciPy - Nelder-Mead')
plt.plot(dim_run, f_hist_sbfa, label='SciPy - BFGS Analytical Derivatives')
plt.plot(dim_run, f_hist_sbfd, label='SciPy - BFGS Finite Difference Derivatives')

# Adding labels and legend
plt.xlabel('Dimension')
plt.yscale('log')
plt.ylabel('Optimal Function Value')
plt.title('Comparison of Optimization Algorithms')
plt.legend(fontsize=8)

# Show the plot
plt.show()

# Plotting error norms of optimal design points
plt.plot(dim_run, x_normerror_ndm, label='Nelder-Mead algorithm')
plt.plot(dim_run, x_normerror_sndm, label='SciPy - Nelder-Mead')
plt.plot(dim_run, x_normerror_sbfa, label='SciPy - BFGS Analytical Derivatives')
plt.plot(dim_run, x_normerror_sbfd, label='SciPy - BFGS Finite Difference Derivatives')

# Adding labels and legend
plt.xlabel('Dimension')
plt.ylabel('Error Norm of Optimal Design Point')
plt.yscale('log')
plt.title('Comparison of Optimization Algorithms')
plt.legend()

# Show the plot
plt.show()
