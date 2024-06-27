import numpy as np
import matplotlib.pyplot as plt
from gradfree import NDM
from funcgrad import *
from scipy.optimize import minimize

## tolerance paramters
tau_x=10**-6
tau_f=10**-6 
l=1

######################## bean function
run = 0

if run == 1:
    x0= np.array([-2, 2])
    func=bean

    X_opt, F_opt, error_x, error_f, traingle_hist = NDM(x0, l, func, tau_x, tau_f)
    print("Nelder-Mead Gradient free Optimizer")
    print("Optimal parameters:", X_opt)
    print("Optimal function value:",F_opt)
    # Set the tolerance
    tolerance = 1e-9
    # Use scipy minimize with BFGS method
    result = minimize(bean_real, x0, method='BFGS', jac=bean_grad, tol=tolerance, options={'disp': True})
    # Print the result
    print("Optimal parameters:", result.x)
    print("Optimal function value:", result.fun)

    plot = 1
    if plot == 1:
        plotc(traingle_hist, X_opt, x0, func, tau_x, tau_f)
        plot_error(error_x, error_f)

##################### bean function with noise
run = 0

if run == 1:
    x0= np.array([-2, 2])
    func=beannoise

    X_opt, F_opt, error_x, error_f, traingle_hist = NDM(x0, l, func, tau_x, tau_f)
    print("Nelder-Mead Gradient free Optimizer")
    print("Optimal parameters:", X_opt)
    F_opt = bean_real(X_opt)
    print("Optimal function value:",F_opt)
    iterations = list(range(len(error_x)))
    print("Number of iteration = ", iterations[-1])
    # Set the tolerance
    tolerance = 1e-9
    # Use scipy minimize with BFGS method
    result = minimize(beannoise_real, x0, method='BFGS', jac=beannoise_grad, tol=tolerance, options={'disp': True})
    # Print the result
    print("Optimal parameters:", result.x)
    F_opt = bean_real(result.x)
    print("Optimal function value:", F_opt)

    plot = 1
    if plot == 1:
        error_f = error_f[:-1:50]
        error_x = error_x[:-1:50]
        #plotc(traingle_hist, X_opt, x0, func, tau_x, tau_f)
        plot_error(error_x, error_f)


################# bean function checkboard
run = 0

if run == 1:
    x0= np.array([1, 1])
    func=beancheckerboard

    X_opt, F_opt, error_x, error_f, traingle_hist = NDM(x0, l, func, tau_x, tau_f)
    print("Nelder-Mead Gradient free Optimizer")
    print("Optimal parameters:", X_opt)
    F_opt = bean_real(X_opt)
    print("Optimal function value:",F_opt)
    iterations = list(range(len(error_x)))
    print("Number of iteration = ", iterations[-1])
    # Set the tolerance
    tolerance = 1e-4
    # Use scipy minimize with BFGS method
    result = minimize(beancheckerboard_real, x0, method='BFGS', jac=beancheckerboard_grad, tol=tolerance, options={'disp': True})
    # Print the result
    print("Optimal parameters:", result.x)
    print("Optimal function value:", result.fun)

    plot = 0
    if plot == 1:
        plotc(traingle_hist, X_opt, x0, func, tau_x, tau_f)
        plot_error(error_x, error_f)


################# problem d) function 
run = 0

if run == 1:
    x0= np.array([3, 2, 1])
    func=three

    X_opt, F_opt, error_x, error_f, traingle_hist = NDM(x0, l, func, tau_x, tau_f)
    print("Nelder-Mead Gradient free Optimizer")
    print("Optimal parameters:", X_opt)
    print("Optimal function value:",F_opt)
    iterations = list(range(len(error_x)))
    print("Number of iteration = ", iterations[-1])
    # Set the tolerance
    tolerance = 1e-12
    # Use scipy minimize with BFGS method
    result = minimize(three_real, x0, method='BFGS',  tol=tolerance, options={'disp': True})
    # Print the result
    print("Optimal parameters:", result.x)
    print("Optimal function value:", result.fun)




