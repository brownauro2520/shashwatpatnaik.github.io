import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
import time
from new_uncon_optimizer import uncon_optimizer


#### EXAMPLE 5.4 #########
#exterior function
def f_ext_quad(x):
    x1 = x[0]
    x2 = x[1]
    f = x1+2*x2+mu_pen/2 * (max(0, 1/4 * x1**2 + x2**2 - 1))**2
    g = np.zeros(2)
    g[0] = 1.0 + mu_pen * 1/2 * x1 * (1/4 * x1**2 + x2**2 - 1)
    g[1] = 2.0 + mu_pen * 2 * x2 * (1/4 * x1**2 + x2**2 - 1)
    return f, g

def f_int_log(x):
    x1 = x[0]
    x2 = x[1]
    
    argument = -1/4 * x1**2 - x2**2 + 1
    if argument < 1e-16:
        return 0,np.zeros(2)

    f = x1+2*x2 - mu_pen * np.log(-1/4 * x1**2 - x2**2 + 1)
    g = np.zeros(2)
    g[0] = 1 - mu_pen * (-1/2 * x1) / (-1/4 * x1**2 - x2**2 + 1)
    g[1] = 2 - mu_pen * (-2 * x2) / (-1/4 * x1**2 - x2**2 + 1)
    return f, g



def f_ext_quad_real(x):
    f=x[0]+(2*x[1])
    return f

def ineq(x):
    G = (((x[0]**2)/4) + (x[1]**2) - 1)
    return G


##### BEAM ##########

def Lbeam_int(x):
    # Define constants
    b = 0.125
    h = 0.25
    P = 1E5
    l = 1
    sig = 2E8
    tao = 1.16E8

    # Define symbolic variables
    t_b, t_w = sp.symbols("t_b t_w")

    # Define expressions
    I = (h**3 * t_w) / 12 + (b * t_b**3) / 6 + (h**2 * b * t_b) / 2
    g1 = ((P * l * h) / (2) - sig*I)/sig
    g2 = ((1.5 * P) / (h) - tao*t_w)/tao

    f = 2 * b * t_b + h * t_w
    fhat = f - mu_pen*sp.log(-g1) - mu_pen*sp.log(-g2)
    ghat = np.array([sp.diff(fhat,t_b),sp.diff(fhat,t_w)])
    fhat_ = fhat.subs({t_b: x[0], t_w: x[1]})
    if not isinstance(fhat_, sp.Float): #sometimes breaks when complex
        fhat_ = abs(fhat_)
    fhatval = np.float64(fhat_)
    ghatval = np.array([np.float64(ghat[0].subs({t_b: x[0], t_w: x[1]})),np.float64(ghat[1].subs({t_b: x[0], t_w: x[1]}))])
    return fhatval, ghatval


def Lbeam_ext(x):
    # Define constants
    b = 0.125
    h = 0.25
    P = 1E5
    l = 1
    sig = 2E8
    tao = 1.16E8

    # Define symbolic variables
    t_b, t_w = sp.symbols("t_b t_w")

    # Define expressions
    I = (h**3 * t_w) / 12 + (b * t_b**3) / 6 + (h**2 * b * t_b) / 2
    g1 = (P * l * h) / (2) - sig*I
    g2 = (1.5 * P) / (h) - tao*t_w
    f = 2 * b * t_b + h * t_w
    
    fhat = f + 0.5*mu_pen*(g1**2+g2**2)
    ghat = np.array([sp.diff(fhat,t_b),sp.diff(fhat,t_w)])
    
    fhatval = np.float64(fhat.subs({t_b: x[0], t_w: x[1]}))
    ghatval = np.array([np.float64(ghat[0].subs({t_b: x[0], t_w: x[1]})),np.float64(ghat[1].subs({t_b: x[0], t_w: x[1]}))])
    return fhatval, ghatval

def beam_func(x):
    b = 125e-3
    h = 250e-3
    area = (2*b*x[0])+(h*x[1])
    return area


#maximum iteration
maxk = 29

#int = 0 vs iter, int = 1 vs mu
def plot_error(x_history, mu_history, x0, mu, rho, int):
    # ploting
    iter =np.linspace(1, len(x_history), len(x_history))
    error = np.linalg.norm(np.array(x_history) - np.array(x_true), axis = 1)
    # Create a semilogy plot
    if int == 0:
        plt.semilogy(iter, error)
        plt.xlabel('Iteration', fontname='Times New Roman', fontsize=12)
    else:
        plt.loglog(mu_history, error)
        plt.xlabel('μ', fontname='Times New Roman', fontsize=12)
    plt.ylabel(r'Error $|\widehat{f^*}-f^*|$', fontname='Times New Roman', fontsize=12)
    plt.title('Convergence Plot with starting guess [%.2f, %0.2f]\nPenaly Parameters μ = %.3f and ρ =%0.2f' %(x0[0], x0[1], mu, rho), fontname='Times New Roman', fontsize=12)
    plt.grid(True)
    plt.show()

def plot(x_history, lim, x0, mu, rho, func, func_ineq):
    # Create a grid of x and y values for the contour plot
    x = np.linspace(lim[0], lim[1], 400)
    y = np.linspace(lim[2], lim[3], 400)
    X, Y = np.meshgrid(x, y)

    # Compute the function values at each point in the grid
    Z = func([X, Y])

    # Create the contour plot
    plt.contour(X, Y, Z, levels=20, cmap='viridis')
    plt.colorbar()
    plt.xlabel(r'$x_0$', fontname='Times New Roman', fontsize=12)
    plt.ylabel(r'$x_1$', fontname='Times New Roman', fontsize=12)

    G = func_ineq([X, Y])
    plt.contourf(X, Y, G, [0, 100000], colors='red', alpha=0.2)
    plt.contour(X, Y, G, levels=[0], colors='red', linestyles='dashed')
    # Plot x_history as a line plot
    x_hist = np.array(x_history)
    plt.plot(x_hist[:, 0], x_hist[:, 1], c='k', marker='o', markersize=5, label='Optimization Path')

    # Mark the optimum point with a green star
    optimum_point = x_history[-1]
    plt.scatter(optimum_point[0], optimum_point[1], c='b', marker='*', s=200, label='Optimum Point')

    plt.legend()
    plt.title('Optimization Path with starting guess [%.2f, %0.2f]\nPenaly Parameters μ = %.3f and ρ =%0.2f' %(x0[0], x0[1], mu, rho), fontname='Times New Roman', fontsize=12)
    plt.show()

def plot_beam(x_history, x0, mu, rho):
    # constant
    b = 125e-3 # m 
    h = 250e-3 # m
    P = 100e3 # N
    l = 1 # m
    sigma_yield = 200e6 # pa
    tau_yield = 116e6 # pa

    # Define functions
    def f(tb, tw):
        return (2 * b * tb) + (h * tw)

    def I(tb, tw):
        return ((h**3 * tw) / 12) + ((tb**3 * b) / 6) + ((h**2 * b * tb) / 2)

    def g1(tb, tw):
        return (P * l * h / (2 * I(tb, tw))) - sigma_yield

    def g2(tb, tw):
        return (1.5 * P / (h * tw)) - tau_yield

    # Create a grid of values for tb and tw
    tb_values = np.linspace(0.005, 0.025, 400)
    tw_values = np.linspace(0.0013, 0.02, 400)

    # Create meshgrid
    tb_mesh, tw_mesh = np.meshgrid(tb_values, tw_values)

    # Compute function values
    f_values = f(tb_mesh, tw_mesh)

    # Create contour plot for f
    plt.figure(figsize=(10, 6))

    #plot the end result
    #plt.scatter(0.0142604, 0.00517241, c='k', marker='*', s=200,  label='optimal point')

    plt.contour(tb_mesh, tw_mesh, f_values, levels=30, cmap='viridis')
    plt.colorbar(label=r'f($t_b$, $t_w$)')

    tb=tb_values
    tw_g1=(12/h**3)*(((P*l*h)/(2*sigma_yield))-(b*(tb**3)/6)-((h**2)*b*tb/2))
    plt.plot(tb,tw_g1,'--b',label='g1')

    tw_g2=1.5*P/(h*tau_yield)+(0*tb_values)
    plt.plot(tb_values,tw_g2,'--r',label='g2')


    G1 = g1(tb_mesh, tw_mesh); G2 = g2(tb_mesh, tw_mesh)

    plt.contourf(tb_mesh, tw_mesh, G1, levels=[0, 1e9], colors ='b', alpha = 0.2)
    plt.contourf(tb_mesh, tw_mesh, G2, levels=[0, 1e10], colors ='r', alpha = 0.2)

    plt.xlim([0.005, 0.025])
    plt.ylim([0.0013, 0.02])

    # Plot x_history as a line plot
    x_hist = np.array(x_history)
    plt.plot(x_hist[:, 0], x_hist[:, 1], c='k', marker='o', markersize=5, label='Optimization Path')

    # Mark the optimum point with a green star
    optimum_point = x_history[-1]
    plt.scatter(optimum_point[0], optimum_point[1], c='b', marker='*', s=200, label='Optimum Point')
    plt.xlabel(r'$t_b$')
    plt.ylabel(r'$t_w$')
    plt.title('Optimization Path with starting guess [%.2f, %0.2f]\nPenaly Parameters μ = %.3f and ρ =%0.2f' %(x0[0], x0[1], mu, rho), fontname='Times New Roman', fontsize=12)
    plt.legend()
    plt.show()

# Interior Penalty Method ------------------------------------------------------
def interior_penalty(x0, rho_pen, func, func_real, epsilon_g, options='quasi_bp'):
    global mu_pen
    mu_pen = 0.001
    k = 1
    x_history = []
    x_history.append(x0)
    mu_history = []
    mu_history.append(mu_pen)
    error = 20
    while error > 1e-3:
        
        x, f, data = uncon_optimizer(func, x_history[-1], epsilon_g, options)
        mu_pen = rho_pen*mu_pen
        mu_history.append(mu_pen)
        x_history.append(x)
        error = np.linalg.norm(np.array(x_history[-1]) - np.array(x_true))
        print(k, x_history[-1], error)
        k += 1
    
    xopt = x_history[-1]
    fopt = func_real(xopt)

    return xopt, fopt, x_history, mu_history


# Exterior Penalty Method ------------------------------------------------------
def exterior_penalty(x0, rho_pen, func, func_real, epsilon_g, options='quasi_bp'):
    global mu_pen
    mu_pen = 0.02
    k = 1
    x_history = []
    x_history.append(x0)
    mu_history = []
    mu_history.append(mu_pen)
    error = 20
    while error > 1e-4:
        
        x, f, data = uncon_optimizer(func, x_history[-1], epsilon_g, options)
        mu_pen = rho_pen*mu_pen
        mu_history.append(mu_pen)
        x_history.append(x)
        error = np.linalg.norm(np.array(x_history[-1]) - np.array(x_true))
        #print(k, x_history[-1], error)
        k += 1
    
    xopt = x_history[-1]
    fopt = func_real(xopt)

    return xopt, fopt, x_history, mu_history


#interior
# max mu 18000 with rho 0.7 but iteration taken were more 
#exterior 
# max mu 9000 with rho at 5 it take less iteration as mu increases 
#both cases mu greatly influce the first step it took, 
#in both cases linesearch algorthm failed to find a minmum and and step was taken to wrong place

# penalty paramters 
rho_int_exterior = 5
mu_int = 0.001
rho_int_interior = 0.99







# choose penalty method and function
name = 'e_example'

if name == 'e_example':
    # intial guess for example
    x0 = [1, 0.3]
    x_true = [-1.41421356237, -0.70710678118]
    #print('Exterior Optimal Point true = -1.41421356237, -0.70710678118' )
    x_star, f_star, x_history, mu_history = exterior_penalty(x0, rho_int_exterior, f_ext_quad, f_ext_quad_real, 1e-4)
    plot_error(x_history, mu_history, x0, mu_int, rho_int_exterior, 0)
    plot_error(x_history, mu_history, x0, mu_int, rho_int_exterior, 1)
    plot(x_history, [-3.5, 3, -2, 1.5], x0, mu_int, rho_int_exterior, f_ext_quad_real, ineq)

elif name == 'i_example':
    # intial guess for example
    x0 = [1, 0.3]
    x_true = [-1.41421356237, -0.70710678118]
    print('Interior Optimal Point true = -1.41421356237, -0.70710678118' )
    x_star, f_star, x_history, mu_history = interior_penalty(x0, rho_int_interior, f_int_log, f_ext_quad_real, 1e-4)
    plot_error(x_history, mu_history, x0, mu_int, rho_int_interior, 0)
    plot_error(x_history, mu_history, x0, mu_int, rho_int_interior, 1)
    plot(x_history, [-3.5, 3, -2, 1.5], x0, mu_int, rho_int_interior, f_ext_quad_real, ineq)

elif name == 'e_beam':
    # intial guess for beam
    x0 = [0.02, 0.008]
    x_true = [0.0142603955, 0.00517241379]
    print('Exterior beam Optimal Point true = 0.0142603955, 0.00517241379' )
    x_star, f_star, x_history, mu_history = exterior_penalty(x0, rho_int_exterior, Lbeam_ext, beam_func, 1e-2)
    plot_error(x_history, mu_history, x0, mu_int, rho_int_exterior, 0)
    plot_error(x_history, mu_history, x0, mu_int, rho_int_exterior, 1)
    plot_beam(x_history, x0, mu_int, rho_int_exterior)

elif name == 'i_beam':
    # intial guess for beam
    x0 = [0.02, 0.008]
    x_true = [0.0142603955, 0.00517241379]
    print('Interior beam Optimal Point true = 0.0142603955, 0.00517241379' )
    x_star, f_star, x_history, mu_history = interior_penalty(x0, rho_int_interior, Lbeam_int, beam_func, 1e-2)
    x_history.append([0.01426, 0.00517])
    #plot_error(x_history, mu_history, x0, mu_int, rho_int_interior, 0)
    #plot_error(x_history, mu_history, x0, mu_int, rho_int_interior, 1)
    plot_beam(x_history, x0, mu_int, rho_int_interior)





