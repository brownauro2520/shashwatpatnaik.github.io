import numpy as np
import matplotlib.pyplot as plt

##############################


# problem 2 part d function
def three(x):
    f=abs(x[0])+2*abs(x[1])+x[2]**2
    dfdx1=1
    dfdx2=2
    dfdx3=2*x[2]
    gradf=[dfdx1, dfdx2, dfdx3]
    return f, gradf

def three_grad(x):
    dfdx1=1
    dfdx2=2
    dfdx3=2*x[2]
    gradf=[dfdx1, dfdx2, dfdx3]
    return gradf

def three_real(x):
    f=abs(x[0])+2*abs(x[1])+x[2]**2
    return f


# bean function
def bean(x):
    f=((1-x[0])**2)+((1-x[1])**2)+(0.5*((2*x[1])-x[0]**2)**2)
    dfdx1=-2.0*x[0]*(-x[0]**2 + 2*x[1]) + 2*x[0] - 2
    dfdx2=-2.0*x[0]**2 + 6.0*x[1] - 2
    gradf=[dfdx1, dfdx2] #grad of function (as a vector)
    return f, gradf

def bean_grad(x):
    dfdx1=-2.0*x[0]*(-x[0]**2 + 2*x[1]) + 2*x[0] - 2
    dfdx2=-2.0*x[0]**2 + 6.0*x[1] - 2
    gradf=[dfdx1, dfdx2] #grad of function (as a vector)
    return gradf

def bean_real(x):
    f=((1-x[0])**2)+((1-x[1])**2)+(0.5*((2*x[1])-x[0]**2)**2)
    return f


# bean function with noise
global magnitude
magnitude = 1e-5

def beannoise(x):
    # Generate random noise using a Gaussian distribution
    noise = np.random.normal(loc=0, scale=magnitude)

    # Add noise to the function value
    f=((1-x[0])**2)+((1-x[1])**2)+(0.5*((2*x[1])-x[0]**2)**2) + noise
    dfdx1=-2.0*x[0]*(-x[0]**2 + 2*x[1]) + 2*x[0] - 2 + noise
    dfdx2=-2.0*x[0]**2 + 6.0*x[1] - 2 + noise
    gradf=[dfdx1, dfdx2] 

    return f, gradf 

def beannoise_grad(x):
    # Generate random noise using a Gaussian distribution
    noise = np.random.normal(loc=0, scale=magnitude)
    # Add noise to the function value
    dfdx1=-2.0*x[0]*(-x[0]**2 + 2*x[1]) + 2*x[0] - 2 + noise
    dfdx2=-2.0*x[0]**2 + 6.0*x[1] - 2 + noise
    gradf=[dfdx1, dfdx2] 
    return gradf

def beannoise_real(x):
    # Generate random noise using a Gaussian distribution
    noise = np.random.normal(loc=0, scale=magnitude)
    # Add noise to the function value
    f=((1-x[0])**2)+((1-x[1])**2)+(0.5*((2*x[1])-x[0]**2)**2) + noise
    return f


# bean function checkboard
global mag
mag = 0.5
def beancheckerboard(x):

    f=((1-x[0])**2)+((1-x[1])**2)+(0.5*((2*x[1])-x[0]**2)**2)+(mag*np.ceil((np.sin(np.pi*x[0])*np.sin(np.pi*x[1]))))
    dfdx1=-2.0*x[0]*(-x[0]**2 + 2*x[1]) + 2*x[0] - 2 +(mag*np.pi*np.sin(np.pi*x[1])*np.cos(np.pi*x[0]))
    dfdx2=-2.0*x[0]**2 + 6.0*x[1] - 2 +(mag*np.pi*np.sin(np.pi*x[0])*np.cos(np.pi*x[1]))
    gradf=[dfdx1, dfdx2] #grad of function (as a vector)
    return f, gradf 

def beancheckerboard_grad(x):

    dfdx1=-2.0*x[0]*(-x[0]**2 + 2*x[1]) + 2*x[0] - 2 +(mag*np.pi*np.sin(np.pi*x[1])*np.cos(np.pi*x[0]))
    dfdx2=-2.0*x[0]**2 + 6.0*x[1] - 2 +(mag*np.pi*np.sin(np.pi*x[0])*np.cos(np.pi*x[1]))
    gradf=[dfdx1, dfdx2] #grad of function (as a vector)
    return gradf 

def beancheckerboard_real(x):

    f=((1-x[0])**2)+((1-x[1])**2)+(0.5*((2*x[1])-x[0]**2)**2)+(mag*np.ceil((np.sin(np.pi*x[0])*np.sin(np.pi*x[1]))))
    return f


# Rosenbrock Function
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


# Ploting functions
def plotc(traingle_hist, X_opt, x0, func, tau_x, tau_f):
    # Plotting the triangles
    simplex_cord = []
    for i in range(len(traingle_hist)):
        triangle_coords = traingle_hist[i]
        for j in range(len(triangle_coords)):

            simplex_cord.append([triangle_coords[j][0], triangle_coords[j][1]])

    simplex_cord = np.array(simplex_cord)
    # Scatter plot of the final optimal point
    plt.scatter(X_opt[0], X_opt[1], color='r', marker='*', s=100, label='Optimal Point', zorder=2)
    plt.scatter(x0[0], x0[1], color='b', marker='D', s=50, label='Starting Point', zorder=3)
    plt.plot(simplex_cord[:, 0], simplex_cord[:, 1], 'o-k', alpha=0.6, label = 'Optimization Path', zorder=1)


    delta=0.01    
        
    #range for bean
    x = np.arange(-4, 4, delta)
    y = np.arange(-2, 4, delta)   

    X, Y = np.meshgrid(x,y)  
    Z = func([X,Y])[0]    

    plt.contour(X, Y, Z, cmap='viridis', levels=15)  
    plt.colorbar() 

    # Setting plot labels and legend
    plt.xlabel('X', fontname='Times New Roman', fontsize=12)
    plt.ylabel('Y', fontname='Times New Roman', fontsize=12)
    plt.title(r'Nelder–Mead algorithm Optimization Path' + '\n' + r'with tolerance $𝜏_x$ = {:.0e} and $𝜏_f$ = {:.0e}'.format(tau_x, tau_f), fontname='Times New Roman', fontsize=12)
    plt.legend()

    # Show the plot
    plt.show()


def plot_error(error_x, error_f):
    iterations = list(range(len(error_x)))
    # Plotting error_x and error_f against the iteration number
    plt.plot(iterations, error_x, '.-r', label=r'$𝜏_x$')
    plt.plot(iterations, error_f, '.-b',label=r'$𝜏_f$')
    plt.yscale('log')
    # Setting plot labels and legend
    plt.xlabel('Iteration', fontname='Times New Roman', fontsize=12)
    plt.ylabel('Error', fontname='Times New Roman', fontsize=12)
    plt.title('Error Convergence Plot' + '\n' + 'Iteration taken to converge = '+ str(iterations[-1]), fontname='Times New Roman', fontsize=12)
    plt.legend()

    # Show the plot
    plt.show()