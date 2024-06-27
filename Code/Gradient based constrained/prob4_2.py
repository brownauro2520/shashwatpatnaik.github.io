import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

# constant
b = 125e-3 # m 
h = 250e-3 # m
P = 100e3 # N
l = 1 # m
sigma_yield = 200e6 # pa
tau_yield = 116e6 # pa

# Define the symbolic variables
t_b, t_w, sigma_1, sigma_2, s_1, s_2 = sp.symbols('t_b t_w sigma_1 sigma_2 s_1 s_2')

# Lagrangian equation
I = ((h**3) / 12)*t_w + (b / 6)*(t_b**3) + ((h**2)*b*0.5)*t_b
lagrangian = 2*b*t_b + h*t_w + sigma_1*(((P*l*h) / (2* I)) - sigma_yield + s_1**2) + sigma_2*(((1.5*P) / (h*t_w)) - tau_yield + s_2**2)

# Differentiate the equation with respect to each variable
dtb = sp.diff(lagrangian, t_b)
dtw = sp.diff(lagrangian, t_w)
dsigma1 = sp.diff(lagrangian, sigma_1)
dsigma2 = sp.diff(lagrangian, sigma_2)
ds1 = sp.diff(lagrangian, s_1)
ds2 = sp.diff(lagrangian, s_2)

con = 1

if con == 1:
    #print('both constraint active')

    # Substitute values into the derivatives
    dtb = dtb.subs({s_1: 0, s_2: 0})
    dtw = dtw.subs({s_1: 0, s_2: 0})
    dsigma1 = dsigma1.subs({s_1: 0, s_2: 0})
    dsigma2 = dsigma2.subs({s_1: 0, s_2: 0})
    

    dtb_both = P*h*l*sigma_1*(-b*h**2 - b*t_b**2)/(b*h**2*t_b + b*t_b**3/3 + h**3*t_w/6)**2 + 2*b
    dtw_both = -P*h**4*l*sigma_1/(6*(b*h**2*t_b + b*t_b**3/3 + h**3*t_w/6)**2) - 1.5*P*sigma_2/(h*t_w**2) + h
    dsigma1_both = P*h*l/(b*h**2*t_b + b*t_b**3/3 + h**3*t_w/6) - sigma_yield
    dsigma2_both = 1.5*P/(h*t_w) - tau_yield

    # Solve the system of equations
    solutions = sp.solve([dtb_both, dtw_both, dsigma1_both, dsigma2_both], (t_b, t_w, sigma_1, sigma_2))

    # Print the solutions
    #print(solutions)

elif con == 2:
    print('neither constraint active')

    # Substitute values into the derivatives
    dtb1 = dtb.subs({sigma_1: 0, sigma_2: 0})
    dtw1 = dtw.subs({sigma_1: 0, sigma_2: 0})
    dsigma11 = dsigma1.subs({sigma_1: 0, sigma_2: 0})
    dsigma21 = dsigma2.subs({sigma_1: 0, sigma_2: 0})


    dtb_both = 0.250000000000000
    dtw_both = 0.250000000000000
    dsigma1_both = s_1**2 - 200000000.0 + 25000.0/(0.0416666666666667*t_b**3 + 0.0078125*t_b + 0.00260416666666667*t_w)
    dsigma2_both = s_2**2 - 116000000.0 + 600000.0/t_w

    # Solve the system of equations
    solutions = sp.solve([dtb_both, dtw_both, dsigma1_both, dsigma2_both], (t_b, t_w, s_1, s_2))
    # Print the solutions
    print(solutions)

elif con == 3:
    print('constarint 1 is active')

    # Substitute values into the derivatives
    dtb1 = dtb.subs({s_1: 0, sigma_2: 0})
    dtw1 = dtw.subs({s_1: 0, sigma_2: 0})
    dsigma11 = dsigma1.subs({s_1: 0, sigma_2: 0})
    dsigma21 = dsigma2.subs({s_1: 0, sigma_2: 0})

    dtb_both = 14400000.0*sigma_1*(-0.125*t_b**2 - 0.0078125)/(t_b**3 + 0.1875*t_b + 0.0625*t_w)**2 + 0.25
    dtw_both = -37500.0*sigma_1/(t_b**3 + 0.1875*t_b + 0.0625*t_w)**2 + 0.25
    dsigma1_both = -200000000.0 + 25000.0/(0.0416666666666667*t_b**3 + 0.0078125*t_b + 0.00260416666666667*t_w)
    dsigma2_both = s_2**2 - 116000000.0 + 600000.0/t_w

    # Solve the system of equations
    solutions = sp.solve([dtb_both, dtw_both, dsigma1_both, dsigma2_both], (t_b, t_w, s_2, sigma_1))
    # Print the solutions
    print(solutions)

else:
    print('constarint 2 is active')

    # Substitute values into the derivatives
    dtb1 = dtb.subs({s_2: 0, sigma_1: 0})
    dtw1 = dtw.subs({s_2: 0, sigma_1: 0})
    dsigma11 = dsigma1.subs({s_2: 0, sigma_1: 0})
    dsigma21 = dsigma2.subs({s_2: 0, sigma_1: 0})

    dtb_both = 0.250000000000000
    dtw_both = -600000.0*sigma_2/t_w**2 + 0.25
    dsigma1_both = s_1**2 - 200000000.0 + 25000.0/(0.0416666666666667*t_b**3 + 0.0078125*t_b + 0.00260416666666667*t_w)
    dsigma2_both = -116000000.0 + 600000.0/t_w

    # Solve the system of equations
    solutions = sp.solve([dtb_both, dtw_both, dsigma1_both, dsigma2_both], (t_b, t_w, s_1, sigma_2))
    # Print the solutions
    print(solutions)



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
tb_values = np.linspace(0.005, 0.1, 400)
tw_values = np.linspace(0.0013, 0.02, 400)

# Create meshgrid
tb_mesh, tw_mesh = np.meshgrid(tb_values, tw_values)

# Compute function values
f_values = f(tb_mesh, tw_mesh)

# Create contour plot for f
plt.figure(figsize=(10, 6))

#plot the end result
plt.scatter(0.0142604, 0.00517241, c='k', marker='*', s=200,  label='optimal point')

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

plt.xlim([0.005, 0.1])
plt.ylim([0.0013, 0.02])


plt.xlabel(r'$t_b$')
plt.ylabel(r'$t_w$')
plt.title('Contour plot of the function f with constrained g1 and g2')
plt.legend()
plt.show()

