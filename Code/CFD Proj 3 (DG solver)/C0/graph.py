import numpy as np
import matplotlib.pyplot as plt

CFL = [0.6,0.45,0.2,0.1]
for i in range(4):
    R = np.loadtxt("R0P" + str(i) + ".txt")
    y = np.linspace(1, R.size, num=R.size)
    # Create a semilogy plot
    plt.semilogy(y, R,  label="P" + str(i) + " CFL " + str(CFL[i]))

# Add labels and title
plt.xlabel('Number of Iteration')
plt.ylabel(r'|$L_1$|')
plt.title(r'$|L_1|_{norm}$ convergence of Airfoil for Coarsest Mesh', fontsize = 16)

# Show the plot
plt.legend()
plt.show()