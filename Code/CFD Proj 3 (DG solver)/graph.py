import numpy as np
import matplotlib.pyplot as plt


for i in range(4):
    R = np.loadtxt("R0P" + str(i) + ".txt")
    y = np.linspace(1, R.size, num=R.size)
    # Create a semilogy plot
    plt.semilogy(y, R,  label="P" + str(i))

# Add labels and title
plt.xlabel('Number of Iteration')
plt.ylabel(r'|$L_1$|')
plt.title('Preservation test of Bump with CFL 0.3 for Mesh0')

# Show the plot
plt.legend()
plt.show()