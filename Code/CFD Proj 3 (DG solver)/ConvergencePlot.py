import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

def residual_convergence(file1, file2, file3,fname):
    # Read in the data from the text files
    case1 = np.loadtxt(file1, dtype=float)


    # Set the default font to Times New Roman using LaTeX interpreter
    rc('font',**{'family':'serif','serif':['Times New Roman']})

    # Create a log-log plot
    plt.loglog(case1, label='p = 0')

    # Add axis labels and a title
    plt.xlabel(r'Iteration', fontsize=14)
    plt.ylabel(r'L$_1$ Residual Norm', fontsize=14)
    plt.title(f'DG Solver Convergence for {fname} Mesh', fontsize=16)

    # Calculate the slope of each line
    x1 = np.arange(1, len(case1) + 1)

    slope1, _ = np.polyfit(x1, np.log(case1), 1)


    print('p = 0', slope1)


    plt.legend()
    plt.savefig(f'{fname}_Convergence.png')
    plt.show()

if __name__ == '__main__':
    residual_convergence('R0S2_8.txt', 'R0S2_9.txt', 'R0S2_10.txt','Test')
