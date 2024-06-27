import numpy as np
import matplotlib.pyplot as plt

exact_cl = -1.921902839
exact_cd = 0.0000000039
clm0 = np.array([-1.9288588900, -1.9221299800, -1.9219362200, -1.9219325800])
clm1 = np.array([-1.9263735000,	-1.9221131100,	-1.9219378900,	-1.9219169600])
clm2 = np.array([-1.9242972100,	-1.92201157702,	-1.9219145009,	-1.921910244])


cdm0 = abs(np.array([0.0027363100,	0.0000736992,	0.0000043399,	0.0000004707]))
cdm1 = abs(np.array([0.0017357190,	0.0000378307,	0.0000028307,	0.0000003027]))
cdm2 = abs(np.array([0.0005462059,	0.0000187114,	0.0000007675,	0.0000000939]))

clp0 = np.array([-1.9288588900, -1.9263735000, -1.9242972100])
clp1 = np.array([-1.9221299800,-1.9221131100, -1.92201157702])
clp2 = np.array([-1.9219362200, -1.9219378900, -1.9219145009])
clp3 = np.array([-1.9219325800, -1.9219169600, -1.921910244])
p = [0,1,2,3]
#l1_norms = abs(exact_cd-cdm2)
#mesh_sizes = np.array([1, 2,3,4])
mesh_sizes = np.array([102,408,1632])

for i in range(4):
    var_name = 'clp{}'.format(i)
    l1_norms = abs(exact_cl - globals()[var_name])
    # Calculate the convergence rates
    experimental_rates = np.log(l1_norms[:-1] / l1_norms[1:]) / np.log(mesh_sizes[:-1] / mesh_sizes[1:])
    plt.semilogy(mesh_sizes, l1_norms, "o-", label="Rate p{} = {:.2f}".format(i, abs(experimental_rates[1])))
    #experimental_rates = np.log(l1_norms[:-1] / l1_norms[1:]) / np.log(mesh_sizes[:-1] / mesh_sizes[1:])
    print(experimental_rates)

# Plot the convergence study


plt.xlabel(r"Number of elements")
plt.ylabel("L1 Norm")
plt.title(r"Convergence Study for $F_y$")
plt.legend()
plt.show()
