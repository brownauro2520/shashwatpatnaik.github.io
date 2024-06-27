import os
import argparse
import matplotlib.pyplot as plt
from pyoptsparse import History
import numpy as np

plt.rcParams["text.usetex"] = False  # Comment out if latex installation is not present
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.size"] = 20
plt.rcParams["xtick.labelsize"] = 16
plt.rcParams["ytick.labelsize"] = 16
plt.rcParams["lines.linewidth"] = 3.0

parser = argparse.ArgumentParser()
parser.add_argument("--histFile", type=str, default="output_almost there no lete/opt.hst")
parser.add_argument("--outputDir", type=str, default="./")
args = parser.parse_args()

optHist = History(args.histFile)
histValues = optHist.getValues()

val = histValues["DVCon1_thickness_constraints_0"]
iterrr = histValues["iter"]
cl = histValues["fc_cl"]
tiii = histValues["time"]
print(tiii[-1], iterrr[-1])
cl = cl * 2.78
val = np.array(val)
val_store = []
for i in range(len(val)):
    val_store.append(np.sum(val[i]))

fig, axes = plt.subplots(2,figsize=(12, 8))
axes[0].plot(iterrr[:210], cl[:210],  label="Cl")
axes[0].axhline(0.14, linestyle="--", color="gray")
axes[0].axhline(0.2, linestyle="--", color="gray")
axes[0].set_ylabel(r"$C_L$", rotation="horizontal", ha="right", fontsize=14)
axes[0].set_xlabel('Iterations', fontsize=14)
axes[0].autoscale(enable=True, tight=True)

axes[1].plot("iter", "fc_cd", data=histValues, label="Cd")
axes[1].set_ylabel(r"$C_D$", rotation="horizontal", ha="right", fontsize=14)
axes[1].set_xlabel('Iterations', fontsize=14)
axes[1].set_xlim(0,210)





# for ax in axes:
#     ax.spines["right"].set_visible(False)
#     ax.spines["top"].set_visible(False)

#     # Drop the rest of the spines
#     ax.spines["left"].set_position(("outward", 12))
#     ax.spines["bottom"].set_position(("outward", 12))

fig.show()
plt.savefig(os.path.join(args.outputDir, "aero_wing_opt_hist.png"))
