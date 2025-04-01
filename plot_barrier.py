# plot the barrier certificates obtained from running the Julia script. Indices of state variables corrected by 1

import numpy as np
import itertools
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle

savetype = "png" # "svg"
def Bar2D(x): # from harder condition, 2D
    return(0.21028*x[0]**2 - 0.00375*x[0]*x[1] + 0.20497*x[1]**2 + 0.75363*x[0] - 1.19902*x[1] - 0.52708)


def plot2d(): # plot barrier above
    # initial, unsafe sets and state space
    # First 2 for x1, last 2 for x2
    Xo = [0.9,1.1, 0.9,1.1] # initial set
    Xu = [1.5,2, 0,0.5] # unsafe set
    X = [-4,2, 0,5] # state space

    n = 100 # number of sampled points
    x1,x2 = np.linspace(X[0],X[1],n+1), np.linspace(X[2],X[3],n+1)
    X1,X2 = np.meshgrid(x1,x2)
    fig, ax = plt.subplots()
    # filled contour for 0-sublevel set of barriers (B(x) <= 0)
    ax.contourf(x1,x2, Bar2D([X1,X2]), levels = 0, colors=['orange','w'],alpha=0.25) 
    ax.add_patch(Rectangle((3,2), 0, 0, color = 'orange',alpha=0.3,label=r"$B(x)\leq 0$")) # added for legend, no rectangle plotted

    # plot initial and unsafe regions
    ax.add_patch(Rectangle((Xo[0], Xo[2]), Xo[1]-Xo[0], Xo[3]-Xo[2],color = '#006400',alpha=0.7,label=r"$Initial$"))
    ax.add_patch(Rectangle((Xu[0], Xu[2]), Xu[1]-Xu[0], Xu[3]-Xu[2],color = 'red',alpha=0.9,label=r"$Unsafe$"))
    ax.axis('scaled')
    ax.grid()
    ax.legend(loc="upper left")
    plt.xlabel(fr"$x_{1}$") # f - evaluate {}, r - latex style
    plt.ylabel(fr"$x_{2}$")
    plt.savefig(f"./media/hard_system_2d.{savetype}", bbox_inches='tight', format=savetype, dpi=1200)
    plt.show()


def Bar5D(x):  # from relaxed condition, 5D
    return(1.0e-5*x[3]**5 - 1.0e-5*x[0]**4 - 1.0e-5*x[0]**2*x[3]**2 - 1.0e-5*x[0]*x[1]*x[3]**2 - 1.0e-5*x[0]*x[2]*x[3]**2 - 2.0e-5*x[0]*x[3]**3 - 1.0e-5*x[0]*x[3]**2*x[4] - 2.0e-5*x[1]**4 - 1.0e-5*x[1]**2*x[2]**2 - 4.0e-5*x[1]**2*x[3]**2 - 1.0e-5*x[1]**2*x[4]**2 + 3.0e-5*x[1]*x[3]**3 - 2.0e-5*x[2]**4 - 4.0e-5*x[2]**2*x[3]**2 - 1.0e-5*x[2]**2*x[4]**2 + 3.0e-5*x[2]*x[3]**3 - 0.00025*x[3]**4 + 3.0e-5*x[3]**3*x[4] - 3.0e-5*x[3]**2*x[4]**2 - 1.0e-5*x[4]**4 - 1.0e-5*x[0]**3 + 1.0e-5*x[0]**2*x[2] + 3.0e-5*x[0]**2*x[3] + 2.0e-5*x[0]*x[1]**2 + 3.0e-5*x[0]*x[2]**2 + 0.00015*x[0]*x[3]**2 + 2.0e-5*x[0]*x[4]**2 + 0.00012*x[1]**3 + 5.0e-5*x[1]**2*x[3] - 1.0e-5*x[1]**2*x[4] - 1.0e-5*x[1]*x[2]*x[3] - 3.0e-5*x[1]*x[3]**2 - 1.0e-5*x[1]*x[3]*x[4] + 0.00012*x[2]**3 + 5.0e-5*x[2]**2*x[3] - 1.0e-5*x[2]**2*x[4] - 3.0e-5*x[2]*x[3]**2 - 1.0e-5*x[2]*x[3]*x[4] + 0.00112*x[3]**3 - 5.0e-5*x[3]**2*x[4] + 3.0e-5*x[3]*x[4]**2 + 7.0e-5*x[4]**3 + 5.0e-5*x[0]**2 - 3.0e-5*x[0]*x[1] - 4.0e-5*x[0]*x[2] - 9.0e-5*x[0]*x[3] - 3.0e-5*x[0]*x[4] - 0.0003*x[1]**2 + 2.0e-5*x[1]*x[2] - 0.0001*x[1]*x[3] + 2.0e-5*x[1]*x[4] - 0.00028*x[2]**2 - 8.0e-5*x[2]*x[3] + 2.0e-5*x[2]*x[4] - 0.00114*x[3]**2 - 5.0e-5*x[3]*x[4] - 0.00021*x[4]**2 - 8.0e-5*x[0] + 5.0e-5*x[1] + 4.0e-5*x[2] + 0.0002*x[3] + 5.0e-5*x[4] - 0.00026)


def plot2dproj(): # plot barrier above as 2D projection of select indices
    # initial, unsafe sets and state space
    choose2 = [[1,4],[3,4],[5,4]]
    for i in range(len(choose2)):
        nonzero_indices = list(choose2[i])
        nonzero_indices = [x - 1 for x in nonzero_indices]
        # zero_indices = list(set(all_indices).difference(nonzero_indices))
        Xomin = [0.9, 0.9, 0.9, 0.9, 0.9] # initial set (x1,x2,x3,x4,x5).
        Xomax = [1.1, 1.1, 1.1, 1.1, 1.1]
        Xumin = [1.5, 0, 0, 2.5, 0] # unsafe set
        Xumax = [2, 0.5, 0.5, 3, 0.5]
        Xmin = [-4, 0, 0, 0, 0] # state space
        Xmax = [2, 5, 5, 3, 6]
        Xo = [Xomin[nonzero_indices[0]],Xomax[nonzero_indices[0]], Xomin[nonzero_indices[1]],Xomax[nonzero_indices[1]]]
        Xu = [Xumin[nonzero_indices[0]],Xumax[nonzero_indices[0]], Xumin[nonzero_indices[1]],Xumax[nonzero_indices[1]]]
        X = [Xmin[nonzero_indices[0]],Xmax[nonzero_indices[0]], Xmin[nonzero_indices[1]],Xmax[nonzero_indices[1]]]

        n = 100 # number of sampled points
        x1,x2 = np.linspace(X[0],X[1],n+1), np.linspace(X[2],X[3],n+1)
        X1,X2 = np.meshgrid(x1,x2)
        fig, ax = plt.subplots()

        # filled contour for 0-sublevel set of barriers (B(x) <= 0)
        nz = np.zeros((n+1,n+1))
        x = [nz,nz,nz,nz,nz]
        x[nonzero_indices[0]] = X1
        x[nonzero_indices[1]] = X2
        ax.contourf(x1,x2, Bar5D(x), levels = 0, colors=['orange','w'],alpha=0.25)
        ax.add_patch(Rectangle((0,0), 0, 0, color = 'orange',alpha=0.4,label=r"$B(x)\leq 0$")) # added for legend, not plotted

        # plot initial and unsafe regions
        ax.add_patch(Rectangle((Xo[0], Xo[2]), Xo[1]-Xo[0], Xo[3]-Xo[2],color = '#006400',alpha=0.7,label=r"$Initial$"))
        ax.add_patch(Rectangle((Xu[0], Xu[2]), Xu[1]-Xu[0], Xu[3]-Xu[2],color = 'red',alpha=0.9,label=r"$Unsafe$"))
        ax.axis('scaled')
        ax.grid()
        ax.legend(loc="best")
        plt.xlabel(fr"$x_{nonzero_indices[0]+1}$")
        plt.ylabel(fr"$x_{nonzero_indices[1]+1}$")
        plt.savefig(f"./media/projected_system_{nonzero_indices[0]+1}_{nonzero_indices[1]+1}.{savetype}", bbox_inches='tight', format=savetype, dpi=1200)
        plt.show()

plot2d() # plot barrier for 2D system
plot2dproj() # plot barrier for 5D system as 2D projection