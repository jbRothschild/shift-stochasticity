import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os, sys
import mte.shiftStochasticity as ss
import mte.compare_mean_var as comp
#import figures.fig_variables as fv

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


from matplotlib.colors import LogNorm, Normalize #color in log the maps
from matplotlib.ticker import LogLocator #log axes

DIR_OUTPUT = 'figures'

# Fig 1

def plot_qsd(figure='Fig1', nbr_lines=10):
    stochasticity = np.linspace(0.01, 0.99, 100) #100
    variability = np.logspace(-1.0, 1.0, 1001) #1001
    cap = np.linspace(1.0, 100, num=100)
    sto = 20
    var = 300
    K = 99

    bbox = (0.65, 0.95, 0.3, 0.03)
    
    cmap_name = 'viridis'
    cmap = mpl.cm.get_cmap(cmap_name)
    cmap_name2 = 'cividis'

    stochasticity_values = np.linspace(0.2, 0.8, nbr_lines)
    variability_values = np.logspace(-0.6, 0.6, nbr_lines)
    normStoch = Normalize(0.0, 1.0)
    normVaria = LogNorm(0.1, 10.0)
    
    # colors = [cmap(nbr) for nbr in np.linspace(0.0, 0.8, num=len(s))]
    
    fig, ax1 = plt.subplots()

    pdf = []
    
    # A : Varying distribution - delta
    for stochas in stochasticity_values:
        pdf = ss.statDistributionAlgo(ss.maxPop(cap[K],
                                                stochas,
                                                variability[var]),
                                      stochas,
                                      cap[K],
                                      variability[var],
                                      10)
        ax1.plot(range(1, len(pdf)+1), pdf, color=cmap(normStoch(stochas)))

    ax1.set_xlabel("Abundance")
    ax1.set_ylabel("Probability")
    ax1.set_yscale("log")
    ax1.set_ylim(10**(-10), 10**(0))
    ax1.get_yaxis().set_major_locator(LogLocator(numticks=5))
    ax1.minorticks_off()
    ax1.set_xlim(1, 200)
    
    cbaxes = inset_axes(ax1, width='100%', height='100%', bbox_to_anchor=bbox,
                        bbox_transform=ax1.transAxes,
                        loc='center') 
    #cbaxes = inset_axes(ax1, width="30%", height="3%", loc=cloc, bbox_to_anchor=bbox)  
    fig.colorbar(mpl.cm.ScalarMappable(norm=normStoch, cmap=cmap),
                 cax=cbaxes, orientation='horizontal', label=r'$q$')
    
    figname = figure + 'A'
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    plt.close()
    
    # B : Varying distribution - q

    cmap = mpl.cm.get_cmap(cmap_name2)
    
    fig, ax2 = plt.subplots()

    for varia in variability_values:
        pdf = ss.statDistributionAlgo(ss.maxPop(cap[K],
                                                stochasticity[sto],
                                                varia),
                                      stochasticity[sto],
                                      cap[K],
                                      varia,
                                      10)
        ax2.plot(range(1, len(pdf) + 1), pdf, color=cmap(normVaria(varia)))

    ax2.set_xlabel("Population")
    ax2.set_ylabel("Probability density function")
    ax2.set_yscale("log")
    ax2.set_ylim(10**(-10), 10**(0))
    ax2.get_yaxis().set_major_locator(LogLocator(numticks=5))
    ax2.minorticks_off()
    ax2.set_xlim(1, 200)
    
    cbaxes = inset_axes(ax2, width='100%', height='100%', bbox_to_anchor=bbox,
                        bbox_transform=ax2.transAxes,
                        loc='center') 
    fig.colorbar(mpl.cm.ScalarMappable(norm=normVaria, cmap=cmap),
                 cax=cbaxes, orientation='horizontal', label=r'$\delta$')
    

    figname = figure + 'B'
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    plt.close()

    # C : Heatmap of coefficient of variation
    
    fig, ax3 = plt.subplots()
    
    mean, var = comp.mean_var(comp.exact_solution, stochasticity, variability, K + 1)
    var[var < 0] = 0
    
    CoV = np.sqrt(var) / mean
    
    print(np.shape(CoV), np.shape(mean), np.shape(var))
    
    cax = ax3.contourf(stochasticity, variability, var, 100, cmap=plt.cm.magma)#, norm=LogNorm(), levels=np.logspace(minimum, maximum, maximum))
    cbar = fig.colorbar(cax, ticks=[np.min(var), (np.max(var)-np.min(var))/2., np.max(var)])
    for c in ax3.collections:
        c.set_edgecolor("face")
    cbar.ax.set_ylabel(r'Population variance')
    ax3.set_xlabel(r"Quadratic prefactor, $q$")
    ax3.set_ylabel(r"Linear prefactor, $\delta$")
    ax3.set_yscale("log")
    ax3.minorticks_off()
    
    figname = figure + 'C'
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
    plt.close()
    
    
    



# Fig 2

# A : Heatmap

# B : Slice of MTE - delta

# C : Slice of MTE - q


# Fig 3

# A : example of different 

# B : heatmaps (small?) of KL divergence between algorithm and different approximations



# Fig 4

# A : carrying capacity

# B : heatmaps of Ratio


if __name__ == "__main__":
    plot_qsd('Fig1', 10)
