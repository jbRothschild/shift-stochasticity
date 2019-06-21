import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import sys, os
sys.path.append(os.path.abspath('..')) #To add the model modules paths
import mte.shiftStochasticity as ss #All our models

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns #niceeee colors

DIR_INPUT = "input"
DIR_OUTPUT = "output"

for dirs in [DIR_INPUT, DIR_OUTPUT]:
    if not os.path.exists(dirs):
        os.makedirs(dirs)

techniques = ['FP QSD', 'FP Gaussian', 'WKB Realspace', 'QSD Algorithm', r"$\tau[1]$", "small n approx.", "Exact Solution"]
colors_techniques = plt.cm.viridis(np.linspace(0.,1.,len(techniques))) #BuPu
lines = ['-',':', '-', ':', '-', ':', '-', ':']
n = 10
colors_gradient = plt.cm.inferno(np.linspace(0,1,n))
colors_gradient2 = plt.cm.YlGn(np.linspace(0,1,n))

capacity = 100.0
stoc = 100
vari = 1001
#stochasticity = np.linspace(0.01, .99, 100)
stochasticity = np.linspace(0.01, 0.99, stoc) #100
variability = np.logspace(-1.0, 1.0, vari) #1001
cap = np.linspace(1.0, capacity, num=capacity)
K = 99

sto = 20
#sto = 70
var = 300
#var = 800


stochasticity_i = np.linspace(0.1, 0.3, n)
variability_i = np.logspace(-1.0, 0.0, n)
print(stochasticity_i)
print(variability_i)


MAX = ss.maxPop(cap[K], stochasticity[0], variability[-1])

# plot params
FS = 16
SHOW = True

# axes
POINTS_BETWEEN_TICKS_X = stoc/5
POINTS_BETWEEN_TICKS_Y = (vari-1)/2

def plot_heatmap(arr, xrange, yrange, fname, label, show=SHOW):
    # TODO change colour scheme, see https://matplotlib.org/examples/color/colormaps_reference.html
    # TODO fix ticks randomly disappearing on colourbar + flip colourbar minor ticks or remove?
    """
    Colours viridis, YlGnBu, terrain, plasma, BuPu
    """
    nlevels = 50
    LOG = False
    # plot setup
    f = plt.figure()
    if LOG:
        imshow_kw = { 'cmap': 'YlGnBu', 'aspect': None, 'vmin': np.min(arr), 'vmax': np.max(arr), 'norm': mpl.colors.LogNorm(nlevels)}
        im = plt.contourf(arr, levels=np.logspace(np.log10(np.min(arr)),np.log10(np.max(arr)),nlevels), **imshow_kw) #doesn't give values of colorbar
    else:
        imshow_kw = { 'cmap': 'inferno', 'aspect': None, 'vmin': np.min(arr), 'vmax': np.max(arr)}
        im = plt.contourf(arr, nlevels, **imshow_kw)


    # axes setup
    ax = plt.gca()
    for c in ax.collections:
        c.set_edgecolor("face")
    # method 1
    ax.set_xticks([i for i, xval in enumerate(xrange) if i % POINTS_BETWEEN_TICKS_X == 0])
    ax.set_yticks([i for i, yval in enumerate(yrange) if i % POINTS_BETWEEN_TICKS_Y == 0])
    # values to have on the axis
    #ax.set_xticklabels([r'$10^{%d}$' % np.log10(xval) for i, cval in enumerate(xrange) if i % POINTS_BETWEEN_TICKS_X==0], fontsize=FS)
    ax.set_xticklabels([r'$%.1f$' % round(xval,1) for i, xval in enumerate(xrange) if i % POINTS_BETWEEN_TICKS_X == 0], fontsize=FS)
    ax.set_yticklabels([r'$10^{%d}$' % np.log10(yval) for i, yval in enumerate(yrange) if i % POINTS_BETWEEN_TICKS_Y == 0], fontsize=FS)

    ax.set_xlabel(r'$q$', fontsize=FS)
    ax.set_ylabel(r'$\delta$', fontsize=FS)

    # create colorbar
    if LOG:
        cbar = ax.figure.colorbar(im, ax=ax, ticks=mpl.ticker.LogLocator(base=10))
    else:
        cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel(label, rotation=-90, va="bottom", fontsize=FS, labelpad=20)
    cbar.ax.tick_params(labelsize=FS)
    # TODO ID why do ticks hide sometimes?
    #for t in cbar.ax.get_yticklabels(): print(t.get_text())

    # contour line for value 1.0
    #plt.contour(arr, levels=[1.0], linestyles=['dashed'])  # use 'dashed' or 'solid' curve
    plt.tight_layout()
    # save
    plt.savefig(DIR_OUTPUT + os.sep + fname + '.pdf')
    plt.savefig(DIR_OUTPUT + os.sep + fname + '.eps')
    if show:
        plt.show()
    return
