import numpy as np
import sys, os
sys.path.append(os.path.abspath('..'))
import mte.shiftStochasticity as ss
from fig_variables import stochasticity, K, variability, stochasticity_i, variability_i
import fig_variables as fv
import matplotlib.pyplot as plt
plt.style.use('parameters.mplstyle')  # particularIMporting

def make_figure_B1():
    """
    slices of q in mle
    """
    figname = 'varyQ'
    curve1 = DATADICT[figname + '_numeric']
    curve2 = DATADICT[figname + '_heuristic']
    # plot
    #plt.figure(figsize=(10, 5))
    plt.figure()
    c2_part1 = plt.plot(curve2['xpts'][0:400], curve2['ypts'][0:400], color=cs['heuristic'], label='heuristic')
    #c2_part2 = plt.plot(curve2['xpts'][400:], curve2['ypts'][400:], color=cs['heuristic'])
    c1 = plt.plot(curve1['xpts'], curve1['ypts'], marker='o', linestyle='None', color=cs['numerical_fisher_sp'], label='numeric')
    plt.title('Mode 1 MLE: Numeric vs Heuristic ($k_p=10$, $t=100$, $k_{off}=1$)')
    plt.xlabel(r'$q$')
    plt.ylabel(r'$mean time  to extinction$')
    plt.legend()
    # save figure
    plt.gca().set_ylim([-5, max(curve1['ypts'])])
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.pdf', transparent=True)
    plt.savefig(DIR_OUTPUT + os.sep + figname + '.eps')
