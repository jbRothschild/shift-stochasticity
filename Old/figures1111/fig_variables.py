import PopDyn_variation/shiftStochasticity as ss #All our models
techniques = ['FP QSD', 'FP Gaussian', 'WKB Realspace', 'QSD Algorithm', r"$\tau[1]$", "small n approx.", "Exact Solution"]
colors_techniques = plt.cm.viridis(np.linspace(0.,1.,len(techniques))) #BuPu
lines = [':', '-', ':', '-', ':', '-', '-']
n = 10
colors_gradient = plt.cm.inferno(np.linspace(0,1,n))
colors_gradient2 = plt.cm.YlGn(np.linspace(0,1,n))

capacity = 100.0
stochasticity = np.linspace(0.01, .99, 100)
variability = np.logspace(-1.0, 1.0, 1000)
cap = np.linspace(1.0, capacity, num=capacity)

stochasticity_i = np.linspace(0.1, 0.3, n)
#stochasticity_i = np.linspace(0.1, 0.9, n)
variability_i = np.logspace(-1.0, 0.0, n)

var = 500
sto = 5
K = 99
MAX = ss.maxPop(cap[K], stochasticity[0], variability[-1])
