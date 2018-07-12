import numpy as np

def simple_propensity(parameters, population):
    [cap, delta, q] = parameters
    birth = (1+delta/2)*population-q*population**2/cap
    death = delta*population/2+(1-q)*population**2/cap
    return np.array([birth, death])

def pop_update():
    return np.array([[1], [-1]])

def initial_population():
    return np.array([10])
