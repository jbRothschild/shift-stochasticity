import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import gillespie_params as gp

def sample_discrete(probs):
    """
    Randomly sample an index with probability given by probs.
    """
    q = np.random.rand() #random number
    i = 0
    p_sum = 0.0
    while p_sum < q: #find index
        p_sum += probs[i]
        i += 1
    return i - 1

def gillespie_draw(propensity_function, params, population):
    # Compute propensities
    propensity = propensity_function(params,population)
    prob = propensity/sum(propensity)

    # Compute time
    time = np.random.exponential(1.0 / sum(propensity))

    # Draw reaction from this distribution
    reaction_index = sample_discrete(prob)

    return reaction_index, time

def gillespie(params, propensity_function, update, times, initial_pop):
    # Initialize output
    update_strat = update()
    output = np.empty((len(times), update_strat.shape[1]), dtype=np.int)

    # Initialize and perform simulation
    i_time = 1
    i = 0
    t = times[0]
    population = initial_pop.copy()
    output[0,:] = population
    n = 0
    while i < len(times):
        while t < times[i_time]:
            # draw the event and time step
            event, dt = gillespie_draw(propensity_function, params, population)

            # Update the population
            population_previous = population.copy()
            population += update_strat[event,:]

            # Increment time
            t += dt

        # Update the index
        i = np.searchsorted(times > t, True)

        # Update the population
        output[i_time:min(i,len(times))] = population_previous

        # Increment index
        i_time = i
    return output

def plot_pop(time_points ,pops, n_simulations):
    # Set up subplots
    fig, ax = plt.subplots(1, pops.shape[2], figsize=(14, 5))

    # Plot state trajectories
    for i in range(n_simulations):
        ax.plot(time_points, pops[i,:,0], '-', lw=0.3, alpha=0.2, color=sns.color_palette()[0])
    # Plot state mean
    ax.plot(time_points, pops[:,:,0].mean(axis=0), '-', lw=6, color=sns.color_palette()[2])
    # Label axes
    ax.set_xlabel('dimensionless time')
    ax.set_ylabel('probability of being found in state')

    plt.tight_layout()
    plt.show()


def main(params, population, number_simulations):
    initial_pop = population
    time = np.linspace(0,400,801)
    np.random.seed(42)
    num_sim = number_simulations
    # Initialize output array
    pops = np.empty((num_sim, len(time), len(initial_pop)))

    # Run the calculations
    for i in range(num_sim):
        pops[i,:,:] = gillespie(params, gp.simple_propensity, gp.pop_update, time, initial_pop)

    plot_pop(time, pops, num_sim)

main([10, 0.5, 0.5], gp.initial_population(), 10)
