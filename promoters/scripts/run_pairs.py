import numpy as np
import pickle
from promoters.models.linear import LinearModel
from promoters.simulation.environment import ConditionSimulation


N = 5000

conditions = ('normal', 'diabetic', 'minute', 'carbon_limited', 'hyper_metabolic')

# define promoter strengths
eta = (1,1,1)

# run pairwise simulations
simulations = {}
for i in range(3):
    for j in range(3):

        # define feedback strengths
        permanent = np.zeros(3)
        permanent[i] = eta[i]
        removed = np.zeros(3)
        removed[j] = eta[j]

        # define model
        model = LinearModel(g1=0.01, g2=0.001)
        model.add_promoters(*permanent)
        model.add_promoters(*removed, perturbed=True)

        # run simulation
        sim = ConditionSimulation(model)
        sim.run(skwargs=dict(N=N, conditions=conditions))
        simulations[(i, j)] = sim

# save simulations
with open('pairwise_sweep.pkl', 'wb') as file:
    pickle.dump(simulations, file)
