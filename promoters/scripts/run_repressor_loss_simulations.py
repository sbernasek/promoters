import numpy as np
import pickle
from promoters.models.linear import LinearModel
from promoters.simulation.environment import ConditionSimulation


N = 10000

conditions = ('normal', 'diabetic', 'hyper_metabolic',)

# define feedback strengths
eta = (5e-4, 1e-4, 5e-4)

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
        model = LinearModel(g1=0.01, g2=0.001, include_activation=True)
        model.add_feedback(*permanent, perturbed=False)
        model.add_feedback(*removed, perturbed=True)

        # run simulation
        sim = ConditionSimulation(model)
        sim.run(skwargs=dict(N=N, conditions=conditions), ckwargs=dict(mode='empirical'))
        simulations[(i, j)] = sim

# save simulations
with open('pairwise_repressors_sweep.pkl', 'wb') as file:
    pickle.dump(simulations, file)
