import numpy as np
from os import mkdir
from promoters.models.linear import LinearModel
from promoters.simulation.environment import ConditionSimulation


N = 10000

conditions = ('normal', 'diabetic', 'hyper_metabolic')

# define feedback strengths
eta = (5e-4/2, 1e-4/2, 5e-4/2)

# run pairwise simulations
for i, name in enumerate(['activation','transcription','translation']):

    # make directory
    mkdir(name)

    # define feedback strengths
    removed = np.zeros(3)
    removed[i] = eta[i]

    # define model
    model = LinearModel(g1=0.01, g2=0.001, include_activation=True)
    model.add_feedback(*removed, perturbed=True)

    # run simulation
    sim = ConditionSimulation(model)
    sim.run(skwargs=dict(N=N, conditions=conditions), ckwargs=dict(mode='empirical'))

    # save simulation
    sim.save(name, saveall=True)
