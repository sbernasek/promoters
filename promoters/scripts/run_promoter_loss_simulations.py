import numpy as np
from os import mkdir
from promoters.models.linear import LinearModel
from promoters.simulation.environment import ConditionSimulation


N = 10000

conditions = ('normal', 'diabetic', 'carbon_limited', 'hyper_metabolic')

# define promoter strengths
eta = (0.5,0.5,0.5)

# run pairwise simulations
for i, name in enumerate(['activation','transcription','translation']):

    # make directory
    mkdir(name)

    # define feedback strengths
    removed = np.zeros(3)
    removed[i] = eta[i]

    # define model
    model = LinearModel(g1=0.01, g2=0.001, include_activation=True)
    model.add_promoters(*removed, perturbed=True)

    # run simulation
    sim = ConditionSimulation(model)
    sim.run(skwargs=dict(N=N, conditions=conditions))

    # save simulation
    sim.save(name, saveall=True)
