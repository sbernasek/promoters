from time import time
from promoters.simulation.environment import ConditionSimulation
from promoters.execution.arguments import RunArguments


# ======================== PARSE SCRIPT ARGUMENTS =============================

args = RunArguments(description='Simulation arguments.')
skwargs = dict(N=args['number_of_trajectories'],
               debug=args['debug'],
               conditions=['normal', 'half_growth'])

ckwargs = dict(horizon=args['horizon'],
               deviations=args['use_deviations'],
               mode=args['comparison_mode'])


# ============================= RUN SCRIPT ====================================

start_time = time()

# run each simulation in batch file
with open(args['path'], 'r') as batch_file:

     # run each simulation
     for path in batch_file.readlines():

          # load simulation
          simulation = ConditionSimulation.load(path.strip())

          # run simulation and comparison
          simulation.run(skwargs=skwargs, ckwargs=ckwargs)

          # save simulation
          simulation.save(path.strip(), saveall=args['save_all'])

# print runtime to standard out
runtime = time() - start_time
print('\nSIMULATION COMPLETED IN {:0.2f}.\n'.format(runtime))
