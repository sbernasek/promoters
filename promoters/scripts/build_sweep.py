from promoters.execution.arguments import SweepArguments
from promoters.sweep.sweep import LinearSweep, HillSweep, TwoStateSweep, SimpleSweep
from promoters.sweep.dense import SimpleDense2D, SimpleDependenceSweep
from promoters.sweep.influence import RepressorInfluenceSweep, PromoterInfluenceSweep


# ======================== PARSE SCRIPT ARGUMENTS =============================

args = SweepArguments(description='Parameter sweep arguments.')
model = args['model']
num_samples = args['number_of_samples']

# ============================= RUN SCRIPT ====================================

# instantiate sweep object
if model == 'simple':
    sweep = SimpleSweep(num_samples=num_samples)
elif model == 'simple_dependence':
    sweep = SimpleDependenceSweep(num_samples=num_samples)
elif model == 'simple_2d':
    sweep = SimpleDense2D(num_samples=num_samples)
elif model == 'repressor_influence':
    sweep = RepressorInfluenceSweep(num_samples=num_samples)
elif model == 'promoter_influence':
    sweep = PromoterInfluenceSweep(num_samples=num_samples)
elif model == 'linear':
    sweep = LinearSweep(num_samples=num_samples)
elif model == 'hill':
    sweep = HillSweep(num_samples=num_samples)
elif model == 'twostate':
    sweep = TwoStateSweep(num_samples=num_samples)
else:
    raise ValueError('{:s} model type not recognized.'.format(model))

# build sweep
sweep.build(
    directory=args['path'],
    batch_size=args['batch_size'],
    num_trajectories=args['number_of_trajectories'],
    saveall=args['save_all'],
    deviations=args['use_deviations'],
    comparison=args['comparison_mode'],
    walltime=args['walltime'],
    allocation=args['allocation'],
    pulse_start=args['pulse_start'],
    pulse_duration=args['pulse_duration'],
    pulse_baseline=args['pulse_baseline'],
    pulse_magnitude=args['pulse_magnitude'],
    pulse_sensitive=args['pulse_sensitive'],
    simulation_duration=args['simulation_duration'])


"""

# 22 / 05 / 10

# LINEAR REPRESSORS
python build_sweep.py -N 5000 -n 1000 -b 25 -s 0 -w 10 -pb 0.0 -ps 50 -sd 250 -d false -cm empirical -m linear

# LINEAR
python build_sweep.py -N 5000 -n 1000 -b 25 -s 0 -w 10 -pb 0.0 -ps 50 -sd 250 -d false -cm promoters -m linear

# TWOSTATE
python build_sweep.py -N 5000 -n 1000 -b 25 -s 0 -w 10 -pb 0.0 -ps 50 -sd 250 -d false -cm promoters -m twostate

# HILL
python build_sweep.py -N 5000 -n 1000 -b 25 -s 0 -w 10 -pb 0.0 -ps 50 -sd 250 -d false -cm promoters -m hill

# BASAL STIMULUS
python build_sweep.py -N 5000 -n 1000 -b 25 -s 0 -w 10 -pb 0.1 -ps 50 -sd 250 -d true -cm promoters -m linear

# PULSE SENSITIVE TO METABOLISM
python build_sweep.py -N 5000 -n 1000 -b 25 -s 0 -w 10 -pb 0.0 -ps 50 -sd 250 -d false -cm promoters -m linear -PS true

"""
