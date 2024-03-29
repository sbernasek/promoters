from os.path import join, abspath, relpath, isdir
from os import mkdir, chmod, pardir
import shutil
from glob import glob
import pickle
import numpy as np
from time import time
from datetime import datetime

from ..simulation.environment import ConditionSimulation


class Batch:
    """
    Class defines a collection of batch job submissions for Quest.

    Attributes:

        path (str) - path to batch directory

        script_name (str) - name of script for running batch

        parameters (iterable) - parameter sets

        simulation_paths (dict) - relative paths to simulation directories

        sim_kw (dict) - keyword arguments for simulation

    Properties:

        N (int) - number of samples in parameter space

    """
    def __init__(self, parameters):
        """
        Instantiate batch of jobs.

        Args:

            parameters (iterable) - each entry is a parameter set that defines a simulation. Parameter sets are passed to the build_model method.

        """
        self.simulation_paths = {}
        self.parameters = parameters
        self.script_name = 'run_batch.py'

    def __getitem__(self, index):
        """ Returns simulation instance. """
        return self.load_simulation(index)

    def __iter__(self):
        """ Iterate over serialized simulations. """
        self.count = 0
        return self

    def __next__(self):
        """ Returns next simulation instance. """
        if self.count < len(self.simulation_paths):
            simulation = self.load_simulation(self.count)
            self.count += 1
            return simulation
        else:
            raise StopIteration

    @property
    def N(self):
        """ Number of samples in parameter space. """
        return len(self.parameters)

    @staticmethod
    def load(path):
        """ Load batch from target <path>. """
        with open(join(path, 'batch.pkl'), 'rb') as file:
            batch = pickle.load(file)
        batch.path = path
        return batch

    @staticmethod
    def build_run_script(path,
                         script_name,
                        num_trajectories,
                        saveall,
                        deviations,
                        comparison):
        """
        Writes bash run script for local use.

        Args:

            path (str) - path to simulation top directory

            script_name (str) - name of run script

            num_trajectories (int) - number of simulation trajectories

            saveall (bool) - if True, save simulation trajectories

            deviations (bool) - if True, use deviation variables

            comparison (str) - type of comparison

        """

        # define paths
        path = abspath(path)
        job_script_path = join(path, 'scripts', 'run.sh')

        # copy run script to scripts directory
        run_script = abspath(__file__).rsplit('/', maxsplit=2)[0]
        run_script = join(run_script, 'scripts', script_name)
        shutil.copy(run_script, join(path, 'scripts'))

        # declare outer script that reads PATH from file
        job_script = open(job_script_path, 'w')
        job_script.write('#!/bin/bash\n')

        # move to batch directory
        job_script.write('cd {:s} \n\n'.format(path))

        # run each batch
        job_script.write('echo "Starting all batches at `date`"\n')
        job_script.write('while read P; do\n')
        job_script.write('echo "Processing batch ${P}"\n')
        job_script.write('python ./scripts/{:s}'.format(script_name)+' ${P} ')
        args = (num_trajectories, saveall, deviations, comparison)
        job_script.write('-N {:d} -s {:d} -d {:d} -cm {:s} \n'.format(*args))
        job_script.write('done < ./batches/index.txt \n')
        job_script.write('echo "Completed all batches at `date`"\n')
        job_script.write('exit\n')

        # close the file
        job_script.close()

        # change the permissions
        chmod(job_script_path, 0o755)

    @staticmethod
    def build_submission_script(path,
                                script_name,
                                num_trajectories=5000,
                                saveall=False,
                                deviations=False,
                                comparison='promoters',
                                walltime=10,
                                allocation='p30653'):
        """
        Writes job submission script for QUEST.

        Args:

            path (str) - path to simulation top directory

            script_name (str) - name of run script

            num_trajectories (int) - number of simulation trajectories

            saveall (bool) - if True, save simulation trajectories

            deviations (bool) - if True, use deviation variables

            comparison (str) - type of comparison

            walltime (int) - estimated job run time

            allocation (str) - project allocation, e.g. p30653 (comp. bio)

        """

        # define paths
        path = abspath(path)
        job_script_path = join(path, 'scripts', 'submit.sh')

        # copy run script to scripts directory
        run_script = abspath(__file__).rsplit('/', maxsplit=2)[0]
        run_script = join(run_script, 'scripts', script_name)
        shutil.copy(run_script, join(path, 'scripts'))

        # determine queue
        if walltime <= 4:
            queue = 'short'
        elif walltime <= 48:
            queue = 'normal'
        else:
            queue = 'long'

        # declare outer script that reads PATH from file
        job_script = open(job_script_path, 'w')
        job_script.write('#!/bin/bash\n')

        # move to batch directory
        job_script.write('cd {:s} \n\n'.format(path))

        # begin outer script for processing batch
        job_script.write('while IFS=$\'\\t\' read P\n')
        job_script.write('do\n')
        job_script.write('b_id=$(echo $(basename ${P}) | cut -f 1 -d \'.\')\n')
        job_script.write('   JOB=`sbatch << EOJ\n')

        # =========== begin submission script for individual batch ============
        job_script.write('#! /bin/bash\n')
        job_script.write('#SBATCH -A {:s} \n'.format(allocation))
        job_script.write('#SBATCH --partition {:s} \n'.format(queue))
        job_script.write('#SBATCH --time={0:02d}:00:00 \n'.format(walltime))
        ##job_script.write('#SBATCH -M sebastian@u.northwestern.edu \n')
        job_script.write('#SBATCH --output ./log/${b_id}/outlog \n')
        job_script.write('#SBATCH --error ./log/${b_id}/errlog \n')
        job_script.write('#SBATCH -J ${b_id} \n')
#        job_script.write('#SBATCH --nodes=1:ppn=1 \n')
        job_script.write('#SBATCH --mem=1gb \n\n')

        # load python module and metabolism virtual environment
        job_script.write('module load python/anaconda3\n')
        job_script.write('source activate promoters_env\n\n')

        # move to batch directory
        job_script.write('cd {:s} \n\n'.format(path))

        # run script
        job_script.write('python ./scripts/{:s}'.format(script_name)+' ${P} ')
        args = (num_trajectories, saveall, deviations, comparison)
        job_script.write('-N {:d} -s {:d} -d {:d} -cm {:s} \n'.format(*args))
        job_script.write('EOJ\n')
        job_script.write('`\n\n')
        # ============= end submission script for individual batch ============

        # print job id
        #job_script.write('echo "JobID = ${JOB} submitted on `date`"\n')
        job_script.write('done < ./batches/index.txt \n')
        job_script.write('echo "All batches submitted as of `date`"\n')
        job_script.write('exit\n')

        # close the file
        job_script.close()

        # change the permissions
        chmod(job_script_path, 0o755)

    def build_batches(self, batch_size=25):
        """
        Creates directory and writes simulation paths for each batch.

        Args:

            batch_size (int) - number of simulations per batch

        """

        # get directories for all batches and logs
        batches_dir = join(self.path, 'batches')
        logs_dir = join(self.path, 'log')

        # create index file for batches
        index_path = join(batches_dir, 'index.txt')
        index = open(index_path, 'w')

        # write file containing simulation paths for each batch
        for i, simulation_path in self.simulation_paths.items():

            # determine batch ID
            batch_id = i // batch_size

            # process new batch
            if i % batch_size == 0:

                # open batch file and append to index
                batch_path = join(batches_dir, '{:d}.txt'.format(batch_id))
                index.write('{:s}\n'.format(relpath(batch_path, self.path)))
                batch_file = open(batch_path, 'w')

                # create log directory for batch
                mkdir(join(logs_dir, '{:d}'.format(batch_id)))

            # write paths to batch file
            batch_file.write('{:s}\n'.format(simulation_path))

            # close batch file
            if i % batch_size == (batch_size - 1):
                batch_file.close()
                chmod(batch_path, 0o755)

        index.close()

        chmod(index_path, 0o755)

    def make_directory(self, directory='./'):
        """
        Create directory for batch of jobs.

        Args:

            directory (str) - destination path

        """

        # assign name to batch
        timestamp = datetime.fromtimestamp(time()).strftime('%y%m%d_%H%M%S')
        name = '{:s}_{:s}'.format(self.__class__.__name__, timestamp)

        # create directory (overwrite existing one)
        path = join(directory, name)
        if not isdir(path):
            mkdir(path)
        self.path = path

        # make subdirectories for simulations and scripts
        mkdir(join(path, 'scripts'))
        mkdir(join(path, 'simulations'))
        mkdir(join(path, 'batches'))
        mkdir(join(path, 'log'))

    def build(self,
              directory='./',
              batch_size=25,
              num_trajectories=5000,
              saveall=False,
              deviations=False,
              comparison='empirical',
              walltime=10,
              allocation='b1022',
              **sim_kw):
        """
        Build directory tree for a batch of jobs. Instantiates and saves a simulation instance for each parameter set, then generates a single shell script to submit each simulation as a separate job.

        Args:

            directory (str) - destination path

            batch_size (int) - number of simulations per batch

            num_trajectories (int) - number of simulation trajectories

            saveall (bool) - if True, save simulation trajectories

            deviations (bool) - if True, use deviation variables

            comparison (str) - type of comparison

            walltime (int) - estimated job run time

            allocation (str) - project allocation

            sim_kw (dict) - keyword arguments for ConditionSimulation

        """

        # create batch directory
        self.make_directory(directory)

        # store parameters (e.g. pulse conditions)
        self.sim_kw = sim_kw
        self.batch_size = batch_size

        # determine the simulation type
        if comparison == 'promoters':
            simulation_type = 'promoters'
        else:
            simulation_type = 'repressors'

        # build simulations
        for i, parameters in enumerate(self.parameters):
            simulation_path = join(self.path, 'simulations', '{:d}'.format(i))
            self.simulation_paths[i] = relpath(simulation_path, self.path)
            self.build_simulation(parameters, simulation_type, simulation_path, **sim_kw)

        # save serialized batch
        with open(join(self.path, 'batch.pkl'), 'wb') as file:
            pickle.dump(self, file, protocol=-1)

        # build parameter file for each batch
        self.build_batches(batch_size=batch_size)

        # build job run script
        self.build_run_script(self.path,
                              self.script_name,
                             num_trajectories,
                             saveall,
                             deviations,
                             comparison)

        # build job submission script
        self.build_submission_script(self.path,
                                     self.script_name,
                                     num_trajectories,
                                     saveall,
                                     deviations,
                                     comparison,
                                     walltime=walltime,
                                     allocation=allocation)

    @classmethod
    def build_simulation(cls, parameters, simulation_type, simulation_path, **kwargs):
        """
        Builds and saves a simulation instance for a set of parameters.

        Args:

            parameters (iterable) - parameter sets

            simulation_type (str) - either "repressors" or "promoters"

            simulation_path (str) - simulation path

            kwargs: keyword arguments for ConditionSimulation

        """

        # build model
        model = cls.build_model(parameters, simulation_type=simulation_type)

        # instantiate simulation
        simulation = ConditionSimulation(model, **kwargs)

        # create simulation directory
        if not isdir(simulation_path):
            mkdir(simulation_path)

        # save simulation
        simulation.save(simulation_path)

    def load_simulation(self, index):
        """
        Load simulation instance from file.

        Args:

            index (int) - simulation index

        Returns:

            simulation (ConditionSimulation)

        """
        simulation_path = join(self.path, self.simulation_paths[index])
        return ConditionSimulation.load(simulation_path)

    def apply(self, func):
        """
        Applies function to entire batch of simulations.

        Args:

            func (function) - function operating on a simulation instance

        Returns:

            output (dict) - {simulation_id: function output} pairs

        """
        f = lambda path: func(ConditionSimulation.load(path))
        return {i: f(p) for i, p in self.simulation_paths.items()}
