import numpy as np
import pandas as pd
from os.path import join, exists

from ..execution.batch import Batch
from ..models.linear import LinearModel
from ..models.hill import HillModel
from ..models.twostate import TwoStateModel
from ..models.simple import SimpleModel
from .sampling import SobolLogSampler, SobolLinearSampler, DenseLinearSampler
from .figure import SweepFigure, SweepHistogram, SweepLines


class Sweep(Batch):
    """
    Class defines a parameter sweep of a given model.

    Attributes:

        base (np.ndarray[float]) - base parameter values

        delta (float or np.ndarray[float]) - log-deviations about base

        labels (list of str) - labels for each parameter

        results (pd.DataFrame) - results with (condition, mode) multiindex

    Inherited attributes:

        path (str) - path to batch directory

        run_script (str) - path to script for running batch

        parameters (np.ndarray[float]) - sampled parameter values

        simulation_paths (dict) - paths to simulation directories

        sim_kw (dict) - keyword arguments for simulation

    Properties:

        N (int) - number of samples in parameter space

    """

    def __init__(self,
                 base,
                 delta=0.5,
                 num_samples=1000,
                 labels=None,
                 pad=.1,
                 logbase=10):
        """
        Instantiate parameter sweep.

        Args:

            base (np.ndarray[float]) - base parameter values

            delta (float or np.ndarray[float]) - log-deviations about base

            num_samples (int) - number of samples in parameter space

            labels (list of str) - labels for each parameter

            pad (float) - extra padding added to delta

            logbase (float) - logarithmic basis for sampling

        """

        self.base = base
        self.delta = delta
        self.pad = pad
        self.labels = labels
        self.results = None

        # sample parameter space
        sampler = SobolLogSampler(base-delta-pad, base+delta+pad, base=logbase)
        parameters = sampler.sample(num_samples)

        # instantiate batch job
        super().__init__(parameters=parameters)

    @staticmethod
    def load(path):
        """ Load sweep from target <path>. """

        sweep = Batch.load(path)

        # load results
        if exists(join(path, 'data.hdf')):
            sweep.results = pd.read_hdf(join(path, 'data.hdf'), 'results')
            sweep.completed = pd.read_hdf(join(path, 'data.hdf'), 'completed')

        return sweep

    def save(self):
        """ Save sweep data. """
        if self.results is not None:
            p = join(self.path, 'data.hdf')
            self.results.to_hdf(p, 'results', mode='a')
            self.completed.to_hdf(p, 'completed', mode='a')

    @staticmethod
    def parse_simulation(simulation, index=-3):
        """ Returns over, under, and total error from <simulation>. """
        if simulation.comparisons is None:
            return False
        else:
            errors = {}
            for c, comparison in simulation.comparisons.items():
                failed = False

                if comparison.__class__.__name__ == 'MultiComparison':

                    if True in comparison.reached_comparison:
                        errors[(c, 'above')] = comparison.above[index]
                        errors[(c, 'below')] = comparison.below[index]
                        errors[(c, 'error')] = comparison.error[index]
                        errors[(c, 'above_threshold')] = comparison.above_threshold[index]
                        errors[(c, 'below_threshold')] = comparison.below_threshold[index]
                        errors[(c, 'threshold_error')] = comparison.threshold_error[index]
                    else:
                        failed = True

                else:

                    if comparison.reached_comparison:
                        errors[(c, 'above')] = comparison.above
                        errors[(c, 'below')] = comparison.below
                        errors[(c, 'error')] = comparison.error
                        errors[(c, 'above_threshold')] = comparison.above_threshold
                        errors[(c, 'below_threshold')] = comparison.below_threshold
                        errors[(c, 'threshold_error')] = comparison.threshold_error
                    else:
                        failed = True

                if failed:
                    fail = None
                    errors[(c, 'above')] = fail
                    errors[(c, 'below')] = fail
                    errors[(c, 'error')] = fail
                    errors[(c, 'above_threshold')] = fail
                    errors[(c, 'below_threshold')] = fail
                    errors[(c, 'threshold_error')] = fail

            return errors

    @property
    def percent_complete(self):
        """ Fraction of simulations that ran to completion. """
        return self.completed.values.sum()/self.N

    def aggregate(self):
        """
        Aggregate results from each completed simulation and compile them into a dataframe.
        """

        # parse simulation results
        results = [self.parse_simulation(sim) for sim in self]

        # store indices of completed simulations
        index = np.arange(self.N)
        self.completed = pd.DataFrame([x!=False for x in results], index=index)

        # compile results dataframe
        error_dicts = [x for x in results if x != False]
        self.results = pd.DataFrame.from_dict(error_dicts, orient='columns')
        self.results.columns = pd.MultiIndex.from_tuples(self.results.columns)

    def slice_by_mode(self, mode='error'):
        """ Returns all results for a specified <mode>. """
        return self.results.swaplevel(axis=1)[mode]

    def slice_by_condition(self, condition='normal'):
        """ Returns all results for a specified <condition>. """
        return self.results[condition]

    def build_figure(self,
                     condition='normal',
                     mode='error',
                     relative=False,
                     **kwargs):
        """
        Returns parameter sweep visualization.

        Args:

            condition (str or tuple) - environmental condition

            mode (str) - comparison metric

            relative (bool) - if True, computes difference relative to normal

            kwargs: keyword arguments for SweepFigure

        """

        # evaluate results
        results = self.results.loc[:, (condition, mode)]
        if relative:
            results = results - self.results.loc[:, ('normal', mode)]

        return SweepFigure(self.parameters[self.completed.values.ravel()],
                           results,
                           labels=self.labels,
                           base=self.base,
                           delta=self.delta,
                           **kwargs)

    def build_histogram(self,
                     condition='normal',
                     mode='error',
                     relative=False):
        """
        Returns 1D histogram visualization of parameter sweep.

        Args:

            condition (str or tuple) - environmental condition

            mode (str) - comparison metric

            relative (bool) - if True, computes difference relative to normal

        """

        # evaluate results
        results = self.results.loc[:, (condition, mode)]
        if relative:
            results = results - self.results.loc[:, ('normal', mode)]

        return SweepHistogram(results.values)

    def build_lines(self,
                     condition='normal',
                     mode='error',
                     relative=False):
        """
        Returns line projection of parameter sweep as a function of where the success threshold is set.

        Args:

            condition (str or tuple) - environmental condition

            mode (str) - comparison metric

            relative (bool) - if True, computes difference relative to normal

        """

        get_error = lambda x: x.error if mode == 'error' else x.threshold_error
        parse_sim = lambda sim: get_error(sim.comparisons[condition])

        # evaluate results then replace None values with NaN
        values = np.array([parse_sim(sim) for sim in self])
        values[(values==None)] = np.nan
        values = values.astype(float)

        if relative:
            parse_normal = lambda sim: get_error(sim.comparisons['normal'])
            reference_values = np.array([parse_normal(sim) for sim in self])
            reference_values[(reference_values==None)] = np.nan
            reference_values = reference_values.astype(float)
            values = reference_values - values

        # get threshold positions
        fractions_of_max = self[0].comparisons['normal'].fraction_of_max

        return SweepLines(values, fractions_of_max)


class SimpleSweep(Sweep):

    """
    Parameter sweep for simple model. Parameters are:

        0: synthesis rate constant
        1: decay rate constant
        2: auxiliary promoter strength

    """

    def __init__(self, base=None, delta=0.5, num_samples=2500):
        """
        Instantiate parameter sweep of a simple model.

        Args:

            base (np.ndarray[float]) - base parameter values

            delta (float or np.ndarray[float]) - log-deviations about base

            num_samples (int) - number of samples in parameter space

        """

        # define parameter ranges, log10(val)
        if base is None:
            base = np.array([0, -3, -1])

        # define parameter labels
        labels = ('k', '\gamma', 'severity')

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters, simulation_type='promoters'):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

            simulation_type (str) - either "repressors" or "promoters"

        Returns:

            model (LinearModel)

        """

        # extract parameters
        k, g, perturbation_severity = parameters


        if simulation_type == 'promoters':

            # instantiate base model
            model = SimpleModel(k=k, g=g, include_activation=False)

            # add promoters (two equivalent sets)
            model.add_promoters(k, perturbed=False)
            model.add_promoters(k*perturbation_severity, perturbed=True)

        else:

            # instantiate base model
            model = SimpleModel(k=k, g=g, include_activation=True)

            # add feedback (two equivalent sets)
            model.add_feedback(perturbation_severity*g, perturbed=True)

        return model


class LinearSweep(Sweep):

    """
    Parameter sweep for linear model. Parameters are:

        0: activation rate constant
        1: transcription rate constant
        2: translation rate constant
        3: deactivation rate constant
        4: mrna degradation rate constant
        5: protein degradation rate constant
        6: activation promoter strength
        7: transcriptional promoter strength
        8: translational promoter strength

    """

    def __init__(self, base=None, delta=0.5, num_samples=1000):
        """
        Instantiate parameter sweep of a linear model.

        Args:

            base (np.ndarray[float]) - base parameter values

            delta (float or np.ndarray[float]) - log-deviations about base

            num_samples (int) - number of samples in parameter space

        """

        # define parameter ranges, log10(val)
        if base is None:
            base = np.array([0, 0, 0, 0, -2, -3, -1])

        # define parameter labels
        labels = ('k_0', 'k_1', 'k_2',
                  '\gamma_0', '\gamma_1', '\gamma_2',
                  'severity',)

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters, simulation_type='promoters'):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

            simulation_type (str) - either "repressors" or "promoters"

        Returns:

            model (LinearModel)

        """

        # extract parameters
        k0, k1, k2, g0, g1, g2, perturbation_severity = parameters

        if simulation_type == 'promoters':

            # instantiate base model
            model = LinearModel(g0=g0, g1=g1, g2=g2, include_activation=False)

            # add promoters (two equivalent sets)
            model.add_promoters(k0, k1, k2, perturbed=False)
            model.add_promoters(k0*perturbation_severity, k1*perturbation_severity, k2*perturbation_severity, perturbed=True)

        else:

            # instantiate base model
            model = LinearModel(k0=k0, k1=k1, k2=k2, g0=g0, g1=g1, g2=g2, include_activation=True)

            # add feedback (two equivalent sets)
            model.add_feedback(perturbation_severity*g0, perturbation_severity*g1, perturbation_severity*g2, perturbed=True)

        return model


class HillSweep(Sweep):

    """
    Parameter sweep of a hill model. Parameters are:

        0: transcription hill coefficient
        1: transcription rate constant
        2: translation rate constant
        3: mrna degradation rate constant
        4: protein degradation rate constant        
        5: transcriptional promoter strength
        6: promoter michaelis constant
        7: promoter hill coefficient
        8: translational promoter strength

    """

    def __init__(self, base=None, delta=0.5, num_samples=1000):
        """
        Instantiate parameter sweep of a Hill model.

        Args:

            base (np.ndarray[float]) - base parameter values

            delta (float or np.ndarray[float]) - log-deviations about base

            num_samples (int) - number of samples in parameter space

        """

        # define parameter ranges, log10(val)
        if base is None:
            base = np.array([0, 0, 0, -2, -3, -1])

        # define parameter labels
        labels = ('H', 'k_R', 'k_P',
                  '\gamma_R', '\gamma_P',
                  'severity',)

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters, simulation_type='promoters'):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

            simulation_type (str) - either "repressors" or "promoters"

        Returns:

            model (HillModel)

        """

        # extract parameters
        n, k1, k2, g1, g2, perturbation_severity = parameters

        if simulation_type == 'promoters':

            # instantiate base model
            model = HillModel(g1=g1, g2=g2, include_activation=False)

            # add promoters (two equivalent sets)
            model.add_promoters(k1, .5, n, k2, perturbed=False)
            model.add_promoters(k1*perturbation_severity, .5, n, k2*perturbation_severity, perturbed=True)

        else:

            # instantiate base model
            model = HillModel(k1=k1, k_m=.5, n=n, k2=k2, g1=g1, g2=g2, include_activation=True)

            # add feedback (two equivalent sets)
            assert False, 'Not Implemented'
            model.add_feedback(k_m, r_n, g1*perturbation_severity, g2*perturbation_severity, perturbed=True)            

        return model


class TwoStateSweep(Sweep):

    """
    Parameter sweep of a twostate model. Parameters are:

        0: activation rate constant
        1: transcription rate constant
        2: translation rate constant
        3: deactivation rate constant
        4: mrna degradation rate constant
        5: protein degradation rate constant
        6: activation promoter strength
        7: transcriptional promoter strength
        8: translational promoter strength

    """

    def __init__(self, base=None, delta=0.5, num_samples=1000):
        """
        Instantiate parameter sweep of a twostate model.

        Args:

            base (np.ndarray[float]) - base parameter values

            delta (float or np.ndarray[float]) - log-deviations about base

            num_samples (int) - number of samples in parameter space

        """

        # define parameter ranges, log10(val)
        if base is None:
            base = np.array([0, 0, 0, -1, -2, -3, -1])

        # define parameter labels
        labels = ('k_G', 'k_R', 'k_P',
                  '\gamma_G', '\gamma_R', '\gamma_P',
                  'severity',
                  )

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters, simulation_type='promoters'):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

            simulation_type (str) - either "repressors" or "promoters"

        Returns:

            model (LinearModel)

        """

        # extract parameters
        k0, k1, k2, g0, g1, g2, perturbation_severity = parameters

        if simulation_type == 'promoters':

            # instantiate base model
            model = TwoStateModel(g0=g0, g1=g1, g2=g2, include_activation=False)

            # add promoters (two equivalent sets)
            model.add_promoters(k0, k1, k2, perturbed=False)
            model.add_promoters(k0*perturbation_severity, k1*perturbation_severity, k2*perturbation_severity, perturbed=True)

        else:

            # instantiate base model
            model = TwoStateModel(k0=k0, k1=k1, k2=k2, g0=g0, g1=g1, g2=g2, include_activation=True)

            # add feedback (two equivalent sets)
            model.add_feedback(perturbation_severity*g0, perturbation_severity*g1, perturbation_severity*g2, perturbed=True)

        return model
