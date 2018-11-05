import numpy as np
import pandas as pd

from ..execution.batch import Batch
from .sampling import LogSampler
from ..models.linear import LinearModel
from ..models.hill import HillModel
from ..models.twostate import TwoStateModel
from .figure import SweepFigure


class Sweep(Batch):
    """
    Class defines a parameter sweep of a given model.

    Attributes:

        base (np.ndarray[float]) - base parameter values

        delta (float or np.ndarray[float]) - log-deviations about base

        labels (list of str) - labels for each parameter

    Inherited attributes:

        path (str) - path to batch directory

        parameters (np.ndarray[float]) - sampled parameter values

        simulation_paths (dict) - paths to simulation directories

        sim_kw (dict) - keyword arguments for simulation

        results (dict) - {simulation_id: results_dict} pairs

    Properties:

        N (int) - number of samples in parameter space

    """

    def __init__(self, base, delta=0.5, num_samples=1000, labels=None, pad=.1):
        """
        Instantiate parameter sweep.

        Args:

            base (np.ndarray[float]) - base parameter values

            delta (float or np.ndarray[float]) - log-deviations about base

            num_samples (int) - number of samples in parameter space

            labels (list of str) - labels for each parameter

            pad (float) - extra padding added to delta

        """

        self.base = base
        self.delta = delta
        self.pad = pad
        self.labels = labels

        # sample parameter space
        sampler = LogSampler(base-delta-pad, base+delta+pad)
        parameters = sampler.sample(num_samples)

        # instantiate batch job
        super().__init__(parameters=parameters)

    @staticmethod
    def parse_simulation(simulation):
        """ Returns over, under, and total error from <simulation>. """
        if simulation.comparisons is None:
            return False
        else:
            errors = {}
            for condition, comparison in simulation.comparisons.items():
                errors[(condition, 'above')] = comparison.above
                errors[(condition, 'below')] = comparison.below
                errors[(condition, 'error')] = comparison.error
            return errors

    def aggregate(self):
        """ Aggregate results from each completed simulation. """

        # parse simulation results
        results = [self.parse_simulation(sim) for sim in self]
        error_dicts = [x for x in results if x != False]
        self.incomplete = sum([x for x in results if x == False])/self.N

        # compile results dataframe
        self.df = pd.DataFrame.from_dict(error_dicts, orient='columns')
        self.df.columns = pd.MultiIndex.from_tuples(self.df.columns)

    def slice_by_mode(self, mode='error'):
        """ Returns all results for a specified <mode>. """
        return self.df.swaplevel(axis=1)[mode]

    def slice_by_condition(self, condition='normal'):
        """ Returns all results for a specified <condition>. """
        return self.df[condition]

    def build_figure(self, condition='normal', mode='error', **kwargs):
        """
        Returns parameter sweep visualization.

        Args:

            condition (str) - environmental condition

            mode (str) - comparison metric

            kwargs: keyword arguments for SweepFigure

        """
        return SweepFigure(self.parameters,
                           self.df.loc[:, (condition, mode)],
                           labels=self.labels,
                           **kwargs)


class LinearSweep(Sweep):

    """
    Parameter sweep for linear model. Parameters are:

        0: activation rate constant
        1: transcription rate constant
        2: translation rate constant
        3: deactivation rate constant
        4: mrna degradation rate constant
        5: protein degradation rate constant
        6: transcriptional feedback strength
        7: post-transcriptional feedback strength
        8: post-translational feedback strength

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
            base = np.array([0, 0, 0, 0, -2, -3, -4, -4, -4])

        # define parameter labels
        labels = ('k_0', 'k_1', 'k_2',
                  '\gamma_0', '\gamma_1', '\gamma_2',
                  '\eta_0', '\eta_1', '\eta_2')

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

        Returns:

            model (LinearModel)

        """

        # extract parameters
        k0, k1, k2, g0, g1, g2, eta0, eta1, eta2 = parameters

        # instantiate base model
        model = LinearModel(k0=k0, k1=k1, k2=k2, g0=g0, g1=g1, g2=g2)

        # add feedback (two equivalent sets)
        model.add_feedback(eta0, eta1, eta2, perturbed=False)
        model.add_feedback(eta0, eta1, eta2, perturbed=True)

        return model


class HillSweep(Sweep):

    """
    Parameter sweep of a hill model. Parameters are:

        0: transcription hill coefficient
        1: transcription rate constant
        2: translation rate constant
        3: mrna degradation rate constant
        4: protein degradation rate constant
        5: repressor michaelis constant
        6: repressor hill coefficient
        7: post-transcriptional feedback strength
        8: post-translational feedback strength

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
            base = np.array([0, 0, 0, -2, -3, 4, 0, -5, -4])

        # define parameter labels
        labels = ('H', 'k_R', 'k_P',
                  '\gamma_R', '\gamma_P',
                  'K_r', 'H_r',
                  '\eta_R', '\eta_P')

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

        Returns:

            model (HillModel)

        """

        # extract parameters
        n, k1, k2, g1, g2, k_m, r_n, eta1, eta2 = parameters

        # instantiate base model
        model = HillModel(k1=k1, k_m=.5, n=n, k2=k2, g1=g1, g2=g2)

        # add feedback (two equivalent sets)
        model.add_feedback(k_m, r_n, eta1, eta2, perturbed=False)
        model.add_feedback(k_m, r_n, eta1, eta2, perturbed=True)

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
        6: transcriptional feedback strength
        7: post-transcriptional feedback strength
        8: post-translational feedback strength

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
            base = np.array([0, 0, 0, -1, -2, -3, -4, -4.5, -4])

        # define parameter labels
        labels = ('k_G', 'k_R', 'k_P',
                  '\gamma_G', '\gamma_R', '\gamma_P',
                  '\eta_G', '\eta_R', '\eta_P')

        # call parent instantiation
        super().__init__(base, delta, num_samples, labels=labels)

    @staticmethod
    def build_model(parameters):
        """
        Returns a model instance defined by the provided parameters.

        Args:

            parameters (np.ndarray[float]) - model parameters

        Returns:

            model (LinearModel)

        """

        # extract parameters
        k0, k1, k2, g0, g1, g2, eta0, eta1, eta2 = parameters

        # instantiate base model
        model = TwoStateModel(k0=k0, k1=k1, k2=k2, g0=g0, g1=g1, g2=g2)

        # add feedback (two equivalent sets)
        model.add_feedback(eta0, eta1, eta2, perturbed=False)
        model.add_feedback(eta0, eta1, eta2, perturbed=True)

        return model
