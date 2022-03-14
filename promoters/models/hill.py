from genessa.models.hill import HillCell
from .mutation import Mutation


class HillModel(HillCell, Mutation):
    """
    Class defines a cell with a single protein coding gene subject to negative feedback. Transcription and feedback terms are based on Hill kinetics.

    Attributes:

        name (str) - name of controlled gene

    Inherited Attributes:

        transcripts (dict) - {name: node_id} pairs

        proteins (dict) - {name: node_id} pairs

        phosphorylated (dict) - {name: node_id} pairs

        nodes (np.ndarray) - vector of node indices

        node_key (dict) - {state dimension: node id} pairs

        reactions (list) - list of reaction objects

        stoichiometry (np.ndarray) - stoichiometric coefficients, (N,M)

        N (int) - number of nodes

        M (int) - number of reactions

        I (int) - number of inputs

    """

    def __init__(self, name='X', k1=1, k_m=1, n=1, k2=1, g1=1, g2=1, include_activation=False):
        """
        Instantiate the Hill model.

        Args:

            name (str) - name of controlled gene

            k1 (float) - maximal transcription rate

            k_m (float) - michaelis constant

            n (float) - hill coefficient

            k2 (float) - translation rate constant

            g1 (float) - transcript decay rate constant

            g2 (float) - protein decay rate constant

            include_activation (bool) - indicates whether or not to include activation

        """

        self.name = name

        # instantiate Hill cell with a single gene
        gene_kw = dict(k=k2, g1=g1, g2=g2)
        super().__init__(genes=(self.name,), I=1, **gene_kw)

        # add transcriptional activation by input
        if include_activation:
            self.add_transcription(self.name, ('IN',), k=k1, k_m=k_m, n=n)

    def add_transcriptional_feedback(self, k_m=1, n=1, **kwargs):
        """
        Adds transcriptional auto-repression.

        Args:

            k_m (float) - michaelis constant

            n (float) - hill coefficient

            kwargs: keyword arguments for repressor

        """
        self.add_transcriptional_repressor(
             actuators=(self.name,),
             target=self.name,
             k_m=k_m,
             n=n,
             **kwargs)

    def add_post_transcriptional_feedback(self,
                                     k=None,
                                     atp_sensitive=2,
                                     carbon_sensitive=2,
                                     ribosome_sensitive=1,
                                     **kwargs):
        """
        Adds linear negative feedback applied to transcript level.

        Args:

            k (float) - rate parameter (feedback strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """

        self.add_linear_feedback(
             sensor=self.name,
             target=self.name,
             mode='transcript',
             k=k,
             atp_sensitive=atp_sensitive,
             carbon_sensitive=carbon_sensitive,
             ribosome_sensitive=ribosome_sensitive,
             **kwargs)

    def add_post_translational_feedback(self,
                                     k=None,
                                     atp_sensitive=2,
                                     carbon_sensitive=2,
                                     ribosome_sensitive=1,
                                     **kwargs):
        """
        Adds linear negative feedback applied to protein level.

        Args:

            k (float) - rate parameter (feedback strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """

        self.add_linear_feedback(
             sensor=self.name,
             target=self.name,
             mode='protein',
             k=k,
             atp_sensitive=atp_sensitive,
             carbon_sensitive=carbon_sensitive,
             ribosome_sensitive=ribosome_sensitive,
             **kwargs)

    def add_feedback(self, k_m, n, eta1, eta2, perturbed=False):
        """
        Add feedback at the gene, transcript, and protein levels.

        Args:

            k_m (float) - repressor michaelis constant

            n (float) - repressor hill coefficient

            eta1 (float) - post-transcriptional feedback strength

            eta2 (float) - post-translational feedback strength

            perturbed (bool) - if True, feedback is sensitive to perturbation

        """
        self.add_transcriptional_feedback(k_m=k_m, n=n, perturbed=perturbed)
        self.add_post_transcriptional_feedback(k=eta1, perturbed=perturbed)
        self.add_post_translational_feedback(k=eta2, perturbed=perturbed)

    def add_transcriptional_promoter(self,
            k,
            k_m,
            n,
            baseline=0.,
            atp_sensitive=True,
            carbon_sensitive=True,
            ribosome_sensitive=False,
            **kwargs
        ):
        """
        Adds a transcriptional promoter.

        Args:

            k (float) - maximum transcription rate

            k_m (float) - michaelis menten constant

            n (float) - hill coefficients

            baseline (float) - baseline transcription rate

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """
        self.add_transcription(
            self.name, 
            ('IN',), 
            k=k, 
            k_m=k_m, 
            n=n,
            baseline=baseline,
            atp_sensitive=atp_sensitive,
            carbon_sensitive=carbon_sensitive,
            ribosome_sensitive=ribosome_sensitive,
            **kwargs
        )

    def add_translational_promoter(
            self,
            k=None,
            atp_sensitive=True,
            carbon_sensitive=True,
            ribosome_sensitive=True,
            **kwargs
        ):
        """
        Adds linear promoter applied to transcript level.

        Args:

            k (float) - rate parameter (promoter strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """
        super().add_translational_promoter(
            gene=self.name,
            k=k,
            atp_sensitive=atp_sensitive,
            carbon_sensitive=carbon_sensitive,
            ribosome_sensitive=ribosome_sensitive,
            **kwargs
        )

    def add_promoters(self, k_1, k_m, n, k_2, baseline=0.0, perturbed=False):
        """
        Add transcriptional and translational promoters

        Args:

            k_1 (float) - transcriptional promoter strength

            k_m (float) - promoter michaelis constant

            n (float) - promoter hill coefficient

            k_2 (float) - translational promoter strength

            baseline (float) - baseline transcription rate

            perturbed (bool) - if True, feedback is sensitive to perturbation

        """
        self.add_transcriptional_promoter(k=k_1, k_m=k_m, n=n, baseline=baseline, perturbed=perturbed)
        self.add_translational_promoter(k=k_2, perturbed=perturbed)
