from genessa.models.twostate import TwoStateCell
from .mutation import Mutation


class TwoStateModel(TwoStateCell, Mutation):
    """
    Class defines a cell with a single protein coding gene. Transcription is based on a two-state model.

    Attributes:

        name (str) - name of controlled gene

    Inherited Attributes:

        off_states (dict) - {name: node_id} pairs

        on_states (dict) - {name: node_id} pairs

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

    def __init__(self, name='X', k0=1, k1=1, k2=1, g0=1, g1=1, g2=1, include_activation=False):
        """
        Instantiate the two state model.

        Args:

            name (str) - name of controlled gene

            k0 (float) - gene activation rate constant (by input)

            k1 (float) - transcription rate constant

            k2 (float) - translation rate constant

            g0 (float) - gene deactivation rate constant

            g1 (float) - transcript decay rate constant

            g2 (float) - protein decay rate constant

            include_activation (bool) - indicates whether or not to include activation

        """

        self.name = name

        # instantiate twostate cell with a single gene
        gene_kw = dict(k0=0, k1=k1, k2=k2, g0=g0, g1=g1, g2=g2)
        super().__init__(genes=(self.name,), I=1, **gene_kw)

        # add transcriptional activation by input
        if include_activation:
            self.add_activation(activator='IN', k=k0)

    def add_transcriptional_feedback(self,
            k=None,
            atp_sensitive=2,
            carbon_sensitive=2,
            ribosome_sensitive=1,
            **kwargs
        ):
        """
        Adds transcriptional auto-repression.

        Args:

            k (float) - rate parameter (feedback strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """

        self.add_transcriptional_repressor(
             actuator=self.name,
             target=self.name,
             k=k,
             atp_sensitive=atp_sensitive,
             carbon_sensitive=carbon_sensitive,
             ribosome_sensitive=ribosome_sensitive,
             **kwargs)

    def add_post_transcriptional_feedback(self,
            k=None,
            atp_sensitive=2,
            carbon_sensitive=2,
            ribosome_sensitive=1,
            **kwargs
        ):
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
            **kwargs
        ):
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

    def add_activation(self,
            activator='IN',
            k=None,
            atp_sensitive=False,
            carbon_sensitive=False,
            ribosome_sensitive=False,
            **kwargs
        ):
        """
        Adds linear promoter applied to activated-DNA level.

        Args:

            k (float) - rate parameter (promoter strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """
        super().add_activation(
            gene=self.name, 
            activator=activator, 
            k=k,
            atp_sensitive=atp_sensitive,
            carbon_sensitive=carbon_sensitive,
            ribosome_sensitive=ribosome_sensitive,
            **kwargs,
        )

    def add_transcriptional_promoter(self,
            k=None,
            atp_sensitive=True,
            carbon_sensitive=True,
            ribosome_sensitive=False,
            **kwargs
        ):
        """
        Adds linear transcriptional promoter.

        Args:

            k (float) - rate parameter (promoter strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            kwargs: keyword arguments for reaction

        """
        super().add_transcriptional_promoter(
            gene=self.name, 
            k=k,
            atp_sensitive=atp_sensitive,
            carbon_sensitive=carbon_sensitive,
            ribosome_sensitive=ribosome_sensitive,
            **kwargs,
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
        Adds linear translational promoter.

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

    def add_feedback(self, eta0, eta1, eta2, perturbed=False):
        """
        Add feedback at the gene, transcript, and protein levels.

        Args:

            eta0 (float) - transcriptional feedback strength

            eta1 (float) - post-transcriptional feedback strength

            eta2 (float) - post-translational feedback strength

            perturbed (bool) - if True, feedback is sensitive to perturbation

        """
        self.add_transcriptional_feedback(k=eta0, perturbed=perturbed)
        self.add_post_transcriptional_feedback(k=eta1, perturbed=perturbed)
        self.add_post_translational_feedback(k=eta2, perturbed=perturbed)

    def add_promoters(self, eta0, eta1, eta2, perturbed=False):
        """
        Add a promoter at the gene, transcript, and protein levels.

        Args:

            eta0 (float) - activation strength

            eta1 (float) - transcriptional promoter strength

            eta2 (float) - translational promoter strength

            perturbed (bool) - if True, promoters are sensitive to perturbation

        """
        self.add_activation(k=eta0, perturbed=perturbed)
        self.add_transcriptional_promoter(k=eta1, perturbed=perturbed)
        self.add_translational_promoter(k=eta2, perturbed=perturbed)
