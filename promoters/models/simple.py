from genessa.models.simple import SimpleCell
from .mutation import Mutation


class SimpleModel(SimpleCell, Mutation):
    """
    Class defines a cell with a single protein state subject to negative feedback. All reaction rates are based on linear propensity functions.

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

    def __init__(self,
                 name='X',
                 k=1,
                 g=1,
                 lambda_g=0,
                 lambda_k=0,
                 include_activation=True):
        """
        Instantiate a simple model of a single protein.

        Args:

            name (str) - name of controlled protein

            k (float) - protein synthesis rate constant

            g (float) - protein decay rate constant

            lambda_g (float) - degradation growth rate dependence

            lambda_k (float) - synthesis growth rate dependence

            include_activation (bool) - indicates whether or not to include activation

        """

        self.name = name

        # instantiate linear cell with a single gene activated by the input
        gene_kw = dict(g=g, growth_dependence=lambda_g)
        super().__init__(genes=(self.name,), I=1, **gene_kw)

        # add synthesis driven by input
        if include_activation:
            self.add_activation(protein=self.name,
                                activator='IN',
                                k=k,
                                growth_dependence=lambda_k)

    def add_post_translational_feedback(self,
            k=None,
            atp_sensitive=2,
            carbon_sensitive=2,
            ribosome_sensitive=1,
            growth_dependence=0,
            **kwargs):
        """
        Adds linear negative feedback applied to protein level.

        Args:

            k (float) - rate parameter (feedback strength)

            atp_sensitive (int) - order of metabolism dependence

            carbon_sensitive (int) - order of carbon availability dependence

            ribosome_sensitive (int) - order of ribosome dependence

            growth_dependence (float) - log k / log growth

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
             growth_dependence=growth_dependence,
             **kwargs)

    def add_feedback(self, eta, perturbed=False, lambda_eta=0):
        """
        Add feedback.

        Args:

            eta (float) - feedback strength

            perturbed (bool) - if True, feedback is sensitive to perturbation

            lambda_eta (float) - feedback growth rate dependence

        """
        self.add_post_translational_feedback(
            k=eta,
            perturbed=perturbed,
            growth_dependence=lambda_eta,
        )

    def add_promoters(self, eta, perturbed=False, lambda_eta=0):
        """
        Add promoter.

        Args:

            eta (float) - promoter strength

            perturbed (bool) - if True, promoter is sensitive to perturbation

            lambda_eta (float) - promoter growth rate dependence

        """
        self.add_activation(
            protein=self.name,
            activator='IN',
            k=eta,
            perturbed=perturbed,
            growth_dependence=lambda_eta
        )
