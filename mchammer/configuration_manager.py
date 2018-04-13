class ConfigurationManager:
    """
    ConfigurationManager store its own state of the 
    configuration that are being sampled in the ensemble.

    ConfigurationManager is responsible for all information and
    handling of the configuration.

    Parameters
    ----------
    occupations : list of int
        The occupation of the configurations. The integers
        refer to elements via their number in the periodic table
    strict_occupations : list of list of int
        the strictest form of the allowed occupations.
    occupation_constraints : list of list of int
        optional occupation constraint to enfore a more stricter species
        occupation than what is allowed from the Calculator.
    sublattices : list of list of int
        the integers refer to lattice site indices in the configuration.
    """

        def __init__(self, occupations, strict_occupations, sublattices,  occupation_constraints=None):
            
            self._occupations = occupations            
            self._sublattices = sublattices

            self.
            # internal representation
            # (strings are probably not preferable in actual implementation;
            # they are used here for clarity)
            # ConfigurationManager._sublattices = [['Pd', 'Au'],
            #                                      ['H',  'X']]
            # ConfigurationManager._sites = {'Pd': [0, 1, 2, 3],
            #                                'Au': [4, 5, 6, 7],
            #                                'H':  [8, 9, 10],
            #                                'X':  [11, 12, 13, 14]}


        @property
        def occupations(self):
            """ The occupations of the configuration."""
            return self._occupations.copy()

        @property
        def sublattices(self):
            """ The sublattices of the configuration. """
            return self._sublattices.copy()

        def get_swap_indices(self, sublattice: int):
            """
            Returns a tuple of two random indices in a specific sublattice.

            Parameters
            ----------
            sublattice : int
                The sublattice to pick indices from.            
            """
            raise NotImplementedError

        def get_flip_index(self, sublattice: int):
            """
            Returns a tuple of a random index in a specific sublattice
            and the element to flip it to.

            Parameters
            ----------
            sublattice : int
                The sublattice to pick indices from.
            """
            raise NotImplementedError
