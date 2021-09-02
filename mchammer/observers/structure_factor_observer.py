from collections.abc import Sequence
from typing import Dict, List, Tuple
import numpy as np
from ase import Atoms
from ase.data import chemical_symbols
from mchammer.observers.base_observer import BaseObserver
from icet.input_output.logging_tools import logger
logger = logger.getChild('structure_factor_observer')


class StructureFactorObserver(BaseObserver):
    r"""This class represents a structure factor observer.

    This observer allows one to compute structure factors along the
    trajectory sampled by a Monte Carlo (MC) simulation. Structure
    factors are convenient for monitoring long-range order. The
    `structure factor <https://en.wikipedia.org/wiki/Structure_factor>`_
    is defined as:

    .. math::

        S(\vec{q}) = \frac{1}{\sum_{j=1}^N f_j^2}
                 \sum_{j,k}^N e^{-i \vec{q} \cdot (\vec{R}_k - \vec{R}_j)}

    In addition to this "total" structure factor, this observer
    calculates pair-specific structure factors, which correspond to parts
    of the summation defined above, with the summation restricted to pairs
    of specific types, e.g., Au-Au, Au-Cu and Cu-Cu in the example below.

    Parameters
    ----------
    structure : Atoms
        prototype for the structures for which the structure factor will
        be computed later; the supercell size (but not its decoration)
        must be identical; this structure is also used to determine the
        the possible pairs if `symbol_pairs=None`

    q_points : List[np.ndarray]
        array of q-points at which to evaluate the structure factor;
        the q-points need to be compatible with the supercell in order
        for the structure factor to be real

    symbol_pairs : Optional[List[Tuple[str, str]]]
        list of symbol pairs for which structure factors will be computed,
        e.g. [('Al', 'Cu'), ('Al', 'Al')]; if `None` (default) use all
        pairs possible based on the input structure

    form_factors : Optional[Dict[str, float]]
        form factors for each atom type; this can be used to (very
        coarsely) simulate X-ray or neutron spectra; note that in general
        the form factors are q-dependent, see, e.g.,
        `here <https://wiki.fysik.dtu.dk/ase/ase/xrdebye.html>`__ and
        `here <https://wiki.fysik.dtu.dk/ase/ase/
        xrdebye.html#ase.utils.xrdebye.XrDebye.get_waasmaier>`__;
        by default (`None`) all form factors are set to 1

    interval : int
        the observation interval, defaults to None meaning that if the
        observer is used in a Monte Carlo simulation, the `Ensemble` object
        will set the interval.

    Raises
    ------
    ValueError
        if q-point is not consistent with metric of input structure

    Attributes
    ----------
    tag : str
        name of observer

    interval : int
        observation interval

    Example
    -------
    The following snippet illustrates how to use the structure factor
    observer in a simulated annealing run of dummy Cu-Au model to
    observe the emergence of a long-range ordered L1_2 structure::

        >>> import numpy as np
        >>> from ase.build import bulk
        >>> from icet import ClusterSpace, ClusterExpansion
        >>> from mchammer.calculators import ClusterExpansionCalculator
        >>> from mchammer.ensembles import CanonicalAnnealing
        >>> from mchammer.observers import StructureFactorObserver

        >>> # parameters
        >>> size = 3
        >>> alat = 4.0
        >>> symbols = ['Cu', 'Au']

        >>> # setup
        >>> prim = bulk('Cu', a=alat, cubic=True)
        >>> cs = ClusterSpace(prim, [0.9*alat], symbols)
        >>> ce = ClusterExpansion(cs, [0, 0, 0.2])

        >>> # make supercell
        >>> supercell = prim.repeat(size)
        >>> ns = int(0.25 * len(supercell))
        >>> supercell.symbols[0:ns] = 'Au'
        >>> np.random.shuffle(supercell.symbols)

        >>> # define q-points to sample
        >>> q_points = []
        >>> q_points.append(2 * np.pi / alat * np.array([1, 0, 0]))
        >>> q_points.append(2 * np.pi / alat * np.array([0, 1, 0]))
        >>> q_points.append(2 * np.pi / alat * np.array([0, 0, 1]))

        >>> # set up structure factor observer
        >>> sfo = StructureFactorObserver(supercell, q_points)

        >>> # run simulation
        >>> calc = ClusterExpansionCalculator(supercell, ce)
        >>> mc = CanonicalAnnealing(supercell, calc,
        ...                         T_start=900, T_stop=500, cooling_function='linear',
        ...                         n_steps=400*len(supercell))
        >>> mc.attach_observer(sfo)
        >>> mc.run()

    After having run this snippet, one can access the structure factors via the data
    container::

        >>> dc = mc.data_container
        >>> print(dc.data)

    The emergence of the ordered low-temperature structure can be monitored by
    following the temperature dependence of any of the pair-specific structure
    factors.
    """

    def __init__(self,
                 structure: Atoms,
                 q_points: List[Sequence],
                 symbol_pairs: List[Tuple[str, str]] = None,
                 form_factors: Dict[str, float] = None,
                 interval: int = None) -> None:
        super().__init__(interval=interval, return_type=dict, tag='StructureFactorObserver')

        self._original_structure = structure.copy()
        self._q_points = q_points
        self._Sq_lookup = self._get_Sq_lookup(structure, q_points)

        if symbol_pairs is not None:
            self._pairs = symbol_pairs
        else:
            self._pairs = [(e1, e2)
                           for e1 in set(structure.symbols)
                           for e2 in set(structure.symbols)]
        self._pairs = sorted(set([tuple(sorted(p)) for p in self._pairs]))
        if len(self._pairs) == 0:
            raise ValueError('The StructureFactorObserver requires '
                             'at least one pair to be defined')
        if len(self._pairs) == 1 and self._pairs[0][0] == self._pairs[0][1]:
            logger.warning(f'Only one pair requested {self._pairs[0]}; '
                           'use either a structure occupied by more than '
                           'one species or use the symbol_pairs argument '
                           'to specify more pairs.')

        self._unique_symbols = set(sorted([s for p in self._pairs for s in p]))
        if form_factors is not None:
            self._form_factors = form_factors
        else:
            self._form_factors = {s: 1 for s in self._unique_symbols}
        self._form_factors = dict(sorted(self._form_factors.items()))

        for symbol in self._unique_symbols:
            if symbol not in chemical_symbols:
                raise ValueError(f'Unknown species {symbol}')
            if symbol not in self._form_factors:
                raise ValueError(f'Form factor missing for {symbol}')

        for symbol, ff in self._form_factors.items():
            if ff == 0:
                raise ValueError(f'Form factor for {symbol} is zero')

    def _get_Sq_lookup(self,
                       structure: Atoms,
                       q_points: List[Sequence]) -> np.ndarray:
        """ Get S(q) lookup data for a given supercell and q-points. """
        n_atoms = len(structure)
        dist_vectors = structure.get_all_distances(mic=True, vector=True)
        Sq_lookup = np.zeros((n_atoms, n_atoms, len(self._q_points)), dtype=np.complex128)
        for i in range(n_atoms):
            for j in range(n_atoms):
                Rij = dist_vectors[i][j]
                Sq_lookup[i, j, :] = np.exp(-1j * np.dot(q_points, Rij))
        return Sq_lookup

    def _get_indices(self,
                     structure: Atoms) -> Dict[str, List[int]]:
        """Returns the indices of all atoms by type considering all unique
        atom types in the structure used to instantiate the observer object.
        """
        indices = dict()
        symbols = structure.get_chemical_symbols()
        for symbol in self._unique_symbols:
            indices[symbol] = np.where(np.array(symbols) == symbol)[0]
        return indices

    def _compute_structure_factor(self,
                                  structure: Atoms) -> Dict[Tuple[str, str], float]:
        indices = self._get_indices(structure)
        norm = np.sum([self._form_factors[sym] ** 2 * len(lst) for sym, lst in indices.items()])
        if norm <= 0:
            prefactor = 0  # occurs if none of the specified atom types are present in the structure
        else:
            prefactor = 1 / norm

        Sq_dict = dict()
        for sym1, sym2 in self._pairs:
            inds1 = indices[sym1]
            inds2 = indices[sym2]
            norm = prefactor * self._form_factors[sym1] * self._form_factors[sym2]
            Sq = np.zeros(len(self._q_points), dtype=np.complex128)
            if len(inds1) != 0 and len(inds2) != 0:
                for i in inds1:
                    Sq += self._Sq_lookup[i, inds2].sum(axis=0) * norm
            if sym1 == sym2:
                # the imaginary part must be zero because both ij and ji appear in the sum
                Sq_dict[(sym1, sym2)] = Sq.real
            else:
                # since we only consider A-B (not B-A; see __init__) we explicitly add ij and ji,
                # which comes down to adding the complex conjugate
                Sq_dict[(sym1, sym2)] = np.real(Sq + Sq.conj())
        return Sq_dict

    def get_observable(self,
                       structure: Atoms) -> Dict[str, float]:
        """
        Returns the structure factors for a given atomic configuration.

        Parameters
        ----------
        structure
            input atomic structure

        Raises
        ------
        ValueError
            if input structure is incompatible with structure used for initialization
        """
        if len(structure) != len(self._original_structure):
            raise ValueError('Input structure incompatible with structure used for initialization\n'
                             f'  n_input: {len(structure)}\n'
                             f'  n_original: {len(self._original_structure)}')
        Sq_dict = self._compute_structure_factor(structure)
        return_dict = dict()
        for pair, Sq in Sq_dict.items():
            for i in range(len(Sq)):
                tag = 'sfo_{}_{}_q{}'.format(*pair, i)
                return_dict[tag] = Sq[i]
                return_dict[f'total_q{i}'] = return_dict.get(f'total_q{i}', 0) + Sq[i]
        return_dict = dict(sorted(return_dict.items()))
        return return_dict

    @property
    def form_factors(self) -> Dict[str, float]:
        """ Form factors used in structure factor calculation """
        return self._form_factors.copy()

    @property
    def q_points(self) -> List[np.ndarray]:
        """ q-points for which structure factor is calculated """
        return self._q_points.copy()

    def __str__(self) -> str:
        """ string representation of observer object. """
        width = 60
        name = self.__class__.__name__
        s = [' {} '.format(name).center(width, '=')]

        fmt = '{:15} : {}'
        s += [fmt.format('tag', self.tag)]
        s += [fmt.format('return_type', self.return_type)]
        s += [fmt.format('interval', self.interval)]
        s += [fmt.format('pairs', self._pairs)]
        s += [fmt.format('form_factors', self.form_factors)]
        for k, qpt in enumerate(self.q_points):
            s += [fmt.format(f'q-point{k}', qpt)]
        return '\n'.join(s)
