import itertools
from collections import OrderedDict


class LabelingGenerator():
    """
    Object used to generate all possible permutations of elements on a given
    set of sites. If concentration restrictions are not specified, the
    generation will be a simple itertools.product loop.

    Parameters
    ----------
    iter_elements : list of tuples of ints
        Specifies the allowed elements on each site
    concentrations : dict
        Keys are elements (integers), values specify concentration ranges.
    tol : float
        Tolerance parameter used when comparing concentrations
    """

    def __init__(self, iter_elements, concentrations, tol=1e-5):
        self.concentrations = concentrations
        self.iter_elements = iter_elements
        self.tol = tol

        if self.concentrations:
            self.site_groups = OrderedDict()
            count = 0
            for iter_element_key in iter_elements:
                if iter_element_key in self.site_groups.keys():
                    self.site_groups[iter_element_key].multiplicity += 1
                else:
                    self.site_groups[iter_element_key] = SiteGroup(
                        iter_element_key, count)
                    count += 1

    def yield_labelings(self, ncells):
        """
        Yield labelings that comply with the concentration restrictions and,
        with `ncells` primitive cells.

        Parameters
        ----------
        ncells : int
            Size of supercell

        Yields
        ------
        tuple of ints
            Labeling
        """

        if self.concentrations:
            for site_group in self.site_groups.values():
                site_group.compute_all_combinations(ncells)
            for product in self.yield_products(ncells):
                for labeling in self.yield_permutations(product, 0):
                    yield self.sort_labeling(labeling, ncells)
        else:
            for labeling in itertools.product(*self.iter_elements * ncells):
                yield labeling

    def yield_products(self, ncells):
        """
        Yield combinations (or rather products in the itertools terminology)
        of decorated site group combinations that comply with the concentration
        restrictions.

        Returns
        -------
        tuple of tuples of ints
            All unique combinations (products) of unordered iter element
            combinations,
        """
        natoms = len(self.iter_elements) * ncells
        combinations = [sg.combinations for sg in self.site_groups.values()]
        for product in itertools.product(*combinations):
            counts = {element: 0 for element in self.concentrations.keys()}
            for element_group in product:
                for element in self.concentrations.keys():
                    counts[element] += element_group.count(element)
            for element, conc_range in self.concentrations.items():
                if counts[element] / natoms > conc_range[0] - self.tol and \
                   counts[element] / natoms < conc_range[1] + self.tol:
                    yield product

    def yield_unique_permutations(self, unique_elements, permutation, position):
        """
        Recursively generate all _unique_ permutations of a set of elements with
        given multiplicites. The algorithm is inspired by the one given at
        https://stackoverflow.com/a/6285203/6729551

        Parameters
        ----------
        unique_elements : dict
            The keys are elements, the values the multiplicity for that element.
        permutation : list of ints
            Permutation in process, should have length equal to the sum of all
            multiplicites.
        position : int
            Position currently processed. Should be the index of the last element
            of `permutation` upon initialization.

        Yields
        ------
        tuple with ints
            Permutation where each element occur according to the multiplicity
            specified in `unique_elements`
        """
        if position < 0:
            # Finish recursion
            yield tuple(permutation)
        else:
            for element, occurrences in unique_elements.items():
                # Add if the multiplicity allows
                if occurrences > 0:
                    permutation[position] = element
                    unique_elements[element] -= 1
                    for perm in self.yield_unique_permutations(unique_elements,
                                                               permutation,
                                                               position - 1):
                        yield perm
                    unique_elements[element] += 1

    def yield_permutations(self, product, position):
        """
        Recursively generate all combinations of unique permutations of the
        tuples in `product`.

        Parameters
        ----------
        product : tuple of tuples of ints
            Elements allowed for each site group
        position : int
            Keeps track of the position where recursion occurs. Set to 0 upon
            initilization.

        Yields
        ------
        list of tuples of ints
            Unique combinations of unqiue permutations, ordered by site group
        """
        unique_elements = {element: product[position].count(element)
                           for element in set(product[position])}
        natoms = len(product[position])
        for permutation in self.yield_unique_permutations(unique_elements,
                                                          [0] * natoms,
                                                          natoms - 1):
            if position == len(product) - 1:
                yield [permutation]
            else:
                for permutation_rest in self.yield_permutations(product,
                                                                position + 1):
                    yield [permutation] + permutation_rest

    def sort_labeling(self, labeling, ncells):
        """
        The elements in labeling are now given in site groups. We want just
        one labeling, ordered by site group, but one primitive cell at a time.
        So if iter_elements given upon init was  `[(0, 1), (2), (0, 1)]` and
        ncells = 2, the current labeling has two tuples, the first
        corrresponding to the `(0, 1)` site group, with 8 elements in total,
        and the second corresponding to the `(2)` site group and with 2
        elements in total. Now, we reorder it so that the elements are ordered
        according to `[(0, 1), (2), (0, 1), (0, 1), (2), (0, 1)]`

        Parameters
        ----------
        labeling : list of tuple of ints
            As given by yield_permutations
        ncells : int
            Size of supercell

        Returns
        -------
        tuple of ints
            Labeling properly sorted, ready to be given to the enumeration code
        """
        sorted_labeling = [0] * len(self.iter_elements * ncells)
        count = 0
        for _ in range(ncells):
            for iter_element in self.iter_elements:

                # Find the site group corresponding to the iter element
                site_group = self.site_groups[iter_element]

                # Find the right element by checking (1) where the
                # proper site_group is in the unosrted labeling and
                # (2) which element is next in turn
                sorted_labeling[count] = labeling[
                    site_group.position][site_group.next_to_add]
                count += 1
                site_group.next_to_add += 1

        # Reset site group counters
        for site_group in self.site_groups.values():
            assert site_group.next_to_add == len(labeling[site_group.position])
            site_group.next_to_add = 0

        return tuple(sorted_labeling)


class SiteGroup():
    """
    Keeps track of a group of sites that have the same allowed elements.
    I.e. a site group could correspond to all sites on which element 0 and
    1 are allowed, and the number of such sites is stored in the `
    multiplicity`.

    Parameters
    ----------
    iter_element : tuple of ints
        The allowed elements on these sites
    position : int
        Helps to keep track of when the first group occured; the first
        site group encountered will have position = 0, the next 1 etc.
    """

    def __init__(self, iter_element, position):
        self.iter_element = iter_element
        self.position = position
        self.multiplicity = 1
        self.next_to_add = 0  # Will keep count when reordering elements

    def compute_all_combinations(self, ncells):
        """
        Compute all combinations (without considering order) of the elements
        in the group.

        Parameters
        ----------
        ncells : int
            Size of supercell
        """
        self.combinations = []
        for combination in \
            itertools.combinations_with_replacement(self.iter_element,
                                                    self.multiplicity * ncells):
            self.combinations.append(combination)
