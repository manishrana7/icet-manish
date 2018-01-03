import numpy as np
from scipy.spatial import ConvexHull as ConvexHullScipy
from scipy.interpolate import griddata
import itertools


class ConvexHull(object):
    '''
    Convex hull of concentration and energies (or other quantities).

    Attributes
    ----------
    dimensions : int
        Number of independent concentrations needed to specify a point in
        concentration space (1 for binaries, 2 for ternaries etc)
    concentrations : ndarray of floats, size (N, dimensions)
        Concentration of structures on the convex hull.
    energies : ndarray of floats, size N
        Energy of structures on the convex hull.
    structures : list of ints
        Indices of structures that constitute the convex hull (indices are
        defined by the order their concentrations and energies are feeded when
        initializing the ConvexHull object).
    '''

    def __init__(self, points):
        """
        Construct convex hull for (free) energy of mixing,
        see http://docs.scipy.org/doc/scipy-dev/reference/\
            generated/scipy.spatial.ConvexHull.html

        Parameters
        ----------
        points : list of tuples
            Concentrations and energies as [((c1, c2), e), ...] where (c1, c2)
            are independent concentrations needed to specify a point in
            concentration space, and e is the corresponding energy. In a
            binary system where only one concentration is needed, the format
            [(c, e), ...] is fine.
        """

        # Prepare data in format suitable for Scipy-ConvexHull
        if isinstance(points[0][0], tuple):
            points = np.array([np.array(point[0] + (point[1],))
                               for point in points])
        else:
            points = np.array([np.array(point) for point in points])

        self.dimensions = len(points[0]) - 1

        # Construct convex hull
        hull = ConvexHullScipy(points)

        # Collect convex hull points in handy arrays
        concentrations = []
        energies = []
        for vertex in hull.vertices:
            if self.dimensions == 1:
                concentrations.append(points[vertex][0])
            else:
                concentrations.append(points[vertex][0:-1])
            energies.append(points[vertex][-1])
        concentrations = np.array(concentrations)
        energies = np.array(energies)

        # If there is just one independent concentration, we'd better sort
        # according to it
        if self.dimensions == 1:
            concentrations, energies = (np.array(t) for t in
                                        zip(*sorted(zip(concentrations,
                                                        energies))))

        self.concentrations = concentrations
        self.energies = energies

        # Remove points that are above "pure element line"
        self.remove_points_above_pure_elements()

    def remove_points_above_pure_elements(self, tol=1e-6):
        '''
        Remove all points on the convex hull that correspond to maximum rather
        than minimum energy.

        Parameters
        ----------
        tol : float
            Tolerance for what energy constitutes a lower one.
        '''

        # Identify the "complex concentration hull", i.e. the extremum
        # concentrations. In the simplest case, these should simply be the
        # pure elements.
        if self.dimensions == 1:
            # Then the ConvexHullScipy function doesn't work, so we just pick
            # the indices of the lowest and highest concentrations.
            vertices = []
            vertices.append(np.argmin(self.concentrations))
            vertices.append(np.argmax(self.concentrations))
            vertices = np.array(vertices)
        else:
            concentration_hull = ConvexHullScipy(self.concentrations)
            vertices = concentration_hull.vertices

        # Remove all points of the convex energy hull that have an energy that
        # is higher than what would be gotten with pure elements at the same
        # concentration. These points are mathematically on the convex hull,
        # but in the physically uninteresting upper part, i.e. they maximize
        # rather than minimize energy.
        to_delete = []
        for i, concentration in enumerate(self.concentrations):
            # The points on the convex concentration hull should always be
            # included, so skip them.
            if i in vertices:
                continue

            # The energy obtained as a linear combination of concentrations on
            # the convex hull is "z coordinate" of the position on a
            # (hyper)plane in the (number of independent concentrations +
            # 1)-dimensional (N-D) space. This plane is spanned by N points.
            # If there are more vertices on the convex hull, we need to loop
            # over all combinations of N vertices.
            for plane in itertools.combinations(vertices,
                                                min(len(vertices),
                                                    self.dimensions + 1)):
                # Calculate energy that would be gotten with pure elements
                # with ascribed concentration.
                energy_pure = griddata(self.concentrations[np.array(plane)],
                                       self.energies[np.array(plane)],
                                       concentration,
                                       method='linear')

                # Prepare to delete if the energy was lowered. griddata gives
                # NaN if the concentration is outside the triangle formed by
                # the three vertices. The result of the below comparison is
                # then False, which is what we want.
                if energy_pure < self.energies[i] - tol:
                    to_delete.append(i)
                    break

        # Finally remove all points
        print(to_delete)
        self.concentrations = np.delete(self.concentrations, to_delete, 0)
        self.energies = np.delete(self.energies, to_delete, 0)

    def get_convex_hull_at_concentrations(self, target_concentrations):
        '''
        Get energy of convex hull at specified concentrations.

        Parameters
        ----------
        target_concentrations : ndarray of floats, shape (N, ndim)
            Concentrations at N target points (ndim is the number of
            independent concentrations needed to specify a point in
            concentration space (2 for binaries, 3 for ternaries etc))

        Returns
        -------
        ndarray of floats, size N
            Energies at the specified target_concentrations. If any
            concentration is outside the allowed range, NaN is returned.
        '''
        if self.dimensions > 1:
            assert len(target_concentrations[0]) == self.dimensions

        hull_energies = griddata(self.concentrations, self.energies,
                                 target_concentrations, method='linear')
        return hull_energies
