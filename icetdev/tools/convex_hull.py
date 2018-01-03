import numpy as np
from scipy.spatial import ConvexHull as ConvexHullScipy
from scipy.interpolate import griddata
import itertools


class ConvexHull:

    def __init__(self, points, tol=1e-3):
        """
        Construct convex hull for (free) energy of mixing,
        see http://docs.scipy.org/doc/scipy-dev/reference/\
            generated/scipy.spatial.ConvexHull.html

        Parameters
        ----------
        points : 
            concentrations and energies
        retfunc : bool
            if True the method returns a function handle otherwise a dictionary

        Returns
        -------
        function handle
            object that provides a linear interpolation of the convex hull
        dict
            convex hull (concentrations, energies)
        """

        # prepare data in format suitable for ConvexHull
        if isinstance(points[0][0], tuple):
            points = np.array([np.array(point[0] + (point[1],))
                               for point in points])
        else:
            points = np.array([np.array(point) for point in points])

        self.dimensions = len(points[0]) - 1

        # construct convex hull
        hull = ConvexHullScipy(points)

        # collect convex hull points
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
        print(self.energies)

        # Remove points that are above "pure element line"
        self.remove_points_above_pure_elements()
        print(self.energies)

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
            # The points on the convex concentration should always be
            # included, so skip them.
            if i in vertices:
                continue

            # If there are 2 or 3 vertices on the convex concentration hull,
            # this loop will simply yield the same array of vertices If there
            # are more vertices, we need to check the mixed energy against all
            # combinations of three "pure elements" (concentrations on the
            # convex concentration hull).
            for triangle in itertools.combinations(vertices,
                                                   min(len(vertices), 3)):
                # Calculate energy that would be gotten with pure elements
                # with ascribed concentration.
                energy_pure = griddata(self.concentrations[vertices],
                                       self.energies[vertices],
                                       concentration,
                                       method='linear')

                # Prepare to delete if the energy was lowered. griddata gives
                # NaN if the concentration is outside the triangle formed by
                # the three vertices. The result of the below comparison is
                # then False, which is what we want.
                if energy_pure < self.energies[i] - tol:
                    to_delete.append(i)

        # Finally remove all points
        print(to_delete)
        self.concentrations = np.delete(self.concentrations, to_delete, 0)
        self.energies = np.delete(self.energies, to_delete, 0)

    def get_convex_hull_at_concentrations(self, target_concentrations):

        if self.dimensions > 1:
            assert len(target_concentrations[0]) == self.dimensions

        hull_energies = griddata(self.concentrations, self.energies,
                                 target_concentrations, method='linear')
        return hull_energies
