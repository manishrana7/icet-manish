import numpy as np
from scipy.spatial import ConvexHull as ConvexHullScipy
from scipy.interpolate import griddata


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

        # Remove points that are above "pure element line"
        #self._emove_points_above_pure_elements()


    def get_convex_hull_at_concentrations(self, target_concentrations):

        if self.dimensions > 1:
            assert len(target_concentrations[0]) == self.dimensions

        hull_energies = griddata(self.concentrations, self.energies,
                                 target_concentrations, method='linear')
        return hull_energies
