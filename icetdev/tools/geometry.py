import numpy as np
from icetdev.lattice_site import LatticeSite


def get_scaled_positions(positions, cell, wrap=True, pbc=[True, True, True]):
    """Get positions relative to unit cell.

    If wrap is True, positions outside the unit cell will be wrapped into
    the cell in those directions with periodic boundary conditions
    so that the scaled coordinates are between zero and one.
    """

    fractional = np.linalg.solve(cell.T, positions.T).T

    if wrap:
        for i, periodic in enumerate(pbc):
            if periodic:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test.
                fractional[:, i] %= 1.0
                fractional[:, i] %= 1.0

    return fractional


def find_lattice_site_from_position_python(structure, position):
    """
    Get lattice neighbor from position.

    This is the Python version of
    `structure.findLatticeSiteFromPosition(position)`

    It is slower but kept for debugging and if further development is needed.
    """

    fractional = np.linalg.solve(structure.cell.T, np.array(position).T).T
    unit_cell_offset = [int(round(x)) for x in fractional]

    residual = np.dot(fractional - unit_cell_offset, structure.cell)
    try:
        index = structure.find_index_of_position(residual)
    except Exception:
        msg = ['error did not find index with pos: {}'.format(residual)]
        msg += ['position in structure are:']
        msg += ['\n' + str(structure.positions)]
        raise Exception(' '.join(msg))

    latNbr = LatticeSite(index, unit_cell_offset)
    return latNbr
