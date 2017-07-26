import numpy as np


def get_scaled_positions(positions, cell, wrap=True, pbc=[True, True, True]):
    """Get positions relative to unit cell.

    If wrap is True, atoms outside the unit cell will be wrapped into
    the cell in those directions with periodic boundary conditions
    so that the scaled coordinates are between zero and one."""

    fractional = np.linalg.solve(cell.T, positions.T).T

    if wrap:
        for i, periodic in enumerate(pbc):
            if periodic:
                # Yes, we need to do it twice.
                # See the scaled_positions.py test.
                fractional[:, i] %= 1.0
                fractional[:, i] %= 1.0

    return fractional
