from ase.build import cut
import numpy as np


def map_configuration_to_reference(input_configuration,
                                   reference_structure,
                                   tolerance_mapping,
                                   vacancy_type=None,
                                   tolerance_cell=0.05,
                                   tolerance_positions=0.01,
                                   verbose=False):
    '''
    Map a relaxed configuration onto a reference configuration.

    Parameters
    ----------
    input_configuration :  ASE Atoms object
        relaxed input configuration
    reference_structure : ASE Atoms object
        reference structure, which can but need not represent the
        primitive cell
    tolerance_mapping : float
        maximum allowed displacement for mapping an atom in the
        relaxed (but rescaled) configuration to the reference
        supercell

        *Note*: A reasonable choice is up to 20-30% of the first
        nearest neighbor distance (`r1`). A value above 50% of `r1`
        will most likely lead to atoms being multiply assigned, which
        will raise an exception.
    vacancy_type : str
        If this parameter is set to a non-zero string unassigned sites
        in the reference structure will be assigned to this type.

        *Note 1*: By default (``None``) the method will print an error
        message and raise an exception if there are *any* unassigned
        sites in the reference structure.

        *Note 2*: `vacancy_type` must be a valid element type as
        enforced by the ASE Atoms class.
    tolerance_cell : float
        tolerance factor applied when computing permutation matrix to
        generate supercell (h = P h_p)
    tolerance_positions : float
        tolerance factor applied when scanning for overlapping
        positions in Angstrom (forwarded to `ase.build.cut`)
    verbose : bool
        turn on verbose output

    Returns
    -------
    ASE Atoms, float, float

        - the ideal supercell most closely matching the input
          configuration
        - the largets deviation of any input coordinate from its
          ideal coordinate
        - the average deviation of the input coordinates from the
          ideal coordinates

    Example
    -------
    The following code snippet illustrates the general usage. To this end, it
    first creates a primitive FCC cell, which is latter used as reference
    structure.  To emulate a relaxed structure obtained from, e.g., a density
    functional theory calculation, the code then creates a 4x4x4 conventional
    FCC supercell, which is populated with two different atom types, has
    distorted cell vectors, and random displacements to the atoms. Finally,
    the present function is used to map the configuration back the ideal
    lattice::

        from ase.build import bulk
        reference = bulk('Au', a=4.09)
        atoms = bulk('Au', cubic=True, a=4.09).repeat(4)
        atoms.set_chemical_symbols(10 * ['Ag'] + (len(atoms) - 10) * ['Au'])
        atoms.set_cell(atoms.cell * 1.02, scale_atoms=True)
        atoms.rattle(0.1)
        mapped_atoms = map_configuration_to_reference(atoms, reference, 1.0)

    '''

    np.set_printoptions(suppress=True, precision=6)

    if np.any(input_configuration.pbc != reference_structure.pbc):
        s = '''The periodic boundary conditions differ
        between input and reference structure'''
        raise Exception(s)

    # The input structure is scaled in order to match the lattice
    # structure of the reference configuration. The reference
    # structure can be a primitive cell, in which case the input
    # configuration would usually be a supercell thereof. This
    # requires some special care here.
    atvol_in = input_configuration.get_volume() / len(input_configuration)
    atvol_ref = reference_structure.get_volume() / len(reference_structure)
    scale = atvol_in / atvol_ref
    if verbose:
        print('Volume per atom ratio of input configuration and reference'
              ' structure: {}'.format(scale))
        print('Number of atoms in reference configuration:'
              ' {}'.format(len(reference_structure)))
        print('Number of atoms in input configuration:'
              ' {}'.format(len(input_configuration)))
        print('Reference cell metric:\n'
              '{}'.format(reference_structure.cell))
        print('Input cell metric:\n'
              '{}'.format(input_configuration.cell))
        print('Reference cell volume per atom:\n'
              ' {}'.format(atvol_ref))
        print('Input cell volume per atom:\n'
              '{}'.format(atvol_in))

    # rescale cell metric of input configuration to match volume
    # per atom of reference structure
    modcell = input_configuration.get_cell()
    if vacancy_type is None:
        modcell *= (1.0 / scale) ** (1.0 / 3.0)
        if verbose:
            print('Modified input cell metric:\n'
                  ' {}'.format(modcell))
            atvol_mod = np.linalg.det(modcell) / len(input_configuration)
            print('Modified input cell volume per atom :'
                  ' {}'.format(atvol_mod))

    # obtain the (in general non-integer) transformation matrix
    # connecting the input configuration to the reference structure
    # L = L_p.P --> P = L_p^-1.L
    P = np.dot(modcell, np.linalg.inv(reference_structure.cell))
    if verbose:
        print('P:\n {}'.format(P))
        print('P_round:\n {}'.format(np.around(P)))
        dP = np.linalg.det(np.around(P))
        print('det(P_round):\n {}'.format(dP))
    # assert that the transformation matrix does not deviate too
    # strongly from the nearest integer matrix
    if np.linalg.norm(P - np.around(P)) / 9 > tolerance_cell:
        s = 'Failed to map configuration to reference'
        s += 'structure (tolerance_cell exceeded).\n'
        s += 'reference:\n {}\n'.format(reference_structure.cell)
        s += 'input:\n {}\n'.format(input_configuration.cell)
        s += 'modcell:\n {}\n'.format(modcell)
        s += 'P:\n {}\n'.format(P)
        s += 'P_round:\n {}\n'.format(np.around(P))
        s += 'Deviation: {}\n'.format(np.linalg.norm(P - np.around(P)) / 9)
        s += 'You can try raising `tolerance_cell`.'
        raise Exception(s)
    # reduce the (real) transformation matrix to the nearest integer one
    P = np.around(P)
    if verbose:
        print('Transformation matrix connecting reference structure'
              ' and idealized input configuration:\n {}'.format(P))
    # scale input configuration to idealized cell metric
    scaled_configuration = input_configuration.copy()
    scaled_configuration.set_cell(np.dot(P, reference_structure.cell),
                                  scale_atoms=True)
    # generate supercell of (presumably primitive) reference structure
    ideal_supercell = cut(reference_structure,
                          P[0], P[1], P[2],
                          tolerance=tolerance_positions)
    if verbose:
        print('ideal_supercell:\n{}'.format(ideal_supercell.cell))
    if (len(ideal_supercell) != len(scaled_configuration) and
            vacancy_type is None):
        s = 'Number of atoms in ideal supercell does not'
        s += ' match input configuration.\n'
        s += 'ideal: {}\n'.format(len(ideal_supercell))
        s += 'input: {}'.format(len(scaled_configuration))
        raise Exception(s)
    if verbose:
        print('Cell metric of input configuration:\n'
              '{}'.format(input_configuration.cell))
        print('Cell metric of rescaled input configuration:\n'
              '{}'.format(scaled_configuration.cell))
        print('Cell metric of ideal supercell:\n'
              '{}'.format(ideal_supercell.cell))

    # map atoms in input configuration to closest site in ideal
    # supercell
    dr_max = 0.0
    dr_sum = 0.0
    dr_sumsq = 0.0
    # per-atom-list for keeping track of mapped atoms
    mapped = [-1] * len(ideal_supercell)
    # distances between ideal and input sites
    drs = [None] * len(ideal_supercell)
    for ideal_site in ideal_supercell:
        if vacancy_type is not None:
            try:
                ideal_site.symbol = vacancy_type
            except:
                s = 'Failed to assign "{}" as vacancy type.\n'
                s += 'Check whether `vacancy_type` represents a'
                s += ' valid element symbol.'.format(vacancy_type)
                raise Exception(s)
        for atom in scaled_configuration:
            # in order to compute the distance the current atom from
            # the input configuration is temporarily added to the
            # ideal supercell. This allows one to simply use the ASE
            # Atoms method for computing the interatomic distance
            ideal_supercell.append(atom)
            dr = ideal_supercell.get_distance(ideal_site.index,
                                              ideal_supercell[-1].index,
                                              mic=True)
            del ideal_supercell[-1]
            if dr < tolerance_mapping:
                if mapped[ideal_site.index] >= 0:
                    s = 'More than one atom from the relaxed'
                    s += ' (and rescaled) configuration have been'
                    s += ' mapped onto the same ideal site.\n'
                    s += ' Try reducing `tolerance_mapping`.'
                    raise Exception(s)
                mapped[ideal_site.index] = atom.index
                drs[ideal_site.index] = dr
                ideal_site.symbol = atom.symbol
                dr_max = max(dr, dr_max)
                dr_sum += dr
                dr_sumsq += dr * dr
                break
        else:
            if vacancy_type is not None:
                continue
            s = 'Failed to assign an atom from the relaxed (and'
            s += ' rescaled) configuration to the ideal lattice.'
            s += ' Try increasing `tolerance_mapping`.\n'
            s += ' {}'.format(ideal_site)
            raise Exception(s)

    dr_avg = dr_sum / len(ideal_supercell)
    dr_sdv = np.sqrt(dr_sumsq / len(ideal_supercell) - dr_avg ** 2)
    if verbose:
        print('Maximum, average and standard deviation of atomic'
              ' displacements: {} {} {}'.format(dr_max, dr_avg, dr_sdv))

        print('{:52} {:}'.format('Input configuration:',
                                 'Scaled configuration:'))
        for k, (input_atom, scaled_atom) in enumerate(
                zip(input_configuration,
                    scaled_configuration)):
            print('{:4}'.format(k), end='')
            print(' {:2}'.format(input_atom.symbol), end='')
            print(' {:12.6f} {:12.6f} {:12.6f}'.format(*input_atom.position),
                  end='')
            print('    -->', end='')
            print(' {:12.6f} {:12.6f} {:12.6f}'.format(*scaled_atom.position),
                  end='')
            print('')

        print('{:52} {}'.format('Ideal supercell:',
                                'Scaled configuration:'))
        for ideal_atom, k, dr in zip(ideal_supercell, mapped, drs):
            print(' {:2}'.format(ideal_atom.symbol), end='')
            print((3 * '  {:12.6f}').format(*ideal_atom.position), end='')
            print('    -->', end='')
            print(' {:4}'.format(k), end='')
            if k >= 0:
                scaled_pos = scaled_configuration[k].position
                print((3 * ' {:12.6f}').format(*scaled_pos), end='')
                print('    --> {:.4}'.format(dr), end='')
            print('')

    for k in set(mapped):
        if k >= 0 and mapped.count(k) > 1:
            s = 'Site {} has been assigned more than once.'.format(k)
            raise Exception(s)

    for element in set(input_configuration.get_chemical_symbols()):
        n1 = input_configuration.get_chemical_symbols().count(element)
        n2 = ideal_supercell.get_chemical_symbols().count(element)
        if n1 != n2:
            s = """Number of atoms of type {} differs between
            input configuration ({}) and ideal
            supercell ({}).\n""".format(element, n1, n2)
            print(s)
            raise Exception(s)

    return ideal_supercell, dr_max, dr_avg
