# molecule/loader.py
#----------------------------------------------------------------------------------------------#
#                               Structure loading for Molektra                                 #
#                                                                                              #
#                                  Reads XYZ, PDB and CIF files                                #
#----------------------------------------------------------------------------------------------#
from ase.io import read
from ase.neighborlist import neighbor_list, natural_cutoffs
import numpy as np

SUPPORTED_EXTENSIONS = (".xyz", ".pdb", ".cif")


def _atoms_to_tuples(atoms):
    return [
        (atom.symbol, float(atom.position[0]),
         float(atom.position[1]), float(atom.position[2]))
        for atom in atoms
    ]


def unwrap(atoms, bond_multiplier=1.2, center=True):
    """Unwraps all molecules/fragments in the cell independently.

    Handles disconnected components (complex + counterions + solvent).
    Uses accumulated periodic image offsets rather than repeated MIC
    lookups, so fragments spanning more than one cell length stay intact.
    """
    if atoms.cell.rank == 0 or not any(atoms.pbc):
        return atoms

    mol = atoms.copy()

    radii = natural_cutoffs(mol, mult=bond_multiplier)
    i_ind, j_ind, offsets = neighbor_list('ijS', mol, cutoff=radii)

    neighbors = {i: [] for i in range(len(mol))}
    for i, j, S in zip(i_ind, j_ind, offsets):
        neighbors[i].append((j, S))

    cell = mol.get_cell()
    shifts = np.zeros((len(mol), 3))
    visited = set()

    for seed in range(len(mol)):
        if seed in visited:
            continue
        visited.add(seed)
        queue = [seed]
        while queue:
            current = queue.pop(0)
            for neighbor, S in neighbors[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    shifts[neighbor] = shifts[current] + S @ cell
                    queue.append(neighbor)

    mol.set_positions(mol.get_positions() + shifts)

    if center:
        mol.center()
    return mol


def load_structure(file_path, supercell=(1, 1, 1), do_unwrap=True):
    atoms = read(file_path)

    if len(atoms) == 0:
        raise ValueError(f"No atoms found in: {file_path}")

    is_periodic = atoms.cell.rank == 3 and any(atoms.pbc)

    if is_periodic:
        if tuple(supercell) != (1, 1, 1):
            atoms = atoms.repeat(supercell)
        if do_unwrap:
            atoms = unwrap(atoms)

    cell = atoms.get_cell().array if is_periodic else None
    return _atoms_to_tuples(atoms), cell


def load_xyz(file_path):
    return load_structure(file_path)


def load_pdb(file_path):
    return load_structure(file_path)