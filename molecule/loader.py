# molecule/loader.py
#----------------------------------------------------------------------------------------------#
#                               Structure loading for Molektra                                 #
#                                                                                              #
#                                  Reads XYZ, PDB and CIF files                                #
#                                                                                              #
# # Copyright (C) 2026 MOLEKTRA                                                                #
# This program is free software: you can redistribute it and/or modify                         #
# it under the terms of the GNU General Public License as published by                         #
# the Free Software Foundation, either version 3 of the License, or                            #
# (at your option) any later version.                                                          #
#----------------------------------------------------------------------------------------------#

from ase.io import read
from ase import Atoms
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

def _read_pdb_tolerant(file_path):
    symbols, positions = [], []
    with open(file_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except ValueError:
                    continue
                elem = line[76:78].strip()
                if not elem:  
                    elem = line[12:16].strip().rstrip('0123456789')[:2]
                symbols.append(elem)
                positions.append([x, y, z])
    if not symbols:
        raise ValueError("No atoms found in PDB")
    return symbols, np.array(positions)


def load_structure(file_path, supercell=(1, 1, 1), do_unwrap=True):
    try:
        atoms = read(file_path)
    except (IndexError, ValueError) as e:
        if file_path.lower().endswith('.pdb'):
            symbols, positions = _read_pdb_tolerant(file_path)
            atoms = Atoms(symbols=symbols, positions=positions)
        else:
            raise
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