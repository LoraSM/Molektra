# gui/coordination_geometry.py

# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

SHAPE_SYMMETRIES = {
    'L-2': {'Ci': 'Ci'}



}


SYMMETRY_LABELS = {
    'E': 'Identity (E)',
    'Ci': 'Inversion (Ci)',
    'Cs': 'Reflection (Cs)',
    'C2': 'C2 Rotation',
    'C3': 'C3 Rotation',
    'C4': 'C4 Rotation',
    'C5': 'C5 Rotation',
    'C6': 'C6 Rotation',
    'C8': 'C8 Rotation',
    'C∞': 'Linear Rotation (C∞)',
    'S2': 'S2 Rotation-Reflection',
    'S4': 'S4 Rotation-Reflection',
    'S6': 'S6 Rotation-Reflection',
    'S8': 'S8 Rotation-Reflection',
}

SHAPE_LABELS = {
    'L-2': 'Linear',
    'vT-2': 'V-shape',
    'vOC-2': 'L-shape',
}

POINT_GROUPS_TO_TEST = [
    'C1', 'Ci', 'Cs',
    'C2', 'C2v', 'C2h',
    'C3', 'C3v', 'C3h',
    'C4', 'C4v', 'C4h',
    'C6', 'C6v', 'C6h',
    'D2', 'D2h', 'D2d',
    'D3', 'D3h', 'D3d',
    'D4', 'D4h', 'D4d',
    'D6', 'D6h',
    'T', 'Td', 'Th',
    'O', 'Oh',
    'I', 'Ih',
]