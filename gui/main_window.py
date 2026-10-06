# Copyright (C) 2026 MOLEKTRA
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

import numpy as np, os, pyqtgraph as pg , csv
import io, matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from cosymlib import Cosymlib
from PyQt5.QtWidgets import (
     QMainWindow, QWidget, QPushButton,QTableWidget,QMessageBox,QHeaderView,QMenu,QLabel,
     QVBoxLayout, QHBoxLayout, QFileDialog,QTableWidgetItem,QColorDialog,QDialog,QSizePolicy,
     QGroupBox,QShortcut,QLineEdit
)
from PyQt5.QtCore import QCoreApplication, QUrl,Qt
QCoreApplication.processEvents()
from matplotlib.path import Path as MplPath
from PyQt5.QtCore import Qt,QSize,QPoint
from PyQt5.QtGui import QPalette, QColor, QVector4D,QDesktopServices,QPixmap,QIcon,QPainter, QPen, QKeySequence
from cosymlib import Geometry
from pyqtgraph.opengl import GLViewWidget
from pyqtgraph.opengl import GLTextItem
from pyqtgraph.opengl import GLMeshItem, MeshData
from pyqtgraph.opengl import GLLinePlotItem
from matplotlib import colors 
from utils.graphics import create_cylinder
from gui.atom_colors import atom_color_dict
from gui.atom_info_window import AtomInfoWindow
from gui.shape_data import SHAPE_REFERENCES
from gui.coordination_geometry import SHAPE_SYMMETRIES, POINT_GROUPS_TO_TEST
from gui.shape_results_window import ShapeResultsWindow
from gui.measurements_window import MeasurementsWindow
from gui.label_window import LabelWindow
from gui.style_window import StyleWindow
from gui.style_window import StyleWindow2
from posym import SymmetryMolecule
from posym.config import Configuration
from molecule.loader import load_structure
from pointgroup import PointGroup

###For selection 
class SelectionOverlay(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setAttribute(Qt.WA_TransparentForMouseEvents) 
        self.mode = None
        self.path = []

    def update_path(self, mode, path):
        self.mode = mode
        self.path = path
        self.update()  

    def clear(self):
        self.mode = None
        self.path = []
        self.update()

    def paintEvent(self, event):
        if not self.path or self.mode is None:
            return
        painter = QPainter(self)
        pen = QPen(QColor(31, 65, 74), 2, Qt.DashLine)
        painter.setPen(pen)

        pts = self.path
        if self.mode == 'rect':
            x0, y0 = pts[0]
            x1, y1 = pts[-1]
            painter.drawRect(int(min(x0, x1)), int(min(y0, y1)),
                             int(abs(x1 - x0)), int(abs(y1 - y0)))
        elif self.mode == 'free':
            for a, b in zip(pts, pts[1:]):
                painter.drawLine(int(a[0]), int(a[1]), int(b[0]), int(b[1]))


#SYMMETRY/SHAPE
def calculate_symmetry_posym(coords, symbols, groups=None):
    Configuration().scan_steps = 10
    print("N átomos:", len(symbols))
    print("Símbolos:", symbols)
    if groups is None:
        groups = POINT_GROUPS_TO_TEST
    if isinstance(groups, str):
        groups = [groups]

    results = {}
    operations = []
    center = np.zeros(3)
    error_reason = None
    coordinates = [list(c) for c in coords]

    for group in groups:
        try:
            sym = SymmetryMolecule(group=group, coordinates=coordinates, symbols=symbols)
            results[group] = round(float(sym.measure_pos), 4)
            operations = sym.get_oriented_operations()
            center = np.array(sym.center)  
        except Exception as e:
            if 'permutation' in str(e).lower():
                error_reason = 'permutation'
            elif 'axis' in str(e).lower():
                error_reason = 'empty'
            else:
                error_reason = str(e)

    return dict(sorted(results.items(), key=lambda x: x[1])), operations, center

def polyhedron_edges(points, tol=1e-3):
        from collections import defaultdict
        from scipy.spatial import ConvexHull
        n = len(points)
        try:
            hull = ConvexHull(points)
            edge_faces = defaultdict(list)
            for f,(a,b,c) in enumerate(hull.simplices):
                for e in ((a,b),(b,c), (a,c)):
                    edge_faces[tuple(sorted(e))].append(f)
            return [e for e,(f1,f2) in edge_faces.items()
                    if np.dot(hull.equations[f1,:3], hull.equations[f2,:3]) < 1 - tol]
        except Exception:
            print("ConvexHull failed, using distance-based edge detection",e)
            d = lambda i,j: np.linalg.norm(points[i]-points[j])
            thresold = min(d(i,j) for i in range(n) for j in range(i+1,n)) * 1.2
            return [(i,j) for i in range(n) for j in range(i+1,n) if d(i,j)<thresold]



def get_ideal_coords(symbols, coords, shape_code, central_atom=1):
    geo = Geometry(symbols=symbols, positions=np.array(coords))
    mol = Cosymlib([geo])

    buffer = io.StringIO()
    buffer.name = "buffer"
    mol.print_shape_structure([shape_code], central_atom=central_atom, output=buffer)
    output = buffer.getvalue()

    ideal_coords = []
    lines = output.strip().split('\n')
    found = False
    for line in lines:
        if shape_code in line:
            found = True
            continue
        if found:
            parts = line.split()
            if len(parts) == 4:
                ideal_coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return np.array(ideal_coords) if ideal_coords else None

def calculate_distortion_path(symbols, coords, shape_code_A, shape_code_B,
                               central_atom=1, n_points=50):
    geo = Geometry(symbols=symbols, positions=np.array(coords))
    mol_sA = geo.get_shape_measure(shape_code_A, central_atom=central_atom)
    mol_sB = geo.get_shape_measure(shape_code_B, central_atom=central_atom)

    ideal_A = get_ideal_coords(symbols, coords, shape_code_A, central_atom)
    ideal_B = get_ideal_coords(symbols, coords, shape_code_B, central_atom)

    if ideal_A is None or ideal_B is None or len(ideal_A) != len(ideal_B):
        return None, None, mol_sA, mol_sB

    path_sA, path_sB = [], []
    for t in np.linspace(0, 1, n_points):
        interp = (1 - t) * ideal_A + t * ideal_B
        try:
            geo_t = Geometry(symbols=symbols, positions=interp)
            path_sA.append(geo_t.get_shape_measure(shape_code_A, central_atom=central_atom))
            path_sB.append(geo_t.get_shape_measure(shape_code_B, central_atom=central_atom))
        except Exception:
            pass

    return path_sA, path_sB, mol_sA, mol_sB

POINT_GROUP_ORDER = {
    'Ih': 120, 'I': 60,
    'Oh': 48, 'O': 24,
    'Td':24, 'Th': 24, 'T': 12,
    'Dinfh': 999, 'Cinfv': 998, 
    'D8h': 32, 'D8d': 32, 'D8': 16,
    'D7h': 28, 'D7d': 28, 'D7': 14,
    'D6h': 24, 'D6d': 24, 'D6': 12,
    'D5h': 20, 'D5d': 20, 'D5': 10,
    'D4h': 16, 'D4d': 16, 'D4': 8,
    'D3h': 12, 'D3d': 12, 'D3': 6,
    'D2h': 8, 'D2d': 8, 'D2': 4,
    'C8v': 16, 'C8h': 16, 'C8': 8,'S8': 8,
    'C7v': 14, 'C7h': 14, 'C7': 7,
    'C6v': 12, 'C6h': 12, 'C6': 6, 'S6': 6,
    'C5v': 10, 'C5h': 10, 'C5': 5,
    'C4v': 8, 'C4h': 8, 'C4': 4, 'S4': 4,
    'C3v': 6, 'C3h': 6, 'C3': 3,
    'C2v': 4, 'C2h': 4, 'C2': 2,
    'Cs': 2, 'Ci': 2, 'C1': 1
}

def get_group_hierarchy(group_string):
    if not group_string:
        return 0
    clean_group = group_string.strip()
    return POINT_GROUP_ORDER.get(clean_group, 1)

def choose_point_group(results, threshold=0.01):
    if not results:
        return None, None
    min_value = min(results.values())
    candidates = [
        (group, value) for group, value in results.items()
        if value - min_value <= threshold
    ]
    best_group, best_value = max(candidates, key=lambda gv: get_group_hierarchy(gv[0]))
    return best_group, best_value

def get_candidate_groups(atom_positions_3d, atom_types, tolerances_ang=(2,4,8,15)):
    positions = np.array(atom_positions_3d)
    candidates = set()
    for tol_ang in tolerances_ang:
        pg = PointGroup(positions=positions, symbols=atom_types, tolerance_ang=tol_ang)
        candidates.add(pg.get_point_group())
    return candidates


def calculate_point_group(atom_positions_3d, atom_types, threshold=0.01, tolerances_ang=(2,4,8,15)):
    try:
        candidates = get_candidate_groups(atom_positions_3d, atom_types, tolerances_ang)
        if len(candidates) == 1:
            return candidates.pop()

        coordinates = [list(c) for c in atom_positions_3d]
        symbols = atom_types
        results, operations, center = calculate_symmetry_posym(
            coordinates, symbols, groups=list(candidates)
        )
        if not results:
            return None
        best_group, best_value = choose_point_group(results, threshold=threshold)
        return best_group
    except Exception as e:
        print(f"Error calculating point group: {e}")
        return None
      
def _symmetry_eigenvector(matrix, target_eigenvalue):
    eigenvalues, eigenvectors = np.linalg.eig(matrix)
    idx = np.argmin(np.abs(np.real(eigenvalues) - target_eigenvalue))
    return np.real(eigenvectors[:, idx])

def get_symmetry_axis(matrix):
    # rotation axis: eigenvector whose eigenvalue is +1
    return _symmetry_eigenvector(matrix, 1.0)

def get_plane_normal(matrix):
    # mirror-plane normal: eigenvector whose eigenvalue is -1
    return _symmetry_eigenvector(matrix, -1.0)

def calc_distance(p1, p2):
    return np.linalg.norm(p2 - p1)

def calc_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    cos_a = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    cos_a = np.clip(cos_a, -1.0, 1.0)
    return np.degrees(np.arccos(cos_a))

def calc_dihedral(p1, p2, p3, p4):
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3
    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)
    n1 /= np.linalg.norm(n1)
    n2 /= np.linalg.norm(n2)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    x = np.dot(n1, n2)
    y = np.dot(m1, n2)
    return np.degrees(np.arctan2(y, x))

def create_dashed_line(p1, p2, n_dashes=20, color=(1, 1, 0, 1)):
    points = []
    for i in range(n_dashes):
        t_start = i / n_dashes
        t_end   = (i + 0.6) / n_dashes  
        points.append(p1 + t_start * (p2 - p1))
        points.append(p1 + t_end   * (p2 - p1))
        points.append([np.nan, np.nan, np.nan])  
    return GLLinePlotItem(
        pos=np.array(points),
        color=color,
        width=5,
        antialias=True,
        mode='lines'
    )

def calculate_shape(central_atom_idx, atom_positions_3d, atom_types, ligand_indices):
    central_pos = atom_positions_3d[central_atom_idx]
    central_sym = atom_types[central_atom_idx]

    print(f"Detected ligands: {len(ligand_indices)}")
    cn = len(ligand_indices)
    coords = [central_pos] + [atom_positions_3d[i] for i in ligand_indices]
    symbols = [central_sym] + [atom_types[i] for i in ligand_indices]

    #-----------GEOMETRY//SHAPE MEASURES----------------
    geo = Geometry(symbols=symbols, positions=coords)
    results = {}
    shapes = SHAPE_REFERENCES.get(cn, {})

    for shape_code, shape_name in shapes.items():
        value = geo.get_shape_measure(shape_code, central_atom=1)
        results[shape_name] = value
        print(f"{shape_name}: {value:.4f}")

    return results,cn

#----------------------------------------------------------------------#

#FUNCTIONS FOR 3D VISUALIZATION
def create_sphere(radius=0.3, color=(1, 0, 0, 1)):
    md = MeshData.sphere(rows=20, cols=40, radius=radius)
    sphere = GLMeshItem(meshdata=md, smooth=True, color=color, shader='shaded', glOptions='opaque')
    return sphere

def covalent_threshold(type1, type2):
    r1 = atom_color_dict.get(type1, {}).get('cov_radius', 1.5)
    r2 = atom_color_dict.get(type2, {}).get('cov_radius', 1.5)
    return (r1 + r2) +0.3

def get_ligands(central_idx, positions, types, extra_bonds=None, excluded_bonds=None):
    central_pos = np.asarray(positions[central_idx])
    central_sym = types[central_idx]
    ligands = [
        i for i, (pos,sym) in enumerate(zip(positions, types))
        if i != central_idx and np.linalg.norm(np.asarray(pos) - central_pos) <= covalent_threshold(central_sym, sym)
    ]
    if extra_bonds:
        for a, b in extra_bonds:
            if a == central_idx and b not in ligands:
                ligands.append(b)
            elif b == central_idx and a not in ligands:
                ligands.append(a)
    if excluded_bonds:
        ligands = [
            i for i in ligands
            if (central_idx, i) not in excluded_bonds
            and (i, central_idx) not in excluded_bonds
        ]
    return ligands




class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self.export_data = {}
        self._measure_registry = {}
        self._measure_counter = 0
        self._bond_pairs = []
        self.manual_bonds = []
        self.excluded_bonds = []
        self._label_items = []
        self._labeled_atoms = []
        self.setAcceptDrops(True)
        self.style_window = StyleWindow(parent=self) 
        self.style_window_bond = StyleWindow2(parent=self)
        self.bond_color = (0.3, 0.3, 0.3, 1)
        self.bond_radius = 0.1
        self.label_window = LabelWindow(parent=self)
        self.measurements_window = MeasurementsWindow(parent=self)
        self.isolate_mode = False
        self.setWindowTitle("Molektra")
        self.setStyleSheet("""
            QMainWindow, QWidget {
                background-color: #f8f9fa;
            }

            QMessageBox {background-color: #ffffff;}
            QMessageBox QLabel {color:black;}

            QPushButton {
                background-color: #ffffff;
                color: #1f2937;
                font-size: 13px;
                font-weight: 600;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 9px 14px;
            }
            QPushButton:hover {
                background-color: #f3f4f6;
                border-color: #d1d5db;
            }
            QPushButton:pressed {
                background-color: #e5e7eb;
            }

            /* Acciones principales, por objectName */
            QPushButton#primary {
                background-color: #2563eb;
                color: #ffffff;
                border: none;
            }
            QPushButton#primary:hover   { background-color: #1d4ed8; }

            QPushButton#success {
                background-color: #059669;
                color: #ffffff;
                border: none;
            }
            QPushButton#success:hover   { background-color: #047857; }

            QPushButton#accent {
                background-color: #7c3aed;
                color: #ffffff;
                border: none;
            }
            QPushButton#accent:hover    { background-color: #6d28d9; }

            QPushButton#danger {
                background-color: #ffffff;
                color: #b91c1c;
                border: 1px solid #fecaca;
            }
            QPushButton#danger:hover    { background-color: #fef2f2; }
        """)
        self.resize(600, 600)
        self.shape_results_window = ShapeResultsWindow(parent=self)
        self.shape_results_window.table.itemSelectionChanged.connect(self.on_shape_selection_changed)
        self.shape_results_window.show_poly_button.clicked.connect(self.toggle_polyhedron)
        self.shape_results_window.distortion_path_button.clicked.connect(self.show_distortion_path)  
        self._polyhedron_items = []

        self.atom_positions_3d = []  
        self.atom_types = []          
        self.selected_atom_index = None  
        self.selected_atoms = []
        self.selection_mode = None
        self._drag_path = []
        
        self.atom_info_window = AtomInfoWindow(parent=self)

        #---------------Layout---------------------#
        container = QWidget()
        main_layout = QHBoxLayout()
        container.setLayout(main_layout)
        self.setCentralWidget(container)

        #---------------3D View ---------------------#
        self.molecule_view = GLViewWidget()
        self.molecule_view.setMinimumSize(400, 600)
        self.molecule_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.molecule_view.opts['distance'] = 20
        self.molecule_view.setBackgroundColor('#f0f0f0')
        palette = self.molecule_view.palette()
        palette.setColor(QPalette.Window, QColor("#e0e0e0"))
        self.molecule_view.setPalette(palette)

        self.molecule_view.setContextMenuPolicy(Qt.CustomContextMenu)
        self.molecule_view.customContextMenuRequested.connect(self.show_context_menu)
        self.molecule_view.mouseMoveEvent = self.on_mouse_move
        self.molecule_view.mousePressEvent = self.on_mouse_press
        self.molecule_view.mouseReleaseEvent = self.on_mouse_release
        self.selection_overlay = SelectionOverlay(self.molecule_view)
        self.selection_overlay.setGeometry(self.molecule_view.rect())
        self.selection_overlay.hide()
        self.last_mouse_pos = None
        self.sym_label_overlay = QLabel("", self.molecule_view)
        self.sym_label_overlay.setStyleSheet("""
        QLabel {
            color: white;
            background-color: rgba(0, 0, 0, 140);
            border-radius: 6px;
            padding: 4px 10px;
            font-size: 14px;
            font-weight: bold;
            }
            """)
        self.sym_label_overlay.move(10, 10)
        self.sym_label_overlay.hide()


        #--------Visualization Toolbar----------------#
        toolbar = QWidget()
        toolbar.setFixedHeight(44)
        toolbar.setStyleSheet('background-color: #FFFFFF;color:#6E6E6E; border-bottom: 1px solid #e5e7eb')
        toolbar_layout = QHBoxLayout()
        toolbar_layout.setContentsMargins(8, 5, 8, 5)
        toolbar_layout.setSpacing(6)
        toolbar.setLayout(toolbar_layout)

        self.delete_button = self.make_toolbar_button("delete.svg", "Delete selected atoms", self.delete_selection, icon_size=20)
        self.boton_ejes     = self.make_toolbar_button("axis.svg",    "Show axes XYZ", self.toggle_ejes_xyz)
        self.hide_h_button  = self.make_toolbar_button("hydrogen.png","Show/Hide Hydrogens", self.toggle_hydrogens,icon_size=25)
        self.coord_env_button = self.make_toolbar_button("cn.png", "Isolate coordination environment", self.toggle_coord_env, icon_size=25)
        self.selection_button = self.make_toolbar_button("selection.png", "Select atoms", self.open_selection_menu, icon_size=25)
        self.reset_button = self.make_toolbar_button("reset.svg", "Reset structure", self.restore_original, icon_size=26)
        self.distance_button = self.make_toolbar_button("distance.png", "Measure distance", self.measure_distance, icon_size=25)
        self.angle_button = self.make_toolbar_button("angle.png", "Measure angle", self.measure_angle, icon_size=40)
        self.dihedral_button = self.make_toolbar_button("dihedral.png", "Measure dihedral angle", self.measure_dihedral, icon_size=28)
        self.add_bond_button = self.make_toolbar_button("add_bond.png", "Add bond (select 2 atoms)", self.add_bond, icon_size=25)
        self.delete_bond_button = self.make_toolbar_button("delete_bond.png", "Delete bond (select 2 atoms)", self.delete_bond, icon_size=25)

        
        toolbar_layout.addWidget(self.hide_h_button)
        toolbar_layout.addWidget(self.distance_button)
        toolbar_layout.addWidget(self.angle_button)
        toolbar_layout.addWidget(self.dihedral_button)
        toolbar_layout.addWidget(self.coord_env_button)
        toolbar_layout.addWidget(self.selection_button)
        toolbar_layout.addWidget(self.boton_ejes) 
        toolbar_layout.addWidget(self.delete_button)  
        toolbar_layout.addWidget(self.reset_button)
        toolbar_layout.addWidget(self.add_bond_button)
        toolbar_layout.addWidget(self.delete_bond_button)
        toolbar_layout.addStretch()

        viewer_container = QWidget()
        viewer_layout = QVBoxLayout()
        viewer_layout.setContentsMargins(0, 0, 0, 0)
        viewer_layout.setSpacing(0)
        viewer_container.setLayout(viewer_layout)
        viewer_layout.addWidget(toolbar)
        viewer_layout.addWidget(self.molecule_view, stretch=1)

        #----------------Shortcuts----------------#
        self._shortcuts = [
        QShortcut(QKeySequence("h"), self, self.toggle_hydrogens),
        QShortcut(QKeySequence("d"), self, self.measure_distance),
        QShortcut(QKeySequence("a"), self, self.measure_angle),
        QShortcut(QKeySequence("x"), self, self.toggle_ejes_xyz),
        QShortcut(QKeySequence("Ctrl+z"), self, self.restore_original),
        QShortcut(QKeySequence("c"), self, self.toggle_coord_env),
        QShortcut(QKeySequence("B"), self, self.add_bond),
        QShortcut(QKeySequence("Shift+B"), self, self.delete_bond),
        QShortcut(QKeySequence("Delete"), self, self.delete_selection),
        QShortcut(QKeySequence("Ctrl+E"),self, self.expand_selection)
        ]

        #---------------Main Toolbar---------------------#
        button_panel = QWidget()
        button_layout = QVBoxLayout()
        button_panel.setFixedWidth(240)
        button_layout.setSpacing(12)
        button_panel.setLayout(button_layout)
        group_style = """
            QGroupBox {
                font-weight: 600;
                font-size: 13px;
                color: #374151;
                border: 1px solid #e5e7eb;
                border-radius: 10px;
                margin-top: 10px;
                padding: 12px 10px 10px 10px;
                background-color: #ffffff;
            }
            QGroupBox::title {
                subcontrol-origin: margin;
                left: 12px;
                padding: 0 4px;
            }
        """
        #---------------GENERAL---------------------#
        general_box = QGroupBox("General")
        general_box.setStyleSheet(group_style)
        g = QVBoxLayout(general_box)
        g.setSpacing(6)

        load_button = QPushButton("Load File")
        load_button.setFixedHeight(35)
        load_button.clicked.connect(self.load_file)
        g.addWidget(load_button)

        #---------------BOTON ATOM---------------------#
        atom_info_button = QPushButton("Atom Information")
        atom_info_button.setFixedHeight(35)
        atom_info_button.clicked.connect(self.show_atom_info)
        g.addWidget(atom_info_button)

                #BOTON RESULTS
        self.export_button = QPushButton("Export")
        self.export_button.setFixedHeight(35)
        self.export_button.clicked.connect(self.toggle_export)
        g.addWidget(self.export_button)

        #BOTON Info 
        self.info_us = QPushButton("About us")
        self.info_us.setFixedHeight(35)
        self.info_us.clicked.connect(self.show_info)
        g.addWidget(self.info_us)
        g.addStretch()

        button_layout.addWidget(general_box)

        #---------------Analysis---------------------#
        analysis_box = QGroupBox("Analysis")
        analysis_box.setStyleSheet(group_style)
        a = QVBoxLayout(analysis_box)
        a.setSpacing(6)

        calc_shape_button = QPushButton("Shape measures")
        calc_shape_button.setFixedHeight(35)
        calc_shape_button.clicked.connect(self.run_shape_calculation)
        calc_shape_button.setStyleSheet("background-color: #D0F5D1; color: black; font-weight: bold;")
        a.addWidget(calc_shape_button)

        #---------------BOTON Point---------------------#
        calc_pointgroup = QPushButton("Get Point Group")
        calc_pointgroup.setFixedHeight(35)
        calc_pointgroup.clicked.connect(self.run_point_group)
        calc_pointgroup.setStyleSheet("background-color: #F5D0F5; color: black; font-weight: bold;")
        a.addWidget(calc_pointgroup)

        #------------------------TABLE----------------------------#
        # Pointgroup #
        self.pointgroup_table = QTableWidget()
        self.pointgroup_table.setColumnCount(1)
        self.pointgroup_table.setHorizontalHeaderLabels(['Point Group'])
        self.pointgroup_table.setMaximumHeight(140)
        self.pointgroup_table.verticalHeader().setVisible(False)
        self.pointgroup_table.setShowGrid(False)
        self.pointgroup_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.pointgroup_table.setAlternatingRowColors(False)
        self.pointgroup_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        self.pointgroup_table.setStyleSheet("""
            QTableWidget {
            background-color: #ffffff;
            border: 1px solid #e5e7eb;
            border-radius: 8px;
            font-size: 13px;
            outline: none;
            }
            QTableWidget::item {
            padding: 8px 14px;
            border: none;
            color: #1f2937;
            border-bottom: 1px solid #f3f4f6;
            }
            QTableWidget::item:selected {
            background-color: #f0f4ff;
            color: #1f2937;
            }
        QHeaderView::section {
            background-color: #f9fafb;
            color: #6b7280;
            font-weight: 600;
            font-size: 12px;
            padding: 7px 14px;
            border: none;
            border-bottom: 1px solid #e5e7eb;
        }
        """)
        self.pointgroup_table.setVisible(False)
        a.addWidget(self.pointgroup_table)

        self.close_pointgroup_button = QPushButton("Close")
        self.close_pointgroup_button.clicked.connect(lambda: (
        self.pointgroup_table.setVisible(False),
        self.close_pointgroup_button.setVisible(False)
        ))
        self.close_pointgroup_button.setVisible(False)
        a.addWidget(self.close_pointgroup_button)



        #-------------Symmetry Measures--------------------#
        calc_sym_posym = QPushButton("Symmetry Measures")
        calc_sym_posym.setFixedHeight(35)
        calc_sym_posym.clicked.connect(self.run_symmetry_posym)
        calc_sym_posym.setStyleSheet("background-color: #D0EAF5; color: black; font-weight: bold;")
        a.addWidget(calc_sym_posym)


        self.btn_reset = QPushButton("Clear")
        self.btn_reset.setFixedHeight(35)
        self.btn_reset.clicked.connect(self.reset_measurements)
        self.btn_reset.setVisible(False)
        a.addWidget(self.btn_reset)
        button_layout.addWidget(analysis_box)


        #----------------Style-----------------#
        style_box = QGroupBox("Style")
        style_box.setStyleSheet(group_style)
        s = QVBoxLayout(style_box)
        s.setSpacing(6)

    
        labels_button = QPushButton("Labels")
        labels_button.setFixedHeight(35)
        labels_button.clicked.connect(self.open_label_window)
        s.addWidget(labels_button)

        style_atoms_button = QPushButton("Atom")
        style_atoms_button.setFixedHeight(35)
        style_atoms_button.clicked.connect(self.open_style_window)
        s.addWidget(style_atoms_button)


        style_bonds_button = QPushButton("Bond")
        style_bonds_button.setFixedHeight(35)
        style_bonds_button.clicked.connect(self.open_style_window_bonds)
        s.addWidget(style_bonds_button)

        button_layout.addWidget(style_box)
        button_layout.addStretch()


        #BOTON EXIT
        self.exit_button = QPushButton("Exit")
        self.exit_button.setFixedHeight(35)
        self.exit_button.clicked.connect(self.close)
        self.exit_button.setStyleSheet("background-color: #F75757; color: black; font-weight: bold;")
        button_layout.addWidget(self.exit_button)

        main_layout.addWidget(viewer_container, stretch=1)
        main_layout.addWidget(button_panel, stretch=0)

        #AXIS
        self.ejes_visibles = False
        self.eje_x, self.eje_y, self.eje_z = self.crear_ejes_xyz()
        


    # MOUSE

    def on_mouse_press(self, event):
        if self.selection_mode in ('rect', 'free'):
            self._drag_path = [np.array([event.pos().x(), event.pos().y()])]
            self.selection_overlay.setGeometry(self.molecule_view.rect())
            self.selection_overlay.show()
            return
        
        if event.button() == Qt.RightButton:
            return

        self._last_mouse_pos = event.pos()
        modifiers = event.modifiers()
        is_ctrl_pressed = modifiers == Qt.ControlModifier
        self.select_atom_closest(event.pos(), ctrl_pressed = is_ctrl_pressed)
        GLViewWidget.mousePressEvent(self.molecule_view, event)

    def on_mouse_release(self, event):
        if self.selection_mode in ('rect', 'free') and self._drag_path:
            self._finish_area_selection()
            self.selection_overlay.clear()
            self.selection_overlay.hide()
            self.selection_mode = None
            self._drag_path = []
            return
        self._last_mouse_pos = None
    
    def on_mouse_move(self, event):
        if self.selection_mode in ('rect', 'free') and self._drag_path:
            self._drag_path.append(np.array([event.pos().x(), event.pos().y()]))
            self.selection_overlay.update_path(self.selection_mode, self._drag_path)
            return
        GLViewWidget.mouseMoveEvent(self.molecule_view, event)


    #------------TOOLBAR----------#
    def make_toolbar_button(self, icon_name, tooltip, callback, icon_size=26):
        btn = QPushButton()
        icon_path = os.path.join(os.path.dirname(__file__), icon_name)
        btn.setIcon(QIcon(icon_path))
        btn.setIconSize(QSize(icon_size, icon_size))
        btn.setToolTip(tooltip)
        btn.setFixedSize(34, 34)
        btn.setStyleSheet("""
            QPushButton {
                background-color: transparent;
                border: none;
                border-radius: 6px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
            """)
        btn.clicked.connect(callback)
        return btn
    
    def _icon(self, name):
        return QIcon(os.path.join(os.path.dirname(__file__), name))


    # SCREEN


    def project_to_screen(self, pos_3d):
        w = self.molecule_view.width()
        h = self.molecule_view.height()

        proj = self.molecule_view.projectionMatrix()
        mv   = self.molecule_view.viewMatrix()
        mvp  = proj * mv  

        v    = QVector4D(float(pos_3d[0]), float(pos_3d[1]), float(pos_3d[2]), 1.0)
        clip = mvp.map(v)

        if clip.w() == 0:
            return None 

        ndc_x = clip.x() / clip.w()
        ndc_y = clip.y() / clip.w()

        screen_x = (ndc_x + 1.0) * 0.5 * w
        screen_y = (1.0 - ndc_y) * 0.5 * h

        return np.array([screen_x, screen_y])


    def show_context_menu(self,pos):
        if not self.selected_atoms:
            return
        menu = QMenu(self)
        menu.setStyleSheet("""
                QMenu {
                background-color: #ffffff; color: #1f2937;
                border: 1px solid #e5e7eb; border-radius: 6px;
            }
            QMenu::item { padding: 6px 20px; }
            QMenu::item:selected { background-color: #f0f4ff; color: #1f2937; }
                           """)
        act_isolate = menu.addAction("Isolate section")
        act_delete  = menu.addAction("Delete atoms selected")
        menu.addSeparator()
        act_shape   = menu.addAction("Shape measures")

        global_pos = self.molecule_view.mapToGlobal(pos)
        chosen = menu.exec_(global_pos)
        if chosen   == act_isolate:
            self.isolate_selection()
        elif chosen == act_delete:
            self.delete_selection()
        elif chosen == act_shape:
            self.run_shape_calculation()


    # ATOM 

    def select_atom_closest(self, pos, threshold=20, ctrl_pressed=False):
        if not self.atom_positions_3d:
            return  

        click = np.array([pos.x(), pos.y()])
        min_dist = float('inf')
        closest_idx = None

        for i, atom_pos in enumerate(self.atom_positions_3d):
            if not self.spheres[i].visible():
                continue
            screen_pos = self.project_to_screen(atom_pos)
            if screen_pos is None:
                continue
            dist = np.linalg.norm(screen_pos - click)
            if dist < min_dist:
                min_dist = dist
                closest_idx = i

        if closest_idx is not None and min_dist < threshold:
            if ctrl_pressed:
                if closest_idx in self.selected_atoms:
                    self.selected_atoms.remove(closest_idx)
                else:
                    self.selected_atoms.append(closest_idx)
            else:
                self.selected_atoms = [closest_idx]
                print(f"Selected atom: {self.atom_types[closest_idx]} (index {closest_idx})")
            self.resaltar_atomo(closest_idx)

        if (hasattr(self, 'shape_results_window')
            and self.shape_results_window.isVisible()
            and len(self.selected_atoms)==1):
            self.run_shape_calculation()



    def _finish_area_selection(self):
        if len(self._drag_path) < 2:
            return
        path = np.array(self._drag_path)

        if self.selection_mode == 'rect':
            x0, y0 = path[0]
            x1, y1 = path[-1]
            xmin, xmax = min(x0, x1), max(x0, x1)
            ymin, ymax = min(y0, y1), max(y0, y1)
            def inside(sx, sy):
                return xmin <= sx <= xmax and ymin <= sy <= ymax
        else:

            poly = MplPath(path)
            def inside(sx, sy):
                return poly.contains_point((sx, sy))

        selected = []
        for i, atom_pos in enumerate(self.atom_positions_3d):
            screen = self.project_to_screen(atom_pos)
            if screen is not None and inside(screen[0], screen[1]):
                selected.append(i)

        self.selected_atoms = selected
        for idx in selected:
            self.resaltar_atomo(idx)
        print(f"Area selection: {len(selected)} atoms")
    def expand_selection(self):
        if not self.selected_atoms or not self.atom_positions_3d:
            return
        new_atoms = []
        for idx in list(self.selected_atoms):
            for neighbor in self._ligands(idx):
                if neighbor not in self.selected_atoms and neighbor not in new_atoms:
                    new_atoms.append(neighbor)
        self.selected_atoms.extend(new_atoms)
        self.resaltar_atomo(self.selected_atoms[0])
            

    # SHAPE

    def run_shape_calculation(self):
        if not self.selected_atoms:
            print("No atoms selected")
            return

        shape_results = {}
        cn = len(self.selected_atoms)

        if len(self.selected_atoms) == 1:
            central_idx = self.selected_atoms[0]
            ligand_indices = self._ligands(central_idx)
            central_pos = self.atom_positions_3d[central_idx]
            central_sym = self.atom_types[central_idx]
            self._last_shape_coords = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
            self._last_shape_symbols = [central_sym] + [self.atom_types[i] for i in ligand_indices]
            self._last_shape_central = 1
            self._last_shape_cn = len(ligand_indices)


            shape_results, cn = calculate_shape(
                central_idx,
                self.atom_positions_3d,
                self.atom_types,
                ligand_indices
            )
        else:
            coords = [self.atom_positions_3d[i] for i in self.selected_atoms]
            symbols = [self.atom_types[i] for i in self.selected_atoms]
            self._last_shape_coords = coords
            self._last_shape_symbols = symbols
            self._last_shape_central = 0
            self._last_shape_cn = cn

            print(f"Calculating shape for {cn} selected atoms: {symbols}")

            try:
                geo = Geometry(symbols=symbols, positions=np.array(coords))
                shape_results = {}
                shapes = SHAPE_REFERENCES.get(cn, {})
                for shape_code, shape_name in shapes.items():
                    value = geo.get_shape_measure(shape_code, central_atom=0)
                    shape_results[shape_name] = value
                    print(f"{shape_name}: {value:.4f}")
            except Exception as e:
                print(f"Error calculating shapes for {cn} selected atoms: {e}")
                shape_results = {}
        if shape_results:
            atom_label = f"Shape_{self.atom_types[self.selected_atoms[0]]}{self.selected_atoms[0]}"
            self.export_data[atom_label] = shape_results

        if not shape_results:
            QMessageBox.information(self, "Shape Measures", f"CN={cn} — no reference shapes defined")
        else:
            atom_symbol = self.atom_types[self.selected_atoms[0]]
            atom_idx    = self.selected_atoms[0]
            self.update_shape_table(shape_results, atom_symbol, atom_idx, cn)


    def load_custom_polyhedron(self):
        if not hasattr(self, '_last_shape_coords'):
            QMessageBox.warning(self, "Custom shape", "Run a shape calculation first.")
            return

        file_path, _ = QFileDialog.getOpenFileName(
            self, "Load custom polyhedron (XYZ)", "", "XYZ Files (*.xyz)")
        if not file_path:
            return
        
        try:
            atoms_tuples, _ = load_structure(file_path)
        except Exception as e:
            QMessageBox.warning(self, "Custom shape", f"Could not read file:\n{e}")
            return

        ref_symbols = [t[0] for t in atoms_tuples]
        ref_coords = [[t[1], t[2], t[3]] for t in atoms_tuples]


        n_expected = len(self._last_shape_coords)
        if len(ref_symbols) != n_expected:
            QMessageBox.warning(self, "Custom shape",
                f"The reference has {len(ref_symbols)} atoms but your fragment has "
                f"{n_expected}. They must match (central + vertices).")
            return

        # CSM custom
        try:
            from cosymlib import Geometry
            from cosymlib.shape import Shape
            problem = Geometry(symbols=self._last_shape_symbols,
                               positions=np.array(self._last_shape_coords))
            
            ref_symbols_reordered = ref_symbols[1:] + [ref_symbols[0]]
            ref_coords_reordered  = ref_coords[1:]  + [ref_coords[0]]

            custom_ref = Geometry(symbols=ref_symbols_reordered, positions=np.array(ref_coords_reordered))
            shape = Shape(problem)
            value = shape.measure(custom_ref, central_atom=1)
        except Exception as e:
            QMessageBox.warning(self, "Custom shape", f"Could not compute measure:\n{e}")
            return


        ref_name = f"Custom ({os.path.basename(file_path)})"

        self.shape_results_window.add_result_row(ref_name, value)
        atom_symbol = self.atom_types[self.selected_atoms[0]]
        atom_idx = self.selected_atoms[0]
        self.export_data.setdefault(f"Shape_{atom_symbol}{atom_idx}", {})[ref_name] = value




    def show_distortion_path(self):
        if not self.selected_atoms:
            QMessageBox.warning(self, "No selection", "Please select an atom first.")
            return
        
        table = self.shape_results_window.table
        selected_rows = list({item.row() for item in table.selectedItems()})
        if len(selected_rows) != 2:
            QMessageBox.warning(self, "Selection", 
                "Select exactly 2 shapes in the table (Ctrl+click).")
            return

        name_A = table.item(selected_rows[0], 0).text()
        name_B = table.item(selected_rows[1], 0).text()

        if len(self.selected_atoms) == 1:
            central_idx = self.selected_atoms[0]
            central_pos = self.atom_positions_3d[central_idx]
            central_sym = self.atom_types[central_idx]
            ligand_indices = self._ligands(central_idx)
            coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
            symbols = [central_sym] + [self.atom_types[i] for i in ligand_indices]
            cn = len(ligand_indices)
            central_atom_param = 1
        else:
            coords  = [self.atom_positions_3d[i] for i in self.selected_atoms]
            symbols = [self.atom_types[i] for i in self.selected_atoms]
            cn = len(self.selected_atoms)
            central_atom_param = 0

        shapes = SHAPE_REFERENCES.get(cn, {})
        code_A = next((c for c, n in shapes.items() if n == name_A), None)
        code_B = next((c for c, n in shapes.items() if n == name_B), None)

        if not code_A or not code_B:
            QMessageBox.warning(self, "Error", "No shape codes were found.")
            return

        try:
            path_sA, path_sB, mol_sA, mol_sB = calculate_distortion_path(
                symbols, coords, code_A, code_B,
                central_atom=central_atom_param
            )
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error calculating distortion path: {e}")
            return

        if path_sA is None:
            QMessageBox.warning(self, "Error", 
                "Ideal structures could not be obtained.")
            return

        # Plot
        dialog = QDialog(self)
        dialog.setWindowTitle(f"Distortion Path: {name_A} ↔ {name_B}")
        dialog.resize(700, 700)
        layout = QVBoxLayout(dialog)

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.grid(False)
        ax.set_facecolor('white')
        fig.patch.set_facecolor('white')

        ax.plot(path_sA, path_sB,
                color='#2563eb', linewidth=1.5, label='Min. distortion path')

        # Real molecule point
        ax.scatter([mol_sA], [mol_sB],
                   color='#dc2626', s=60, zorder=5,
                   label=f'Molecule  ({mol_sA:.2f}, {mol_sB:.2f})')

        # Labels
        lim = max(max(path_sA), max(path_sB), mol_sA, mol_sB)*1.05
        ax.set_xlim(0,lim)
        ax.set_ylim(0,lim)

        ax.set_aspect('equal', adjustable='box')

        ax.annotate(name_A, xy=(path_sA[0], path_sB[0]),
                    xytext=(6, 6), textcoords='offset points',
                    fontsize=9, color='#1d4ed8', fontweight='bold')
        ax.annotate(name_B, xy=(path_sA[-1], path_sB[-1]),
                    xytext=(6, -14), textcoords='offset points',
                    fontsize=9, color='#1d4ed8', fontweight='bold')

        ax.set_xlabel(f'S({name_A})', fontsize=11)
        ax.set_ylabel(f'S({name_B})', fontsize=11)
        ax.set_title('Minimum Distortion Path', fontsize=12, fontweight='bold')
        ax.legend(fontsize=9)

        fig.tight_layout()

        canvas = FigureCanvas(fig)
        layout.addWidget(canvas)

        button_row = QHBoxLayout()

        close_btn = QPushButton("Close")
        close_btn.setFixedSize(100, 34) 
        close_btn.clicked.connect(dialog.close)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #FBE1BC;
                color: black;
                font-size: 13px;
                font-weight: 600;
                border: 1px solid #e5e7eb;
                border-radius: 8px;
                padding: 7px 12px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
        """)

        export_btn = QPushButton("Results ")
        export_btn.setFixedSize(120, 34)
        export_btn.clicked.connect(lambda: self._show_path_export_menu(export_btn, fig, path_sA, path_sB, mol_sA, mol_sB, name_A, name_B))

        button_row.addWidget(close_btn)
        button_row.addStretch()                  
        button_row.addWidget(export_btn)

        layout.addLayout(button_row)

        
        dialog.show()

    def _show_path_export_menu(self, button, fig, path_x, path_y, mol_x, mol_y, name_A, name_B):
        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu { background-color: #ffffff; color: #1f2937;
                    border: 1px solid #e5e7eb; border-radius: 6px; }
            QMenu::item { padding: 6px 20px; }
            QMenu::item:selected { background-color: #f0f4ff; color: #1f2937; }
        """)
        act_img = menu.addAction("Export image")
        act_csv = menu.addAction("Export distortion path (CSV)")
        pos = button.mapToGlobal(button.rect().bottomLeft())
        chosen = menu.exec_(pos)

        if chosen == act_img:
            path, _ = QFileDialog.getSaveFileName(self, "Save image", "distortion_path.png",
                                                  "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
            if path:
                fig.savefig(path, dpi=300, bbox_inches='tight')

        elif chosen == act_csv:
            path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "distortion_path.csv",
                                                  "CSV (*.csv)")
            if path:
                import csv
                with open(path, 'w', newline='') as f:
                    w = csv.writer(f, delimiter=';')
                    w.writerow([f"S({name_A})", f"S({name_B})"])
                    for x, y in zip(path_x, path_y):
                        w.writerow([f"{x:.4f}", f"{y:.4f}"])
                    w.writerow([])
                    w.writerow(["Molecule", f"{mol_x:.4f}", f"{mol_y:.4f}"])


    def run_symmetry_calculation(self):
        if not self.selected_atoms:
            QMessageBox.warning(self, "No selection", "Select an atom")
            return

        shape_label = self.shape_results_window.get_selected_shape()
        if not shape_label:
            QMessageBox.warning(self, "No shape selected", 
                            "Select one shape in the table")
            return

        try:
            if len(self.selected_atoms) == 1:
                central_idx = self.selected_atoms[0]
                central_pos = self.atom_positions_3d[central_idx]
                central_sym = self.atom_types[central_idx]

                ligand_indices = self._ligands(central_idx)
                coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
                symbols = [central_sym]  + [self.atom_types[i] for i in ligand_indices]
            else:
                coords  = [self.atom_positions_3d[i] for i in self.selected_atoms]
                symbols = [self.atom_types[i] for i in self.selected_atoms]

            geo   = Geometry(symbols=symbols, positions=np.array(coords))
            cn_atoms = len(ligand_indices) if len(self.selected_atoms) == 1 else len(self.selected_atoms)
            shape_code = next(
                (code for code, name in SHAPE_REFERENCES.get(cn_atoms, {}).items() 
                if name == shape_label),
                None
                )
            symmetries = SHAPE_SYMMETRIES.get(shape_code, {})
            if not symmetries:
                QMessageBox.warning(self, "No symmetries", 
                    f"There is no symmetry for '{shape_label}'")
                return
            results = {}
            for sym_code, sym_desc in symmetries.items():
                value = geo.get_symmetry_measure(sym_code,central_atom=1)
                results[sym_code] = value
                print(f"{sym_desc}: {value:.4f}")
            
            msg = f"Symmetry measures for {shape_label}:\n\n"
            for sym_code, value in results.items():
                msg += f"  {sym_code}: {value:.4f}\n"
            QMessageBox.information(self, "Symmetry Measures", msg)

            self.export_data[f"Symmetry_{shape_label}"] = results

        except Exception as e:
            QMessageBox.warning(self, "Error", f"Error calculating symmetry: {e}")

    def update_shape_table(self, shape_results, atom_symbol, atom_idx, cn): 
        self.shape_results_window.update_results(shape_results, atom_symbol, atom_idx, cn)
    
    def toggle_polyhedron(self):
        if self._polyhedron_items:
            for item in self._polyhedron_items:
                self.molecule_view.removeItem(item)
            self._polyhedron_items = []
            self.shape_results_window.show_poly_button.setText("Show Polyhedra")
            self.shape_results_window.poly_color_button.setVisible(False)
            return

        shape_code = self.shape_results_window.get_selected_shape()
        if not shape_code or not self.selected_atoms:
            return

        if not hasattr(self,'poly_color'):
            self.poly_color = (0.1, 0.2, 0.6, 1.0)
        self.draw_polyhedron(shape_code)
        self.shape_results_window.show_poly_button.setText("Hide Polyhedra")
        self.shape_results_window.poly_color_button.setVisible(True)
        try:
            self.shape_results_window.poly_color_button.clicked.disconnect()
        except TypeError:
            pass
        self.shape_results_window.poly_color_button.clicked.connect(self.change_poly_color)

    def on_shape_selection_changed(self):
        if not self._polyhedron_items:
            return
        shape_code = self.shape_results_window.get_selected_shape()
        if not shape_code or not self.selected_atoms:
            return
        for item in self._polyhedron_items:
            self.molecule_view.removeItem(item)
        self._polyhedron_items = []
        self.draw_polyhedron(shape_code)

    def toggle_coord_env(self):
        if getattr(self, 'isolated_mode', False):
            self.show_all()
        else:
            self.isolate_coordination_environment()

    def open_selection_menu(self):
        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #ffffff;
                color: #0F141C;
            }
            QMenu::item { padding: 6px 20px; }
            QMenu::item:selected { background-color: #f0f4ff; color: #0F141C; }
        """)
                           
        act_rect = menu.addAction("Rectangular selection")
        act_free = menu.addAction("Free selection")
        pos = self.selection_button.mapToGlobal(self.selection_button.rect().bottomLeft())
        chosen = menu.exec_(pos)
        if chosen == act_rect:
            self.selection_mode = 'rect'
        elif chosen == act_free:
            self.selection_mode = 'free'






    def draw_polyhedron(self, shape_code):
        for item in self._polyhedron_items:
            self.molecule_view.removeItem(item)
        self._polyhedron_items = []

        if len(self.selected_atoms) == 1:
            central_idx = self.selected_atoms[0]
            central_pos = self.atom_positions_3d[central_idx]
            central_sym = self.atom_types[central_idx]
            ligand_indices = self._ligands(central_idx)
            coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
            symbols = [central_sym]  + [self.atom_types[i] for i in ligand_indices]
            central_param = 1
            has_center = True
        else:
            coords  = [self.atom_positions_3d[i] for i in self.selected_atoms]
            symbols = [self.atom_types[i] for i in self.selected_atoms]
            central_param = 0
            has_center = False

        ideal = get_ideal_coords(symbols, coords, shape_code, central_atom=central_param)
        if ideal is None or len(ideal) == 0:
            return

       
        ligand_pts = ideal[1:] if has_center else ideal
        n = len(ligand_pts)

        if n < 2:
            return
        edges = polyhedron_edges(ligand_pts)
        print('Vertex:', len(ligand_pts), "Aristas", len(edges))
        for i,j in polyhedron_edges(ligand_pts):
            line = GLLinePlotItem(
                pos=np.array([ligand_pts[i], ligand_pts[j]]),
                color=self.poly_color,
                width=6,
                antialias=True
            )
            self.molecule_view.addItem(line)
            self._polyhedron_items.append(line)


    def change_poly_color(self):
        initial = QColor(
            int(self.poly_color[0] * 255),
            int(self.poly_color[1] * 255),
            int(self.poly_color[2] * 255)
        )
        qcolor = QColorDialog.getColor(initial=initial, parent=self, title="Polyhedron Color")
        if not qcolor.isValid():
            return

        self.poly_color = (
            qcolor.red()   / 255.0,
            qcolor.green() / 255.0,
            qcolor.blue()  / 255.0,
            1.0
        )
        shape_code = self.shape_results_window.get_selected_shape()
        if shape_code:
            self.draw_polyhedron(shape_code)



        #----------------------Atom selection----------------------
    def run_point_group(self):
        if not self.atom_positions_3d or not self.atom_types:
            print('No molecule')
            return


        if len(self.selected_atoms) == 1:
            central = self.selected_atoms[0]
            indices = [central] + self._ligands(central)
        elif len(self.selected_atoms) > 1:
            indices = self.selected_atoms
        else:
            indices = list(range(len(self.atom_types)))



        symbols = [self.atom_types[i] for i in indices]
        coords  = [self.atom_positions_3d[i] for i in indices]

        point_group = calculate_point_group(coords, symbols)
        if point_group:
            self.current_point_group = str(point_group)
            print(f'Point group: {self.current_point_group}')
            self.pointgroup_table.setRowCount(0)
            self.pointgroup_table.setRowCount(1)
            item_label = QTableWidgetItem("Point Group")
            item_label.setFlags(item_label.flags() & ~Qt.ItemIsEditable)
            item_value = QTableWidgetItem(str(point_group))
            item_value.setFlags(item_value.flags() & ~Qt.ItemIsEditable)
            item_value.setTextAlignment(Qt.AlignCenter)
            self.pointgroup_table.setItem(0, 0, item_value)
            self.pointgroup_table.setVisible(True)
            self.close_pointgroup_button.setVisible(True)
            self.export_data["Point_Group"] = {"Point Group":str(point_group) }
        else:
            self.pointgroup_table.setVisible(False)


    def resaltar_atomo(self, idx):
        for i in range(len(self.spheres)):
            props = atom_color_dict.get(self.atom_types[i], {'color': 'gray', 'size': 0.3})
            if i in self.selected_atoms:
                self.spheres[i].setColor((0.969, 0.192, 0.027, 0.35)) 
                self.spheres[i].setGLOptions('translucent')
            else:
                self.spheres[i].setColor(colors.to_rgba(props['color']))
                self.spheres[i].setGLOptions('opaque')
        if self.selected_atoms and self.atom_info_window.isVisible():
            first_idx = self.selected_atoms[0]
            sym = self.atom_types[first_idx]
            pos = self.atom_positions_3d[first_idx]
            cn = len(self._ligands(first_idx))

            self.atom_info_window.update_atom(
                symbol=self.atom_types[first_idx],
                index=first_idx,
                x=pos[0],
                y=pos[1],
                z=pos[2],
                cn = cn,
            )

    def _ligands(self,central_idx):
        return get_ligands(
            central_idx, self.atom_positions_3d, self.atom_types,
            extra_bonds=self.manual_bonds, excluded_bonds=self.excluded_bonds
        )

    def set_labeled_atoms(self, indices):
        self._labeled_atoms = indices
        self.draw_labels()

    def open_label_window(self):
        if not self.atom_types:
            QMessageBox.warning(self, "No molecule", "Load a molecule first.")
            return
        self.label_window.populate()
        self.label_window.show()

    def open_style_window(self):
        if not self.atom_types:
            QMessageBox.warning(self, "No molecule", "Load a molecule first.")
            return
        self.style_window.populate()
        self.style_window.show()

    def open_style_window_bonds(self):
        if not self.atom_types:
            QMessageBox.warning(self, "No molecule", "Load a molecule first.")
            return
        self.style_window_bond.populate()
        self.style_window_bond.show()

    def draw_labels(self):
        from pyqtgraph.opengl import GLTextItem
        for item in self._label_items:
            self.molecule_view.removeItem(item)
        self._label_items = []

        for idx in self._labeled_atoms:
            if idx >= len(self.atom_positions_3d):
                continue
            pos = self.atom_positions_3d[idx]
            label = GLTextItem(
                pos=np.array([pos[0], pos[1], pos[2] + 0.4]),
                text=f"{self.atom_types[idx]}{idx}",
                color=(20, 20, 20, 255),
            )
            self.molecule_view.addItem(label)
            self._label_items.append(label)
    
    def delete_selection(self):
        if not self.selected_atoms:
            QMessageBox.information(self, "Delete", "No atoms selected.")
            return

        to_delete = set(self.selected_atoms)
        remaining = [
            (self.atom_types[i], *self.atom_positions_3d[i])
            for i in range(len(self.atom_types))
            if i not in to_delete
        ]

        self.selected_atoms = []
        self.selected_atom_index = None
        self.display_molecule(remaining)


    def set_bond_color(self, rgba):
        c = list(rgba)
        if len(c) == 3:
            c.append(1.0)
        self.bond_color = tuple(float(x) for x in c)
        self.display_molecule(
            [(self.atom_types[i], *self.atom_positions_3d[i])
             for i in range(len(self.atom_types))]
        )

    def set_bond_radius(self, radius):
        self.bond_radius = radius
        self.display_molecule([(self.atom_types[i], *self.atom_positions_3d[i])
             for i in range(len(self.atom_types))])


    def add_bond(self):
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(self, "Add bond", "Select exactly 2 atoms.")
            return
        
        i, j = sorted(self.selected_atoms[:2])

        if (i,j) in self.excluded_bonds:
            self.excluded_bonds.remove((i,j))
        

        if (i, j) in self._bond_pairs:
            QMessageBox.information(self, "Add bond", "That bond already exists.")
            return

        pos1 = self.atom_positions_3d[i]
        pos2 = self.atom_positions_3d[j]
        bond_color = getattr(self, 'bond_color', (0.3, 0.3, 0.3, 1))
        radius = getattr(self, 'bond_radius', 0.1)
        cylinder = create_cylinder(pos1, pos2, radius=radius, color=bond_color)
        self.molecule_view.addItem(cylinder)
        self.bonds.append(cylinder)
        self._bond_pairs.append((i, j))
        self.manual_bonds.append((i, j))
        print(f"Bond added: {self.atom_types[i]}{i} - {self.atom_types[j]}{j}")





    def delete_bond(self):
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(self, "Remove bond", "Select exactly 2 atoms.")
            return
        i, j = sorted(self.selected_atoms[:2])

        if (i, j) not in self._bond_pairs:
            QMessageBox.information(self, "Remove bond", "There is no bond between those atoms.")
            return
        if (i,j) in self._bond_pairs:
            idx = self._bond_pairs.index((i, j))
            self.molecule_view.removeItem(self.bonds[idx])
            del self.bonds[idx]
            del self._bond_pairs[idx]
        if (i, j) in self.manual_bonds:
            self.manual_bonds.remove((i, j))

        if (i,j) not in self.excluded_bonds:
            self.excluded_bonds.append((i,j))

        print(f"Bond removed: {self.atom_types[i]}{i} - {self.atom_types[j]}{j}")

    def restore_original(self):
        if not getattr(self, '_original_atoms', None):
            QMessageBox.information(self, 'Reset', 'Load a molecule first.')
            return
        self.selected_atoms = []
        self.display_molecule(self._original_atoms)
    def clear_symmetry_elements(self):
        if hasattr(self, '_sym_items'):
            for item in self._sym_items:
                self.molecule_view.removeItem(item)
        self._sym_items = []
        if hasattr(self, 'sym_label_overlay'):
            self.sym_label_overlay.hide()

    def reset_measurements(self):
        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
            self._measure_items = []
        self.sym_label_overlay.hide()
        self.btn_reset.setVisible(False)


    def clear_all_analysis(self):
        self.reset_measurements()
        self.clear_symmetry_elements()
        self.measurements_window.clear_all()
        for items in getattr(self, '_measure_registry', {}).values():
            for it in items:
                self.molecule_view.removeItem(it)
        self._measure_registry = {}
        self._measure_counter = 0
        for item in self._polyhedron_items:
            self.molecule_view.removeItem(item)
        self._polyhedron_items = []
        self.manual_bonds = []
        self.excluded_bonds = []
        self.selected_atoms = []
        self.selected_atom_index = None
        self.export_data = {}
        self.current_point_group = None
        self.isolated_mode = False
        self._hydrogens_hidden = False
        self.pointgroup_table.setVisible(False)
        for item in getattr(self, '_label_items', []):
            self.molecule_view.removeItem(item)
        self._label_items = []
        self._labeled_atoms = []
        self.close_pointgroup_button.setVisible(False)
        if self.shape_results_window.isVisible():
            self.shape_results_window.hide()
        if self.atom_info_window.isVisible():
            self.atom_info_window.hide()


    def draw_symmetry_operation(self, operation, center,label_text=None):
        matrix = np.real(operation.matrix_representation)
        op_id = (operation.label, tuple(np.round(matrix.flatten(), 3)))

        if getattr(self, '_current_sym_op', None) == op_id:
            self.clear_symmetry_elements()
            self.sym_label_overlay.hide()
            self._current_sym_op = None
            return

        self.clear_symmetry_elements()
        self._current_sym_op = op_id
        label = operation.label
        length = 3.0

        if label in ('E',):
            return

        elif label == 'i':
            md = MeshData.sphere(rows=20, cols=40, radius=0.95)
            sphere = GLMeshItem(meshdata=md, smooth=True, color=(0.969, 0.192, 0.027, 0.2), glOptions='translucent')
            sphere.translate(*center)
            self.molecule_view.addItem(sphere)
            self._sym_items.append(sphere)

        elif label.startswith('C') or label.startswith('S'):

            axis = get_symmetry_axis(matrix)
            axis = axis / np.linalg.norm(axis)
            p1 = center - axis * length
            p2 = center + axis * length
            line = GLLinePlotItem(
                pos=np.array([p1, p2]),
                color=(1, 0.3, 0, 1),
                width=4,
                antialias=True,
                glOptions='opaque' 
            )
            self.molecule_view.addItem(line)
            self._sym_items.append(line)
            arrow_length = 0.4
            arrow_radius = 0.12
            n_sides = 16
            tip = p2+axis*arrow_length
            base_center = p2
            
            arbitrary = np.array([1, 0, 0]) if abs(axis[0]) < 0.9 else np.array([0, 1, 0])
            u = np.cross(axis, arbitrary)
            u /= np.linalg.norm(u)
            v = np.cross(axis, u)

            angles = np.linspace(0, 2 * np.pi, n_sides, endpoint=False)
            base_verts = np.array([base_center + arrow_radius * (np.cos(a) * u + np.sin(a) * v)
                                for a in angles])

            all_verts = np.vstack([base_verts, tip])  
            tip_idx = n_sides
            faces = np.array([[i, (i + 1) % n_sides, tip_idx] for i in range(n_sides)])

            md = MeshData(vertexes=all_verts, faces=faces)
            cone = GLMeshItem(meshdata=md, smooth=True, color=(1, 0.3, 0, 1),
                          shader='shaded', glOptions='opaque')
            self.molecule_view.addItem(cone)
            self._sym_items.append(cone)
        elif label.startswith('s') or label.startswith('σ'):
            
            normal = get_plane_normal(matrix)
            normal = normal / np.linalg.norm(normal)
            self._draw_plane(center, normal, color=(0.2, 0.5, 1.0, 0.25))
            
        display = label_text if label_text else self.pretty_label(label)
        self.sym_label_overlay.setText(display)
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()

    def describe_orientation(self, operation):
        label = operation.label
        matrix = np.real(operation.matrix_representation())
        pretty = self.pretty_label(label)
        if label.startswith('C') or label.startswith('S'):
            v = get_symmetry_axis(matrix)
            v = v / np.linalg.norm(v)
            return f"{pretty} · Axis: ({v[0]:.2f}, {v[1]:.2f}, {v[2]:.2f})"
        elif label.startswith('s') or label.startswith('σ'):
            n = get_plane_normal(matrix)
            n = n / np.linalg.norm(n)
            return f"{pretty} · Plane normal: ({n[0]:.2f}, {n[1]:.2f}, {n[2]:.2f})"
        elif label == 'i':
            return "Inversion center"
        elif label == 'E':
            return "Identity"
        return ""


    def _draw_plane(self, center, normal, size=3.0, color=(0.2, 0.5, 1.0, 0.3)):
        arbitrary = np.array([1, 0, 0]) if abs(normal[0]) < 0.9 else np.array([0, 1, 0])
        v1 = np.cross(normal, arbitrary)
        v1 /= np.linalg.norm(v1)
        v2 = np.cross(normal, v1)

        corners = np.array([
            center + size * (v1+v2),
            center + size * (v1-v2),
            center + size * (-v1-v2),
            center + size * (-v1+v2)
        ])

        faces = np.array([[0, 1, 2], [0, 2, 3]])


        md = MeshData(vertexes=corners, faces=faces)
        mesh = GLMeshItem(meshdata=md, smooth=False, color=color, glOptions='translucent')
        self.molecule_view.addItem(mesh)
        self._sym_items.append(mesh)
    #---------------------------ATOM INFO ------------------------------#
    def show_atom_info(self):
        if len(self.selected_atoms) == 0:
            QMessageBox.warning(self, "No selection", "There are no atom selected")
            return
        if len(self.selected_atoms) > 1:
            QMessageBox.warning(self, "Multiple selection","More than one atom has been selected." )
            return
        
        idx = self.selected_atoms[0]
        sym = self.atom_types[idx]
        central_pos = self.atom_positions_3d[idx]

        cn = len(self._ligands(idx))


        self.atom_info_window.update_atom(
        symbol=sym,
        index=idx,
        x=central_pos[0],
        y=central_pos[1],
        z=central_pos[2],
        cn = cn,
    )


    def set_element_color(self, elem, rgb):
        rgba = (rgb[0], rgb[1], rgb[2], 1.0)
        for i, sym in enumerate(self.atom_types):
            if sym == elem and i not in self.selected_atoms:
                self.spheres[i].setColor(rgba)

        if elem in atom_color_dict:
            atom_color_dict[elem]['color'] = rgb
        else:
            atom_color_dict[elem] = {'color': rgb, 'size': 0.15}

    def get_element_size(self, elem):
        return atom_color_dict.get(elem, {}).get('size', 0.15)

    def set_element_size(self, elem, size):
        atom_color_dict.setdefault(elem, {'color': 'gray', 'size': 0.15})
        atom_color_dict[elem]['size'] = size

        self.display_molecule(
            [(self.atom_types[i], *self.atom_positions_3d[i])
             for i in range(len(self.atom_types))]
        )

    def set_background_color(self, hex_color):
        self.molecule_view.setBackgroundColor(hex_color)
    #---------------------------MEASURES------------------------------#
    def toggle_measurements(self):
        if not self.atom_positions_3d:
            QMessageBox.warning(self, "No molecule", "Load a molecule first.")
            return

        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;
                color: white;
                border: 1px solid black;
                font-weight: bold;
            }
            QMenu::item:selected { background-color: #555555; }
            QMenu::item { padding: 6px 20px; }
        """)

        action_dist     = menu.addAction("Distance (2 atoms)")
        action_angle    = menu.addAction("Angle (3 atoms)")
        action_dihedral = menu.addAction("Dihedral Angle (4 atoms)")

        pos = self.measure_button.mapToGlobal(self.measure_button.rect().bottomLeft())
        chosen = menu.exec_(pos)

        if chosen == action_dist:
            self.measure_distance()
        elif chosen == action_angle:
            self.measure_angle()
        elif chosen == action_dihedral:
            self.measure_dihedral()


    def _add_persistent_measure(self, points, text, color, label_pos=None):
        items = []
        for a, b in zip(points, points[1:]):
            dash = create_dashed_line(a, b, color=color)
            self.molecule_view.addItem(dash)
            items.append(dash)

        pos = label_pos if label_pos is not None else np.mean(points, axis=0)
        label = GLTextItem(pos=np.array(pos), text=text, color=(20, 20, 20, 255))
        self.molecule_view.addItem(label)
        items.append(label)

        for it in items:
            it.setVisible(False)   

        mid = self._measure_counter
        self._measure_counter += 1
        self._measure_registry[mid] = items
        return mid


    def set_measure_visible(self, mid, visible):
            for it in self._measure_registry.get(mid, []):
                it.setVisible(visible)

    def remove_measure(self, mid):
        for it in self._measure_registry.get(mid, []):
            self.molecule_view.removeItem(it)
        self._measure_registry.pop(mid, None)

    def measure_distance(self):
        if len(self.selected_atoms) != 2:
            QMessageBox.warning(self, "Selection", "Select exactly 2 atoms.")
            return

        idx1, idx2 = self.selected_atoms[0], self.selected_atoms[1]
        p1 = self.atom_positions_3d[idx1]
        p2 = self.atom_positions_3d[idx2]
        d  = calc_distance(p1, p2)
        s1 = f"{self.atom_types[idx1]}{idx1}"
        s2 = f"{self.atom_types[idx2]}{idx2}"

        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
        self._measure_items = []

        
        dash = create_dashed_line(p1, p2, color=(0.969, 0.192, 0.027, 0.65))
        self.molecule_view.addItem(dash)
        self._measure_items.append(dash)

        
        self.sym_label_overlay.setText(f"{s1} — {s2}:  {d:.4f} Å")
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()
        self.btn_reset.setVisible(True)
        mid = self._add_persistent_measure([p1, p2], f"{s1}–{s2}: {d:.3f} Å",
                                           color=(0.969, 0.192, 0.027, 0.65))
        self.measurements_window.add_distance(f"{s1} — {s2}", d, mid)
        self.measurements_window.show()

    def measure_angle(self):
        if len(self.selected_atoms) != 3:
            QMessageBox.warning(self, "Selection", "Select exactly 3 atoms.\nAngle is calculated at the middle atom.")
            return

        idx1, idx2, idx3 = self.selected_atoms[0], self.selected_atoms[1], self.selected_atoms[2]
        p1 = self.atom_positions_3d[idx1]
        p2 = self.atom_positions_3d[idx2] 
        p3 = self.atom_positions_3d[idx3]
        a  = calc_angle(p1, p2, p3)
        s1 = f"{self.atom_types[idx1]}{idx1}"
        s2 = f"{self.atom_types[idx2]}{idx2}"
        s3 = f"{self.atom_types[idx3]}{idx3}"

        
        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
        self._measure_items = []

        
        dash1 = create_dashed_line(p1, p2, color=(0.969, 0.192, 0.027, 0.65))
        dash2 = create_dashed_line(p3, p2, color=(0.969, 0.192, 0.027, 0.65))
        self.molecule_view.addItem(dash1)
        self.molecule_view.addItem(dash2)
        self._measure_items.extend([dash1, dash2])

        
        v1 = (p1 - p2) / np.linalg.norm(p1 - p2)
        v2 = (p3 - p2) / np.linalg.norm(p3 - p2)
        arc_radius = min(np.linalg.norm(p1 - p2), np.linalg.norm(p3 - p2)) * 0.3
        n_arc = 30
        arc_points = []
        for i in range(n_arc + 1):
            t = i / n_arc
            vec = (1 - t) * v1 + t * v2
            vec = vec / np.linalg.norm(vec)
            arc_points.append(p2 + arc_radius * vec)

        arc_line = GLLinePlotItem(
            pos=np.array(arc_points),
            color=(0.969, 0.192, 0.027, 0.65),
            width=5,
            antialias=True
        )
        self.molecule_view.addItem(arc_line)
        self._measure_items.append(arc_line)

        self.sym_label_overlay.setText(f"{s1} — {s2} — {s3}:  {a:.2f}°")
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()
        self.btn_reset.setVisible(True)
        mid = self._add_persistent_measure([p1, p2, p3], f"{s1} — {s2} — {s3}: {a:.2f}°",
                                           color=(0.969, 0.192, 0.027, 0.65))
        self.measurements_window.add_angle(f"{s1} — {s2} — {s3}", a, mid)
        self.measurements_window.show()

    def measure_dihedral(self):
        if len(self.selected_atoms) != 4:
            QMessageBox.warning(self, "Selection", "Select exactly 4 atoms.")
            return

        idx1, idx2, idx3, idx4 = self.selected_atoms[0], self.selected_atoms[1], self.selected_atoms[2], self.selected_atoms[3]
        p1 = self.atom_positions_3d[idx1]
        p2 = self.atom_positions_3d[idx2]
        p3 = self.atom_positions_3d[idx3]
        p4 = self.atom_positions_3d[idx4]
        d = calc_dihedral(p1, p2, p3, p4)

        s1 = f"{self.atom_types[idx1]}{idx1}"
        s2 = f"{self.atom_types[idx2]}{idx2}"
        s3 = f"{self.atom_types[idx3]}{idx3}"
        s4 = f"{self.atom_types[idx4]}{idx4}"

        if hasattr(self, '_measure_items'):
            for item in self._measure_items:
                self.molecule_view.removeItem(item)
        self._measure_items = []

        for a, b in [(p1, p2), (p2, p3), (p3, p4)]:
            dash = create_dashed_line(a, b, color=(0.122, 0.467, 0.957, 0.65))
            self.molecule_view.addItem(dash)
            self._measure_items.append(dash)

        axis = p3 - p2
        axis_norm = np.linalg.norm(axis)

        if axis_norm > 1e-6:
            axis_u = axis / axis_norm

            u = (p1 - p2)
            u = u - np.dot(u, axis_u) * axis_u
            u_norm = np.linalg.norm(u)

            v = (p4 - p3)
            v = v - np.dot(v, axis_u) * axis_u
            v_norm = np.linalg.norm(v)

            if u_norm > 1e-6 and v_norm > 1e-6:
                u = u / u_norm
                v = v / v_norm

                arc_radius = min(np.linalg.norm(p1 - p2), np.linalg.norm(p3 - p2)) * 0.35
                arc_center = p2

                n_arc = 40
                arc_points = []
                for i in range(n_arc + 1):
                    t = i / n_arc
                    vec = (1 - t) * u + t * v
                    vec_norm = np.linalg.norm(vec)
                    if vec_norm > 1e-6:
                        vec = vec / vec_norm
                    arc_points.append(arc_center + arc_radius * vec)

                arc_line = GLLinePlotItem(
                    pos=np.array(arc_points),
                    color=(0.122, 0.467, 0.957, 0.65),
                    width=5,
                    antialias=True
                )
                self.molecule_view.addItem(arc_line)
                self._measure_items.append(arc_line)

        self.sym_label_overlay.setText(f"{s1} — {s2} — {s3} — {s4}:  {d:.2f}°")
        self.sym_label_overlay.adjustSize()
        self.sym_label_overlay.show()
        self.btn_reset.setVisible(True)
        mid = self._add_persistent_measure([p1, p2, p3, p4], f"{s1} — {s2} — {s3} — {s4}: {d:.2f}°",
                                           color=(0.969, 0.192, 0.027, 0.65))
        self.measurements_window.add_dihedral(f"{s1} — {s2} — {s3} — {s4}", d, mid)
        self.measurements_window.show()


    #---------------------------ISOLATE------------------------------#
    def toggle_isolate(self):
        if not hasattr(self, 'isolated_mode'):
            self.isolated_mode = False


        if self.isolated_mode:
            self.show_all()
            return


        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;
                color: white;
                border: 1px solid black;
                font-weight: bold;
            }
            QMenu::item:selected {
                background-color: #555555;
            }
            QMenu::item {
                padding: 6px 20px;
            }
        """)

        action_coord = menu.addAction("Coordination Environment")
        action_sel   = menu.addAction("Selection")


        pos = self.isolate_button.mapToGlobal(
            self.isolate_button.rect().bottomLeft()
        )
        chosen = menu.exec_(pos)

        if chosen == action_coord:
            self.isolate_coordination_environment()
        elif chosen == action_sel:
            self.isolate_selection()

    def isolate_atoms(self, visible_indices):
        visible_set = set(visible_indices)
        for i, sphere in enumerate(self.spheres):
            sphere.setVisible(i in visible_set)
        for bond_data, bond_item in zip(self._bond_pairs, self.bonds):
            i, j = bond_data
            bond_item.setVisible(i in visible_set and j in visible_set)
        self.isolated_mode = True

    def isolate_selection(self):
        if not self.selected_atoms:
            print("No atoms selected")
            return
        self.isolate_atoms(self.selected_atoms)


    def isolate_coordination_environment(self):
        if not self.selected_atoms:
            print("No atoms selected")
            return

        central_idx = self.selected_atoms[0]
        central_pos = self.atom_positions_3d[central_idx]
        central_sym = self.atom_types[central_idx]

    
        ligand_indices = self._ligands(central_idx)

        visible_indices = [central_idx] + ligand_indices
        self.isolate_atoms(visible_indices)
        

    def show_all(self):
        for sphere in self.spheres:
            sphere.setVisible(True)
        for bond in self.bonds:
            bond.setVisible(True)
        self.isolated_mode = False
        

    def toggle_hydrogens(self):
        if not hasattr(self, 'spheres') or not self.spheres:
            return

        self._hydrogens_hidden = not getattr(self, '_hydrogens_hidden', False)

        for i, (sphere, atom_type) in enumerate(zip(self.spheres, self.atom_types)):
            if atom_type == 'H':
                sphere.setVisible(not self._hydrogens_hidden)

    
        for bond_data, bond_item in zip(self._bond_pairs, self.bonds):
            i, j = bond_data
            if self.atom_types[i] == 'H' or self.atom_types[j] == 'H':
                bond_item.setVisible(not self._hydrogens_hidden)

        self.hide_h_button.setToolTip(
            "Show Hydrogens" if self._hydrogens_hidden else "Hide Hydrogens"
        )


    #-----------------------XYZ AXIS------------------------ 

    def crear_ejes_xyz(self):
        eje_x = GLLinePlotItem(pos=np.array([[3, -3, -3], [4, -3, -3]]), color=(1, 0, 0, 1), width=3, antialias=True)
        eje_y = GLLinePlotItem(pos=np.array([[3, -3, -3], [3, -2, -3]]), color=(0, 1, 0, 1), width=3, antialias=True)
        eje_z = GLLinePlotItem(pos=np.array([[3, -3, -3], [3, -3, -2]]), color=(0, 0, 1, 1), width=3, antialias=True)
        return eje_x, eje_y, eje_z

    def toggle_ejes_xyz(self):
        if self.ejes_visibles:
            for eje in [self.eje_x, self.eje_y, self.eje_z]:
                if eje is not None and eje in self.molecule_view.items:
                    self.molecule_view.removeItem(eje)
            self.eje_x = self.eje_y = self.eje_z = None
            self.ejes_visibles = False
        else:
            self.eje_x, self.eje_y, self.eje_z = self.crear_ejes_xyz()
            for eje in [self.eje_x, self.eje_y, self.eje_z]:
                self.molecule_view.addItem(eje)
            self.ejes_visibles = True

    #----------------------Load and visualization--------------------
    def display_molecule(self, atoms):
        if hasattr(self, 'spheres'):
            for s in self.spheres:
                self.molecule_view.removeItem(s)
        if hasattr(self, 'bonds'):
            for b in self.bonds:
                self.molecule_view.removeItem(b)

        self.spheres = []
        self.bonds = []
        self._bond_pairs = []   
        self.atom_positions_3d = []
        self.atom_types = []
        self.selected_atom_index = None

        for atom_type, x, y, z in atoms:
            props = atom_color_dict.get(atom_type, {'color': 'gray', 'size': 0.15})
            color = colors.to_rgba(props['color'])
            sphere = create_sphere(radius=props['size']*2, color=color)
            sphere.translate(x, y, z)
            self.molecule_view.addItem(sphere)
            self.spheres.append(sphere)
            self.atom_positions_3d.append(np.array([x, y, z]))
            self.atom_types.append(atom_type)

        bond_count = 0
        for i, (type1, x1, y1, z1) in enumerate(atoms):
            for j, (type2, x2, y2, z2) in enumerate(atoms):
                if i < j:
                    pos1 = np.array([x1, y1, z1])
                    pos2 = np.array([x2, y2, z2])
                    dist = np.linalg.norm(pos1 - pos2)
                    threshold = covalent_threshold(type1, type2)
                    
                    if dist < threshold:
                        bond_count += 1
                        bond_color = getattr(self, 'bond_color', (0.3, 0.3, 0.3, 1))
                        print("bond_color =", repr(bond_color), "tipos:", [type(x) for x in bond_color])
                        radius = getattr(self, 'bond_radius', 0.1)
                        cylinder = create_cylinder(pos1, pos2, radius=radius, color=bond_color)
                        self._bond_pairs.append((i, j))  
                        self.bonds.append(cylinder)      
                        self.molecule_view.addItem(cylinder)
        if self.atom_positions_3d:
            centroid = np.mean(self.atom_positions_3d,axis=0)
            self.molecule_view.opts['center'] = pg.Vector(centroid[0], centroid[1], centroid[2])
            max_extent = max(np.linalg.norm(pos-centroid) for pos in self.atom_positions_3d)
            self.molecule_view.opts['distance'] = max(max_extent * 3,10)

        

    def load_file(self):
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getOpenFileName(
            self, "Open Molecule File", "",
            "Structure Files (*.xyz *.pdb *.cif);;All Files (*)",
            options=options
        )
        if not file_path:
            return
        self.load_path(file_path)

    def load_path(self, file_path):
        try:
            atoms,cell = load_structure(file_path)
        except Exception as e:
            QMessageBox.warning(self, "Load error", f"Could not read file:\n{e}")
            return
        self.cell = cell
        self._original_atoms = list(atoms)
        self.clear_all_analysis()
        self.display_molecule(atoms)
        self.btn_reset.setVisible(False)


    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls():
            for url in event.mimeData().urls():
                path = url.toLocalFile()
                if path.lower().endswith(('.xyz', '.pdb', '.cif')):
                    event.acceptProposedAction()
                    return
        event.ignore()

    def dropEvent(self, event):
        for url in event.mimeData().urls():
            path = url.toLocalFile()
            if path.lower().endswith(('.xyz', '.pdb', '.cif')):
                self.load_path(path)
                break

    def close_application(self):
        self.close()

    def export_image(self, fmt="png"):
        from PyQt5.QtWidgets import QComboBox

        fmt_labels = {
            "png":  "PNG Files (*.png)",
            "tiff": "TIFF Files (*.tiff)",
            "jpg":  "JPEG Files (*.jpg)",
        }

        dialog = QDialog(self)
        dialog.setWindowTitle("Export Image")
        dialog.setFixedSize(300, 200)
        layout = QVBoxLayout(dialog)

        combo_style = """
            QComboBox {padding: 4px; font-size: 13px;}
            QComboBox QAbstractItemView {
                selection-background-color: #e0e7ff;
                selection-color: #1f2937;
            }

        """

        layout.addWidget(QLabel("Resolution:"))
        combo = QComboBox()
        combo.addItems([
            "1x  (screen size)",
            "2x  (publication)",
            "4x  (high resolution)",
        ])
        combo.setCurrentIndex(1)
        combo.setStyleSheet(combo_style)
        layout.addWidget(combo)

        layout.addWidget(QLabel("Background:"))
        bg_combo = QComboBox()
        bg_combo.addItems(["Current ", "White", "Black"])
        bg_combo.setStyleSheet(combo_style)
        layout.addWidget(bg_combo)

        btn_row = QHBoxLayout()
        ok_btn = QPushButton("Export")
        ok_btn.setStyleSheet("background-color: white; color: black; font-weight: bold;")
        cancel_btn = QPushButton("Cancel")
        cancel_btn.setStyleSheet(ok_btn.styleSheet())
        ok_btn.clicked.connect(dialog.accept)
        cancel_btn.clicked.connect(dialog.reject)
        btn_row.addWidget(ok_btn)
        btn_row.addWidget(cancel_btn)
        layout.addLayout(btn_row)

        if dialog.exec_() != QDialog.Accepted:
            return

        factor = [1, 2, 4][combo.currentIndex()]
        bg_choice = bg_combo.currentIndex()

        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save Image", f"molecule.{fmt}",
            f"{fmt_labels[fmt]};;All Files (*)"
        )
        if not file_path:
            return

        if bg_choice == 1:
            self.molecule_view.setBackgroundColor('w')
        elif bg_choice == 2:
            self.molecule_view.setBackgroundColor('k')

        original_size = self.molecule_view.size()
        w = original_size.width()  * factor
        h = original_size.height() * factor
        self.molecule_view.resize(w, h)
        self.molecule_view.repaint()


        img = self.molecule_view.grabFramebuffer()

        self.molecule_view.resize(original_size)
        if bg_choice in (1, 2):
            self.molecule_view.setBackgroundColor('#f8f9fa')
        self.molecule_view.repaint()

        quality = 95 if fmt == "jpg" else -1
        if img.save(file_path, fmt.upper(), quality):
            print(f"Image saved: {file_path}  ({img.width()}x{img.height()} px)")
        else:
            QMessageBox.warning(self, "Error", f"Could not save image to:\n{file_path}")


    def export_f(self):
        if not self.export_data:
            QMessageBox.information(self, "Export", "No results to export yet.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results", "results.csv", "CSV Files (*.csv)"
        )
        if not file_path:
            return
        try:
            with open(file_path, "w", newline="", encoding="utf-8") as fh:
                writer = csv.writer(fh, delimiter=';')
                for section, entries in self.export_data.items():
                    writer.writerow([section])
                    for prop, value in entries.items():
                        writer.writerow([prop, value])
                    writer.writerow([])
            print(f"Results exported: {file_path}")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not export results:\n{e}")

    def toggle_export(self):
        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu {
                background-color: #2b2b2b;
                color: white;
                border: 1px solid black;
                font-weight: bold;
            }
            QMenu::item:selected { background-color: #555555; }
            QMenu::item { padding: 10px 20px; }
        """)

        action_results = menu.addAction("Export Results (.csv)")
        action_geometry = menu.addAction("Save geometry (.xyz)")
        menu.addSeparator()
        action_png  = menu.addAction("Export Image as PNG")
        action_tiff = menu.addAction("Export Image as TIFF")
        action_jpg  = menu.addAction("Export Image as JPEG")

        pos = self.export_button.mapToGlobal(self.export_button.rect().bottomLeft())
        chosen = menu.exec_(pos)

        if chosen == action_results:
            self.export_f()
        elif chosen == action_geometry:
            self.save_geometry()
        elif chosen == action_png:
            self.export_image(fmt="png")
        elif chosen == action_tiff:
            self.export_image(fmt="tiff")
        elif chosen == action_jpg:
            self.export_image(fmt="jpg")

    def save_geometry(self):
        if not self.atom_positions_3d:
            QMessageBox.information(self, "Save geometry", "No molecule loaded.")
            return

        visible = [
            i for i in range(len(self.atom_types))
            if self.spheres[i].visible()
        ]
        if not visible:
            QMessageBox.information(self, "Save geometry", "No visible atoms to save.")
            return
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Save geometry", "fragment.xyz", "XYZ Files (*.xyz)"
        )
        if not file_path:
            return

        try:
            with open(file_path, "w") as f:
                f.write(f"{len(visible)}\n")
                f.write("Exported from Molektra\n")
                for i in visible:
                    sym = self.atom_types[i]
                    x, y, z = self.atom_positions_3d[i]
                    f.write(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")
            print(f"Geometry saved: {file_path} ({len(visible)} atoms)")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Could not save geometry:\n{e}")


    def run_symmetry_posym(self):
        if not hasattr(self, 'current_point_group') or not self.current_point_group:
            QMessageBox.warning(self, "Missing step", "Ep! First get the point group!")
            return
        try:
            if len(self.selected_atoms) == 1:
                central_idx = self.selected_atoms[0]
                central_pos = self.atom_positions_3d[central_idx]
                central_sym = self.atom_types[central_idx]
                ligand_indices = self._ligands(central_idx)
                coords  = [central_pos] + [self.atom_positions_3d[i] for i in ligand_indices]
                symbols = [central_sym]  + [self.atom_types[i] for i in ligand_indices]
                central_label = f"{central_sym}{central_idx}"
            elif len(self.selected_atoms) > 1:
                coords = [self.atom_positions_3d[i] for i in self.selected_atoms]
                symbols = [self.atom_types[i] for i in self.selected_atoms]
                central_label = 'selection'
            else:
                coords  = [self.atom_positions_3d[i] for i in range(len(self.atom_types))]
                symbols = [self.atom_types[i] for i in range(len(self.atom_types))]
                central_label = 'whole_molecule'

            self._last_sym_coords = coords
            self._last_sym_symbols = symbols


            if not self.current_point_group:
                QMessageBox.warning(self, "No point group", "Could not determine the molecule's point group.")
                return

            results, operations, center = calculate_symmetry_posym(
                coords, symbols, groups=[self.current_point_group])

            if not results:
                QMessageBox.information(self, "Error",'No compatible point group found')
                return

            self.show_posym_results(results,operations,center)


            self.export_data[f"CSM_{central_label}"] = results

            if operations:
                ops_data = {}
                counter = {}
                for op in operations:
                    counter[op.label] = counter.get(op.label, 0) + 1
                    n_total = sum(1 for o in operations if o.label == op.label)
                    name = f"{op.label}_{counter[op.label]}" if n_total > 1 else op.label
                    matrix = np.round(np.real(op.matrix_representation), 4)
                    ops_data[name] = " | ".join(
                        "[" + ", ".join(f"{v:.4f}" for v in row) + "]" for row in matrix
                    )
                self.export_data[f"SymmetryOps_{central_label}"] = ops_data

        except Exception as e:
            QMessageBox.warning(self, "Error posym", f"Error calculating CSM: {e}")

    def pretty_label(self, label):
        subscripts   = str.maketrans('0123456789', '₀₁₂₃₄₅₆₇₈₉')
        superscripts = str.maketrans('0123456789', '⁰¹²³⁴⁵⁶⁷⁸⁹')

        if label in ('sh', 'sv', 'sd'):
            return 'σ' + label[1:]
        if label.startswith('s') and not label.startswith('S'):
            return 'σ' + label[1:]

        if '^' in label and '_' in label:
            base = label.split('^')[0]                  
            power = label.split('^')[1].split('_')[0]    
            order = label.split('_')[-1]                  
            return f"{base}{order.translate(subscripts)}{power.translate(superscripts)}"

        if label and label[0] in ('C', 'S'):
            rest = label[1:]
            i = 0
            while i < len(rest) and rest[i].isdigit():
                i += 1
            digits = rest[:i]
            suffix = rest[i:]         
            if digits:
                return label[0] + digits.translate(subscripts) + suffix

        return label  



    def add_group_section(self, results, operations):
        for group, value in results.items():
            box = QGroupBox(f"{group}")
            box.setStyleSheet("""
                QGroupBox { font-weight: 600; color: #374151;
                    border: 1px solid #e5e7eb; border-radius: 10px;
                    margin-top: 10px; padding: 12px 10px 10px 10px;
                    background-color: #ffffff; }
                QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 4px; }
            """)
            v = QVBoxLayout(box)

            lbl = QLabel(f"CShM (S) = {value:.4f}")
            lbl.setStyleSheet("font-size: 13px; font-weight: bold;")
            v.addWidget(lbl)

            grouped = {}
            for op in operations:
                if op.label == 'E':
                    continue
                grouped.setdefault(op.label, []).append(op)

            button_row = QHBoxLayout()
            for label, ops in grouped.items():
                if len(ops) == 1:

                    btn = QPushButton(self.pretty_label(label))
                    btn.setFixedHeight(32)
                    btn.clicked.connect(
                        lambda _, o=ops[0]: self.draw_symmetry_operation(
                            o, self._sym_center, label_text=self.describe_orientation(o)))
                    button_row.addWidget(btn)
                else:

                    btn = QPushButton(f"{self.pretty_label(label)} ({len(ops)})")
                    btn.setFixedHeight(32)
                    btn.clicked.connect(
                        lambda _, lb=label, os=ops, b=btn: self._show_ops_menu(lb, os, b))
                    button_row.addWidget(btn)
            v.addLayout(button_row)

            self.posym_layout.addWidget(box)

    def _show_ops_menu(self, label, ops, button):
        menu = QMenu(self)
        menu.setStyleSheet("""
            QMenu { background-color: #ffffff; color: #1f2937;
                    border: 1px solid #e5e7eb; border-radius: 6px; }
            QMenu::item { padding: 6px 20px; }
            QMenu::item:selected { background-color: #f0f4ff; color: #1f2937; }
        """)
        actions = {}
        for i, op in enumerate(ops, 1):
            v = self.describe_orientation(op)
            act = menu.addAction(f"{label} #{i}  —  {v}")
            actions[act] = op

        pos = button.mapToGlobal(button.rect().bottomLeft())
        chosen = menu.exec_(pos)
        if chosen in actions:
            op = actions[chosen]
            self.draw_symmetry_operation(op, self._sym_center,
                                         label_text=self.describe_orientation(op))


    def add_comparison_group(self):
        group = self.compare_pg_input.text().strip()
        if not group or not hasattr(self, '_last_sym_coords'):
            return
        try:
            results, operations, center = calculate_symmetry_posym(
                self._last_sym_coords, self._last_sym_symbols, groups=[group])
        except Exception as e:
            QMessageBox.warning(self, "Compare", f"Could not measure {group}: {e}")
            return
        if not results:
            QMessageBox.information(self, "Compare",
                f"'{group}' couldn't be measured (invalid or incompatible).")
            return
        self.add_group_section(results, operations)
        self.export_data.setdefault("CSM_comparison", {}).update(results)
        self.compare_pg_input.clear()



    def describe_orientation(self, operation):
        label = operation.label
        matrix = np.real(operation.matrix_representation)
        if label.startswith('C') or label.startswith('S'):
            v = get_symmetry_axis(matrix); v = v/np.linalg.norm(v)
            return f"{label}  ·  axis ({v[0]:.2f}, {v[1]:.2f}, {v[2]:.2f})"
        elif label.startswith('s') or label.startswith('σ'):
            n = get_plane_normal(matrix); n = n/np.linalg.norm(n)
            return f"{label}  ·  plane ⊥ ({n[0]:.2f}, {n[1]:.2f}, {n[2]:.2f})"
        elif label == 'i':
            return "inversion center"
        return label

    def show_posym_results(self, results, operations,center):
        self.posym_window = QDialog(self)
        self.posym_window.setWindowTitle("Symmetry Measures")
        self.posym_window.setStyleSheet("""
            QDialog { background-color: #f8f9fa; }
            QLabel { color: #1f2937; }
            QPushButton {
                background-color: #ffffff; color: #1f2937; font-weight: 600;
                border: 1px solid #e5e7eb; border-radius: 8px; padding: 6px 10px;
            }
            QPushButton:hover { background-color: #f3f4f6; }
            QLineEdit {
                background-color: #ffffff; color: #1f2937;
                border: 1px solid #e5e7eb; border-radius: 6px; padding: 6px;
            }
        """)
        self.posym_layout = QVBoxLayout(self.posym_window)
        self._sym_center = center

        self.add_group_section(results, operations)
        compare_row = QHBoxLayout()
        self.compare_pg_input = QLineEdit()
        self.compare_pg_input.setPlaceholderText("Compare with another group (e.g. Oh)")
        compare_btn = QPushButton("Compare")
        compare_btn.setFixedSize(110, 60)
        compare_btn.clicked.connect(self.add_comparison_group)
        compare_row.addWidget(self.compare_pg_input)
        compare_row.addWidget(compare_btn)
        self.posym_layout.addLayout(compare_row)

        self.posym_window.show()


    def show_info(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("About Us")
        dialog.setFixedSize(420, 320)
        dialog.setStyleSheet("background-color: #1e1e2e;")
        layout = QVBoxLayout(dialog)
        layout.setSpacing(10)
        layout.setContentsMargins(24, 20, 24, 20)
        #-------LOGO------
        logo_label = QLabel()
        logo_path = os.path.join(os.path.dirname(__file__), "molektra.png")  
        if os.path.exists(logo_path):
            pix = QPixmap(logo_path).scaled(110, 110, Qt.KeepAspectRatio, Qt.SmoothTransformation)
            logo_label.setPixmap(pix)
        else:
            logo_label.setText("X")
            logo_label.setStyleSheet("font-size: 48px;")
        logo_label.setAlignment(Qt.AlignCenter)
        layout.addWidget(logo_label)

        #---NAME---
        title = QLabel("Molektra")
        title.setAlignment(Qt.AlignCenter)
        title.setStyleSheet("color: white; font-size: 16px; font-weight: bold;")
        layout.addWidget(title)

        subtitle = QLabel("A visual tool for shape and symmetry analysis of molecules")
        subtitle.setAlignment(Qt.AlignCenter)
        subtitle.setWordWrap(True)
        subtitle.setStyleSheet("color: #9ca3af; font-size: 11px;")
        layout.addWidget(subtitle)

        sep = QLabel()
        sep.setFixedHeight(1)
        sep.setStyleSheet("background-color: #374151;")
        layout.addWidget(sep)

    #------------LINKS------------------#
        links = [
            ("GitHub",  "https://github.com/LoraSM/Molektra",  "#6366f1"),
            ("Zenodo",  "https://zenodo.org/",        "#0ea5e9"),
            ("Paper",   "https://doi.org/",                 "#10b981"),
        ]

        btn_row = QHBoxLayout()
        btn_row.setSpacing(10)
        for label, url, color in links:
            btn = QPushButton(label)
            btn.setFixedHeight(34)
            btn.setStyleSheet(f"""
                QPushButton {{
                    background-color: {color};
                    color: white;
                    font-weight: bold;
                    border-radius: 6px;
                    border: none;
                    font-size: 12px;
                }}
                QPushButton:hover {{ opacity: 0.85; }}
            """)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(lambda _, u=url: QDesktopServices.openUrl(QUrl(u)))
            btn_row.addWidget(btn)
        layout.addLayout(btn_row)

        layout.addStretch()

        #About us
        credits = QLabel("© 2026 IQTC - Universitat de Barcelona · Version 1.0")
        credits.setAlignment(Qt.AlignCenter)
        credits.setStyleSheet("color: #6b7280; font-size: 10px;")
        layout.addWidget(credits)

        close_btn = QPushButton("Close")
        close_btn.setFixedHeight(32)
        close_btn.setStyleSheet("""
            QPushButton {
                background-color: #374151; color: white;
                font-weight: bold; border-radius: 4px;
                border: 2px solid #4b5563;
            }
            QPushButton:hover { background-color: #4b5563; }
        """)
        close_btn.clicked.connect(dialog.close)
        layout.addWidget(close_btn)

        dialog.exec_()